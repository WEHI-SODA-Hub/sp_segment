#!/usr/bin/env python3
"""
Score the membrane combination comparison.

Reads the final cell masks (cellmeasurement/<sample>_mask.tiff) written by the
pipeline for every method x segmenter sample, and reports, per sample:

  cells            number of cells
  median_um2       median cell area in um^2
  nuclei_in_cell   % of DAPI nuclei whose centroid falls inside a cell
  <marker>_px      % of the marker's positive pixels that fall inside a cell
  <marker>_cells   number of cells whose mean rescaled marker value is >= 0.3

A pixel is marker-positive when the marker, rescaled with the same rule as
make_inputs.py (median -> 0, 99.5th percentile -> 1), is >= 0.3 and fewer than
3 of the 5 markers are >= 0.3 there. Pixels bright in 3 or more markers are
treated as autofluorescence (common in brain tissue) and excluded. Nuclei are
blobs found by Laplacian-of-Gaussian detection on the DAPI channel.

Also writes contrast.csv: for each marker, the pixels where only that marker is
bright (it is above its own 99th percentile and every other marker is below its
90th), and the difference between their mean and the mean of that marker's dim
pixels (below its median), divided by the standard deviation of the dim
pixels. It is computed on each membrane input and on the marker alone. The
value does not change under a linear rescale, so it is the same after each
segmenter's own percentile normalisation.

Also writes an overlay figure of a 320 x 320 pixel region.

Usage:
  python analyse.py --tiff <background-subtracted.ome.tif> --results <outdir> \
      --inputs <make_inputs outdir> --out <dir>
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage.feature import blob_log
from skimage.segmentation import find_boundaries

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
from tiff_channels import get_channel_names, read_channels  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_inputs import MARKERS, NUCLEAR, SEGMENTERS, rescale  # noqa: E402

METHODS = ["raw_max", "norm_max", "norm_mean", "single"]
POSITIVE = 0.3
PIXEL_UM2 = 0.28 ** 2
REGION = 320


def short(marker):
    return marker.replace("VCW_", "").replace("_diluted", "").replace("_AIF1", "")


def nuclei_centroids(dapi):
    """Nucleus centres: LoG blobs of 2-6 um radius (sigma 5-15 px at 0.28 um/px)."""
    image = rescale(dapi).clip(0, 1)
    blobs = blob_log(image, min_sigma=5, max_sigma=15, num_sigma=6, threshold=0.05)
    return blobs[:, :2].astype(int)


def score(mask, centroids, positives, rescaled):
    inside = mask > 0
    ids = np.unique(mask)
    ids = ids[ids > 0]
    row = {
        "cells": len(ids),
        "median_um2": round(float(np.median(np.bincount(mask.ravel())[ids])) * PIXEL_UM2, 1),
        "nuclei_in_cell": round(100 * float(inside[centroids[:, 0], centroids[:, 1]].mean()), 1),
    }
    for marker, pos, values in zip(MARKERS, positives, rescaled):
        row[f"{short(marker)}_px"] = round(100 * float(inside[pos].mean()), 1)
        means = ndi.mean(values, labels=mask, index=ids)
        row[f"{short(marker)}_cells"] = int((np.asarray(means) >= POSITIVE).sum())
    return row


def contrast_table(raw, inputs):
    rows = []
    for t, marker in enumerate(MARKERS):
        x = raw[t]
        others = [i for i in range(len(MARKERS)) if i != t]
        others_dim = np.all([raw[i] < np.percentile(raw[i], 90) for i in others], axis=0)
        pos = others_dim & (x > np.percentile(x, 99))
        neg = others_dim & (x < np.percentile(x, 50))

        def contrast(img):
            img = img.astype(np.float64)
            return round(float((img[pos].mean() - img[neg].mean()) / img[neg].std()), 2)

        rows.append({"marker": short(marker), "pixels": int(pos.sum()),
                     **{m: contrast(inputs[m]) for m in METHODS}, "marker_alone": contrast(x)})
    return rows


def pick_region(positive):
    """Top-left corner of the REGION window with most positive pixels."""
    density = ndi.uniform_filter(positive.astype(np.float64), size=REGION, mode="constant")
    half = REGION // 2
    inner = density[half:-half, half:-half]
    y, x = np.unravel_index(np.argmax(inner), inner.shape)
    return y, x


def figure(path, dapi, rescaled, inputs, masks, corner):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    y, x = corner
    win = (slice(y, y + REGION), slice(x, x + REGION))

    def unit(a):
        lo, hi = np.percentile(a, [1, 99.5])
        return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)

    # Same background in every segmentation panel: CD45 red, GFAP green, DAPI blue.
    background = np.dstack([
        np.clip(rescaled[MARKERS.index("VCW_CD45")][win], 0, 1),
        np.clip(rescaled[MARKERS.index("VCW_GFAP_diluted")][win], 0, 1),
        unit(dapi[win]) * 0.7,
    ])

    fig, axes = plt.subplots(len(METHODS), 1 + len(SEGMENTERS),
                             figsize=(3 * (1 + len(SEGMENTERS)), 3 * len(METHODS)))
    for r, method in enumerate(METHODS):
        axes[r, 0].imshow(unit(inputs[method][win]), cmap="gray")
        axes[r, 0].set_ylabel(method, fontsize=11)
        for c, seg in enumerate(SEGMENTERS, start=1):
            panel = background.copy()
            panel[find_boundaries(masks[(method, seg)][win], mode="inner")] = (1, 1, 1)
            axes[r, c].imshow(panel)
    for c, title in enumerate(["membrane input"] + SEGMENTERS):
        axes[0, c].set_title(title, fontsize=11)
    for ax in axes.ravel():
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiff", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    names = get_channel_names(args.tiff)
    dapi = read_channels(args.tiff, [names.index(NUCLEAR)])[0]
    raw = read_channels(args.tiff, [names.index(m) for m in MARKERS])
    rescaled = np.stack([rescale(p) for p in raw])
    bright = rescaled >= POSITIVE
    autofluorescence = bright.sum(axis=0) >= 3
    positives = bright & ~autofluorescence
    centroids = nuclei_centroids(dapi)
    print(f"{len(centroids)} DAPI nuclei; {int(autofluorescence.sum())} autofluorescent pixels; "
          "positive pixels per marker: "
          + ", ".join(f"{short(m)}={int(p.sum())}" for m, p in zip(MARKERS, positives)))

    masks, rows = {}, []
    for method in METHODS:
        for seg in SEGMENTERS:
            mask = tifffile.imread(Path(args.results) / "cellmeasurement" / f"{method}_{seg}_mask.tiff")
            masks[(method, seg)] = mask
            rows.append({"method": method, "segmenter": seg,
                         **score(mask, centroids, positives, rescaled)})

    with open(out / "scores.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out / 'scores.csv'}")
    for row in rows:
        print(row)

    inputs = {m: tifffile.imread(Path(args.inputs) / f"{m}.ome.tif")[1] for m in METHODS}
    contrast = contrast_table(raw, inputs)
    with open(out / "contrast.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(contrast[0]))
        writer.writeheader()
        writer.writerows(contrast)
    print(f"wrote {out / 'contrast.csv'}")
    for row in contrast:
        print(row)
    corner = pick_region(positives[MARKERS.index("VCW_CD45")])
    figure(out / "overlays.png", dapi, rescaled, inputs, masks, corner)


if __name__ == "__main__":
    main()
