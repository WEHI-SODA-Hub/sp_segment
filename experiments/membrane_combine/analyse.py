#!/usr/bin/env python3
"""
Score the membrane combination comparison.

Reads the final cell masks (cellmeasurement/<region>_<method>_<segmenter>_mask.tiff)
written by the pipeline and the <region>_source.ome.tif files written by
make_inputs.py, and writes:

scores.csv, one row per region, method and segmenter:
  cells           number of cells
  median_um2      median cell area in um^2
  nuclei          number of DAPI nuclei (Laplacian-of-Gaussian blobs, radius
                  about 2-6 um)
  nuclei_in_cell  % of those nuclei whose centre is inside a cell
  gt_cells        number of ground-truth cells in the region (if configured)
  gt_in_cell      % of ground-truth centroids inside a cell
  gt_alone        % of ground-truth centroids inside a cell that contains no
                  other ground-truth centroid

gt_by_type.csv: gt_in_cell and gt_alone per ground-truth label, pooled over
regions, per method and segmenter.

contrast.csv, per region and marker: the pixels where only that marker is
bright (above its own 99th percentile, every other marker below its 90th),
and the difference between their mean and the mean of that marker's dim pixels
(below its median), divided by the standard deviation of the dim pixels. It is
computed on each membrane input and on the marker alone, and does not change
under a linear rescale.

overlays_<region>.png: cell outlines for every method and segmenter on a
400 x 400 pixel window.

Usage:
  python analyse.py --config <config.json> --inputs <make_inputs outdir> \
      --results <pipeline outdir> --out <dir>
"""

import argparse
import csv
import io
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage.feature import blob_log
from skimage.segmentation import find_boundaries

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
from combine_channels import source_physical_size  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_inputs import METHODS, rescale  # noqa: E402

WINDOW = 400


def nuclei_centres(dapi, pixel_um):
    """Nucleus centres: LoG blobs of radius about 2-6 um (radius = sigma * sqrt(2))."""
    image = rescale(dapi).clip(0, 1)
    blobs = blob_log(image, min_sigma=2.0 / pixel_um / np.sqrt(2),
                     max_sigma=6.0 / pixel_um / np.sqrt(2), num_sigma=6, threshold=0.05)
    return blobs[:, :2].astype(int)


def load_ground_truth(path, image, y0, x0, size):
    """Ground-truth (row, col, label) inside the region, in region coordinates."""
    lines = [line for line in open(path) if not line.startswith("#")]
    out = []
    for row in csv.DictReader(io.StringIO("".join(lines))):
        if row["Image"] != image:
            continue
        r, c = float(row["CentroidY"]) - y0, float(row["CentroidX"]) - x0
        if 0 <= r < size and 0 <= c < size:
            out.append((int(r), int(c), row["Label"]))
    return out


def score(mask, pixel_um2, nuclei, gt):
    ids = np.unique(mask)
    ids = ids[ids > 0]
    inside = mask > 0
    row = {
        "cells": len(ids),
        "median_um2": round(float(np.median(np.bincount(mask.ravel())[ids])) * pixel_um2, 1)
        if len(ids) else 0.0,
        "nuclei": len(nuclei),
        "nuclei_in_cell": round(100 * float(inside[nuclei[:, 0], nuclei[:, 1]].mean()), 1),
    }
    hits = []
    if gt:
        labels_at = np.array([mask[r, c] for r, c, _ in gt])
        counts = defaultdict(int)
        for lab in labels_at[labels_at > 0]:
            counts[lab] += 1
        hits = [(lab > 0, lab > 0 and counts[lab] == 1) for lab in labels_at]
        row["gt_cells"] = len(gt)
        row["gt_in_cell"] = round(100 * float(np.mean([h[0] for h in hits])), 1)
        row["gt_alone"] = round(100 * float(np.mean([h[1] for h in hits])), 1)
    return row, hits


def contrast_table(region, markers, raw, inputs):
    rows = []
    for t, marker in enumerate(markers):
        x = raw[t].astype(np.float64)
        others = [i for i in range(len(markers)) if i != t]
        others_dim = np.all([raw[i] < np.percentile(raw[i], 90) for i in others], axis=0)
        pos = others_dim & (x > np.percentile(x, 99))
        neg = others_dim & (x < np.percentile(x, 50))
        if pos.sum() < 20:
            continue

        def contrast(img):
            img = img.astype(np.float64)
            return round(float((img[pos].mean() - img[neg].mean()) / img[neg].std()), 2)

        rows.append({"region": region, "marker": marker, "pixels": int(pos.sum()),
                     **{m: contrast(inputs[m]) for m in METHODS}, "marker_alone": contrast(x)})
    return rows


def pick_window(points, shape):
    """Top-left corner of the WINDOW x WINDOW window holding the most points."""
    density = np.zeros(shape)
    for r, c in points:
        density[r, c] += 1
    density = ndi.uniform_filter(density, size=WINDOW, mode="constant")
    half = WINDOW // 2
    inner = density[half:shape[0] - half, half:shape[1] - half]
    return np.unravel_index(np.argmax(inner), inner.shape)


def figure(path, dapi, markers_raw, inputs, masks, segmenters, corner, gt):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    y, x = corner
    win = (slice(y, y + WINDOW), slice(x, x + WINDOW))

    def unit(a):
        lo, hi = np.percentile(a, [1, 99.5])
        return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)

    # Same background in every segmentation panel: rescaled-max membrane in
    # red, DAPI in blue. Ground-truth centroids are yellow crosses.
    membrane = np.max([np.clip(rescale(p[win]), 0, 1) for p in markers_raw], axis=0)
    background = np.dstack([membrane, membrane * 0.25, unit(dapi[win])])
    gt_win = [(r - y, c - x) for r, c, _ in gt if y <= r < y + WINDOW and x <= c < x + WINDOW]

    fig, axes = plt.subplots(len(METHODS), 1 + len(segmenters),
                             figsize=(3 * (1 + len(segmenters)), 3 * len(METHODS)))
    for r, method in enumerate(METHODS):
        axes[r, 0].imshow(unit(inputs[method][win]), cmap="gray")
        axes[r, 0].set_ylabel(method, fontsize=11)
        for c, seg in enumerate(segmenters, start=1):
            panel = background.copy()
            panel[find_boundaries(masks[(method, seg)][win], mode="inner")] = (1, 1, 1)
            axes[r, c].imshow(panel)
            if gt_win:
                axes[r, c].scatter([p[1] for p in gt_win], [p[0] for p in gt_win],
                                   s=6, c="yellow", marker="x", linewidths=0.8)
    for c, title in enumerate(["membrane input"] + list(segmenters)):
        axes[0, c].set_title(title, fontsize=11)
    for ax in axes.ravel():
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    inputs_dir, out = Path(args.inputs), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    markers = cfg["markers"]
    segmenters = cfg.get("segmenters", ["cellpose", "mesmer", "cellsam"])

    scores, contrasts = [], []
    by_type = defaultdict(lambda: [0, 0, 0])  # (method, seg, label) -> n, in_cell, alone
    for region in cfg["regions"]:
        name = region["name"]
        source_path = inputs_dir / f"{name}_source.ome.tif"
        source = tifffile.imread(source_path)
        dapi, raw = source[0], source[1:]
        pixel_um = (source_physical_size(source_path) or (0.28,))[0]
        nuclei = nuclei_centres(dapi, pixel_um)
        gt = []
        if cfg.get("ground_truth"):
            gt = load_ground_truth(cfg["ground_truth"], Path(region["tiff"]).name,
                                   region.get("y", 0), region.get("x", 0),
                                   region.get("size", max(dapi.shape)))
        inputs = {m: tifffile.imread(inputs_dir / f"{name}_{m}.ome.tif")[1] for m in METHODS}
        masks = {}
        for method in METHODS:
            for seg in segmenters:
                mask = tifffile.imread(Path(args.results) / "cellmeasurement"
                                       / f"{name}_{method}_{seg}_mask.tiff")
                masks[(method, seg)] = mask
                row, hits = score(mask, pixel_um ** 2, nuclei, gt)
                scores.append({"region": name, "method": method, "segmenter": seg, **row})
                for (_, _, label), (in_cell, alone) in zip(gt, hits):
                    for key in ((method, seg, label), (method, seg, "all")):
                        by_type[key][0] += 1
                        by_type[key][1] += bool(in_cell)
                        by_type[key][2] += bool(alone)
        contrasts += contrast_table(name, markers, raw, inputs)
        points = [(r, c) for r, c, _ in gt] or [tuple(p) for p in nuclei]
        corner = pick_window(points, dapi.shape)
        figure(out / f"overlays_{name}.png", dapi, raw, inputs, masks, segmenters, corner, gt)
        print(f"{name}: {len(nuclei)} nuclei, {len(gt)} ground-truth cells")

    write_csv(out / "scores.csv", scores)
    if contrasts:
        write_csv(out / "contrast.csv", contrasts)
    if by_type:
        write_csv(out / "gt_by_type.csv", [
            {"method": m, "segmenter": s, "label": lab, "n": n,
             "gt_in_cell": round(100 * k / n, 1), "gt_alone": round(100 * a / n, 1)}
            for (m, s, lab), (n, k, a) in sorted(by_type.items())
        ])


if __name__ == "__main__":
    main()
