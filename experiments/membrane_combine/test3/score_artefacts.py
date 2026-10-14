"""Count cells without a nucleus (likely artefacts) for each Cellpose membrane input.

Per cell of each mask:
  nuclei     LoG nucleus centres inside the cell (same detector as analyse.py)
  dapi_frac  share of the cell's pixels above the crop's DAPI Otsu threshold
  contrast   mean raw-max membrane in the cell / mean in a 3 px ring outside it
             (ring pixels that belong to other cells are excluded; NaN if < 10 px)
  prod_iou   best IoU with a production cell
A cell is "no nucleus" when nuclei == 0 and dapi_frac < 0.1.

python score_artefacts.py <outdir>
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage.draw import polygon as draw_polygon
from skimage.filters import gaussian, threshold_otsu
from skimage.segmentation import find_boundaries

REPO = Path("/vast/scratch/users/mckay.m/sp_segment")
sys.path[:0] = [str(REPO / "bin"), str(REPO / "experiments/membrane_combine")]
from analyse import nuclei_centres  # noqa: E402

T = Path("/vast/scratch/users/mckay.m/sp_segment_tests")
INPUTS = T / "membrane_combine_csl/inputs"
MAIN = T / "membrane_combine_all/results/cellmeasurement"
STROMAL = T / "membrane_combine_stromal/results/cellmeasurement"
PROD = Path("/vast/scratch/users/mckay.m/csl_composite_check/membrane_test/production")
CROPS = ["b6_immune", "b6_macro", "b2_endo", "b2_vsmc", "b2_tcell"]
PIXEL_UM = 0.28


def mask_sources(crop):
    """(label, mask path) for every Cellpose mask of this crop."""
    out = [(m, MAIN / f"{crop}_{m}_cellpose_mask.tiff")
           for m in ("raw_max", "norm_max", "norm_mean", "single")]
    for tag, alone in (("fap", "FAP alone"), ("asma", "aSMA alone")):
        names = {"raw_max": "FAP+aSMA max", "norm_max": "FAP+aSMA rescaled max", "single": alone}
        for m, label in names.items():
            p = STROMAL / f"{crop}_{tag}_{m}_cellpose_mask.tiff"
            if p.exists():
                out.append((label, p))
    return out


def production_mask(crop, shape):
    feats = json.loads((PROD / f"{crop}_production.geojson").read_text())
    feats = feats["features"] if isinstance(feats, dict) else feats
    mask = np.zeros(shape, np.int32)
    for i, f in enumerate(feats, 1):
        g = f["geometry"]
        polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
        for poly in polys:
            ring = np.asarray(poly[0])
            rr, cc = draw_polygon(ring[:, 1], ring[:, 0], shape)
            mask[rr, cc] = i
    return mask


def best_iou(mask, ref):
    """Best IoU of each label in mask with any label in ref."""
    n = int(mask.max())
    area = np.bincount(mask.ravel(), minlength=n + 1)
    ref_area = np.bincount(ref.ravel(), minlength=int(ref.max()) + 1)
    both = (mask > 0) & (ref > 0)
    pairs, inter = np.unique(np.stack([mask[both], ref[both]]), axis=1, return_counts=True)
    best = np.zeros(n + 1)
    for (a, b), k in zip(pairs.T, inter):
        best[a] = max(best[a], k / (area[a] + ref_area[b] - k))
    return best


def cell_table(mask, dapi_fg, nuclei, membrane):
    labels = np.arange(1, int(mask.max()) + 1)
    present = np.bincount(mask.ravel(), minlength=len(labels) + 1)[1:] > 0
    labels = labels[present]
    area = ndi.sum(np.ones_like(mask), mask, labels)
    dapi_frac = ndi.mean(dapi_fg, mask, labels)
    nuc = np.bincount(mask[nuclei[:, 0], nuclei[:, 1]], minlength=int(mask.max()) + 1)[labels]
    inner = ndi.mean(membrane, mask, labels)
    # ring: 3 px dilation of each cell, minus all cells
    grown = ndi.grey_dilation(mask, size=7)
    ring = np.where((mask == 0) & (grown > 0), grown, 0)
    ring_n = ndi.sum(np.ones_like(ring), ring, labels)
    ring_mean = ndi.mean(membrane, ring, labels)
    contrast = np.where(ring_n >= 10, inner / np.maximum(ring_mean, 1e-6), np.nan)
    return labels, area, dapi_frac, nuc, contrast


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    summary, cells_rows = [], []
    for crop in CROPS:
        source = tifffile.imread(INPUTS / f"{crop}_source.ome.tif")
        dapi, membrane = source[0].astype(np.float32), source[1:].max(0).astype(np.float32)
        nuclei = nuclei_centres(source[0], PIXEL_UM)
        sm = gaussian(dapi, 1, preserve_range=True)
        thr = float(np.expm1(threshold_otsu(np.log1p(sm))))
        dapi_fg = (sm > thr).astype(np.float32)
        prod = production_mask(crop, dapi.shape)
        np.save(out / f"{crop}_production_mask.npy", prod)
        sources = mask_sources(crop) + [("production", None)]
        for label, path in sources:
            mask = prod if path is None else tifffile.imread(path).astype(np.int32)
            labels, area, dfrac, nuc, contrast = cell_table(mask, dapi_fg, nuclei, membrane)
            iou = best_iou(mask, prod)[labels] if path is not None else np.ones(len(labels))
            no_nuc = (nuc == 0) & (dfrac < 0.1)
            low_c = contrast < 1.2
            n = len(labels)
            summary.append({
                "crop": crop, "input": label, "cells": n,
                "median_area_um2": round(float(np.median(area)) * PIXEL_UM ** 2, 1),
                "no_nucleus": int(no_nuc.sum()),
                "no_nucleus_pct": round(100 * no_nuc.mean(), 1),
                "no_nucleus_area_pct": round(100 * area[no_nuc].sum() / area.sum(), 1),
                "no_nucleus_low_contrast": int((no_nuc & low_c).sum()),
                "nucleated": int((~no_nuc).sum()),
                "match_production_pct": round(100 * (iou >= 0.5).mean(), 1),
                "no_nucleus_not_in_production": int((no_nuc & (iou < 0.5)).sum()),
                "nuclei_in_crop": len(nuclei), "dapi_otsu": round(thr, 1),
            })
            for i in range(n):
                cells_rows.append({
                    "crop": crop, "input": label, "label": int(labels[i]),
                    "area_um2": round(float(area[i]) * PIXEL_UM ** 2, 1),
                    "nuclei": int(nuc[i]), "dapi_frac": round(float(dfrac[i]), 3),
                    "contrast": round(float(contrast[i]), 2), "prod_iou": round(float(iou[i]), 2),
                    "no_nucleus": bool(no_nuc[i]),
                })
            if path is not None:
                np.save(out / f"{crop}_{label.replace(' ', '_').replace('+', '')}_nonuc.npy",
                        labels[no_nuc])
        print(crop, "done", flush=True)
    for name, rows in (("summary.csv", summary), ("cells.csv", cells_rows)):
        with open(out / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    # totals over crops for the inputs present in all five
    tot = {}
    for r in summary:
        t = tot.setdefault(r["input"], [0, 0, 0, 0])
        t[0] += r["cells"]; t[1] += r["no_nucleus"]; t[2] += 1
        t[3] += r["no_nucleus_not_in_production"]
    for k, (c, nn, ncrop, nnp) in tot.items():
        print(f"{k:28s} crops={ncrop} cells={c:6d} no_nucleus={nn:5d} ({100 * nn / c:.1f}%) "
              f"not_in_production={nnp}")


if __name__ == "__main__":
    main()
