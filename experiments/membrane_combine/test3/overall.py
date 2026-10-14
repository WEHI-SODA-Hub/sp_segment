"""One table per method: cells, low-DAPI cells, labelled cells found, and the shape scores.

low_dapi: cells whose 90th-percentile DAPI < 0.2 x the median production cell's (per crop, fixed scale)
gt_in_cell / gt_alone: all 484 labelled cells (as analyse.py)
shape columns joined from out/shape_summary.csv (118 cells with a usable marker reference)
python overall.py   (run after score_shape.py; reuses its load_masks)
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
from score_shape import CFG, INP, OUT, PX, load_ground_truth, load_masks  # noqa: E402


def main():
    tot = defaultdict(lambda: np.zeros(5))   # cells, low_dapi, gt_n, gt_in, gt_alone
    for region in CFG["regions"]:
        crop = region["name"]
        dapi = tifffile.imread(INP / f"{crop}_source.ome.tif")[0].astype(np.float32)
        gt = load_ground_truth(CFG["ground_truth"], Path(region["tiff"]).name,
                               region["y"], region["x"], region["size"])
        masks = load_masks(crop, dapi)
        prod = masks["production"]
        pl = np.unique(prod)[1:]
        ref = np.median(ndi.labeled_comprehension(dapi, prod, pl, lambda v: np.percentile(v, 90), float, 0))
        for name, mask in masks.items():
            labs = np.unique(mask)
            labs = labs[labs > 0]
            p90 = np.array(ndi.labeled_comprehension(dapi, mask, labs, lambda v: np.percentile(v, 90), float, 0))
            at = np.array([mask[r, c] for r, c, _ in gt])
            counts = np.bincount(at[at > 0], minlength=int(mask.max()) + 1) if (at > 0).any() else None
            alone = sum(1 for a in at if a > 0 and counts[a] == 1)
            tot[name] += [len(labs), (p90 < 0.2 * ref).sum(), len(gt), (at > 0).sum(), alone]
        print(crop, "done", flush=True)
    shape = {r["method"]: r for r in csv.DictReader(open(OUT / "shape_summary.csv"))}
    rows = []
    for name, (c, low, n, k, a) in tot.items():
        s = shape.get(name, {})
        rows.append({"method": name, "cells": int(c), "low_dapi": int(low),
                     "low_dapi_pct": round(100 * low / c, 1),
                     "gt_in_cell_pct": round(100 * k / n, 1), "gt_alone_pct": round(100 * a / n, 1),
                     "shape_mean_iou": s.get("mean_iou", ""), "shape_iou_found": s.get("mean_iou_found", ""),
                     "shape_recall_found": s.get("median_recall_found", ""),
                     "shape_precision": s.get("median_precision", ""),
                     "median_cell_area_um2": s.get("median_cell_area_um2", "")})
    rows.sort(key=lambda r: -float(r["shape_mean_iou"] or 0))
    with open(OUT / "overall.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{'method':38s} {'cells':>6s} {'lowDAPI%':>8s} {'GTin':>6s} {'GTalone':>7s} {'IoU':>6s} {'IoUfnd':>6s} {'recall':>6s} {'prec':>6s} {'area':>6s}")
    for r in rows:
        print(f"{r['method']:38s} {r['cells']:6d} {r['low_dapi_pct']:8.1f} {r['gt_in_cell_pct']:6.1f} {r['gt_alone_pct']:7.1f} "
              f"{r['shape_mean_iou']:>6s} {r['shape_iou_found']:>6s} {r['shape_recall_found']:>6s} {r['shape_precision']:>6s} {r['median_cell_area_um2']:>6s}")


if __name__ == "__main__":
    main()
