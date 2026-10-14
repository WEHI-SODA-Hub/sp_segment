"""Score cell outlines against a marker-based reference outline for each labelled cell.

Reference outline of a labelled cell (centroid from ground_truth_8_10_25.csv):
  1. its type's marker (MARKERS), Gaussian-smoothed (sigma 1 px), thresholded at
     max(Otsu in an 86 x 86 px window around the cell, Otsu of the whole crop), both on log1p
  2. union with DAPI foreground (crop Otsu on log1p), holes filled
  3. pixels within 12 um of the centroid and nearer to it than to any other nucleus
     (LoG nucleus centres more than 3 um away, and other labelled centroids)
  4. the connected piece that contains the centroid
Cells whose reference is < 20% marker-positive are skipped (marker too weak to give a shape).

For each segmentation, the cell under the centroid is compared with the reference:
  iou, recall = |cell & ref| / |ref|, precision = |cell & ref| / |cell|,
  capture = marker signal in cell & ref / marker signal in ref. A missed cell scores 0.
Writes out/shape_cells.csv, out/shape_summary.csv, out/shape_by_type.csv and
out/ref_gallery_<type>.png (reference cyan, norm_mean red, raw_max yellow).
python score_shape.py
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import tifffile  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402
from skimage.filters import gaussian, threshold_otsu  # noqa: E402
from skimage.segmentation import expand_labels, find_boundaries  # noqa: E402

REPO = Path("/vast/scratch/users/mckay.m/sp_segment")
sys.path[:0] = [str(REPO / "bin"), str(REPO / "experiments/membrane_combine")]
from analyse import load_ground_truth, nuclei_centres  # noqa: E402

T = Path("/vast/scratch/users/mckay.m/sp_segment_tests")
PIPE = T / "membrane_combine_all/results/cellmeasurement"
VAR = T / "cell_shape/masks"
INP = T / "membrane_combine_csl/inputs"
PROD = T / "membrane_artefacts/out"
OUT = T / "cell_shape/out"
CFG = json.loads((REPO / "experiments/membrane_combine/config_csl.json").read_text())
CHANNELS = [CFG["nuclear"]] + CFG["markers"]
MARKERS = {"Macrophages": "R2_CSL_CD68", "Monocytes": "R1_CSL_CD45", "Neutrophils": "R1_CSL_S100A9",
           "Tcell": "R1_CSL_CD3", "Endothelial": "R2_CSL_VECadherin", "VSMCs": "R1_CSL_aSMA",
           "Fibroblasts": "R1_CSL_FAP"}
PX = 0.28
RADIUS = int(round(12 / PX))   # 43 px
HALF = RADIUS + 2


def log_otsu(values):
    v = values[values > 0]
    return float(np.expm1(threshold_otsu(np.log1p(v)))) if v.size > 10 else 0.0


def load_masks(crop, dapi):
    masks = {}
    for inp in ["raw_max", "norm_max", "norm_mean", "single"]:
        for seg in ["cellpose", "mesmer", "cellsam"]:
            masks[f"{inp} {seg} (pipeline)"] = tifffile.imread(PIPE / f"{crop}_{inp}_{seg}_mask.tiff")
    masks["production"] = np.load(PROD / f"{crop}_production_mask.npy")
    for inp in ["norm_mean", "raw_max"]:
        for v in ["base", "d40", "d50", "cp-1", "cp-2"]:
            p = VAR / f"{crop}_{inp}_{v}.tiff"
            if p.exists():
                masks[f"{inp} cellpose {v}"] = tifffile.imread(p)
    for v in ["raw_max_base_check", "raw_max_f2_cp0", "raw_max_f2_cp-2", "raw_max_f2_cp-4", "raw_max_f2_cp-6",
              "raw_max_clahe", "raw_max_clahe_f2_cp-2", "raw_max_sopa_default", "3ch", "3ch_f2_cp-2",
              "3ch_f2_cp-4", "3ch_clahe_f2_cp-2", "norm_mean_f2_cp-2"]:
        p = VAR / f"{crop}_{v}.tiff"
        if p.exists():
            masks[v] = tifffile.imread(p)
    p = VAR / f"{crop}_nuclei.tiff"
    if p.exists():
        nuc = tifffile.imread(p)
        masks["nuclei only"] = nuc
        for um in (1.5, 3.0, 4.5):
            masks[f"nuclei + {um} um"] = expand_labels(nuc, distance=int(round(um / PX)))
    # hybrid: norm_mean outlines + raw-max cells with a clear nucleus where norm_mean has no cell
    nm = masks["norm_mean cellpose (pipeline)"].astype(np.int64)
    raw = masks["raw_max cellpose (pipeline)"].astype(np.int64)
    prod = masks["production"]
    pl = np.unique(prod)[1:]
    ref = np.median(ndi.labeled_comprehension(dapi, prod, pl, lambda v: np.percentile(v, 90), float, 0))
    rl = np.unique(raw)[1:]
    p90 = ndi.labeled_comprehension(dapi, raw, rl, lambda v: np.percentile(v, 90), float, 0)
    cover = ndi.mean((nm > 0).astype(float), raw, rl)
    add = rl[(p90 >= 0.5 * ref) & (cover < 0.2)]
    hyb = nm.copy()
    fill = np.isin(raw, add) & (nm == 0)
    hyb[fill] = raw[fill] + nm.max()
    masks["hybrid: norm_mean + raw_max fill-in"] = hyb
    return masks


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cells, examples = [], defaultdict(list)
    for region in CFG["regions"]:
        crop = region["name"]
        src = tifffile.imread(INP / f"{crop}_source.ome.tif").astype(np.float32)
        dapi = src[0]
        nuclei = nuclei_centres(src[0].astype(np.uint16), PX)
        gt = load_ground_truth(CFG["ground_truth"], Path(region["tiff"]).name,
                               region["y"], region["x"], region["size"])
        dsm = gaussian(dapi, 1, preserve_range=True)
        dfg = dsm > log_otsu(dsm)
        smooth, glob = {}, {}
        for lab, marker in MARKERS.items():
            if marker not in smooth:
                smooth[marker] = gaussian(src[CHANNELS.index(marker)], 1, preserve_range=True)
                glob[marker] = log_otsu(smooth[marker])
        masks = load_masks(crop, dapi)
        areas = {k: np.bincount(v.ravel().astype(np.int64)) for k, v in masks.items()}
        H, W = dapi.shape
        gt_pts = np.array([(r, c) for r, c, _ in gt])
        for gi, (r, c, lab) in enumerate(gt):
            if lab not in MARKERS:
                continue
            marker = MARKERS[lab]
            y0, y1 = max(r - HALF, 0), min(r + HALF + 1, H)
            x0, x1 = max(c - HALF, 0), min(c + HALF + 1, W)
            m = smooth[marker][y0:y1, x0:x1]
            thr = max(log_otsu(m), glob[marker])
            pos = m > thr
            union = ndi.binary_fill_holes(pos | dfg[y0:y1, x0:x1])
            yy, xx = np.mgrid[y0:y1, x0:x1]
            d0 = np.hypot(yy - r, xx - c)
            others = [p for p in nuclei if np.hypot(p[0] - r, p[1] - c) > 3 / PX]
            others += [tuple(p) for j, p in enumerate(gt_pts) if j != gi]
            others = np.array([p for p in others
                               if y0 - HALF <= p[0] < y1 + HALF and x0 - HALF <= p[1] < x1 + HALF])
            dmin = np.full(d0.shape, np.inf)
            for p in others:
                dmin = np.minimum(dmin, np.hypot(yy - p[0], xx - p[1]))
            region_ok = union & (d0 <= RADIUS) & (d0 <= dmin)
            lbl, _ = ndi.label(region_ok)
            ry, rx = r - y0, c - x0
            k = lbl[ry, rx]
            if k == 0:
                idx = np.argwhere(lbl > 0)
                if len(idx) == 0:
                    continue
                near = idx[np.argmin(np.hypot(idx[:, 0] - ry, idx[:, 1] - rx))]
                if np.hypot(*(near - [ry, rx])) > 3 / PX:
                    continue
                k = lbl[tuple(near)]
            refm = lbl == k
            pos_frac = float(pos[refm].mean())
            # share of the reference outline set by the marker/DAPI edge rather than by the
            # 12 um limit or the midline to a neighbouring nucleus
            edge = find_boundaries(refm, mode="inner")
            cut = (d0 >= RADIUS - 1.5) | (np.abs(d0 - dmin) <= 1.5)
            natural = float((~cut[edge]).mean()) if edge.any() else 0.0
            if pos_frac < 0.2 or natural < 0.5:
                cells.append({"crop": crop, "label": lab, "method": "_skipped", "ref_area_um2":
                              round(refm.sum() * PX ** 2, 1), "pos_frac": round(pos_frac, 2),
                              "natural_edge": round(natural, 2)})
                continue
            sig = smooth[marker][y0:y1, x0:x1]
            for name, mask in masks.items():
                lab_id = mask[r, c]
                if lab_id == 0:
                    row = dict(iou=0.0, recall=0.0, precision=np.nan, capture=0.0, cell_area_um2=0.0)
                else:
                    cellm = mask[y0:y1, x0:x1] == lab_id
                    area = int(areas[name][lab_id])
                    inter = np.count_nonzero(cellm & refm)
                    row = dict(iou=inter / (area + refm.sum() - inter), recall=inter / refm.sum(),
                               precision=inter / area, capture=float(sig[cellm & refm].sum() / sig[refm].sum()),
                               cell_area_um2=area * PX ** 2)
                cells.append({"crop": crop, "label": lab, "method": name, "x": c, "y": r,
                              "ref_area_um2": round(refm.sum() * PX ** 2, 1), "pos_frac": round(pos_frac, 2),
                              **{k2: round(v, 3) for k2, v in row.items()}})
            if len(examples[lab]) < 12:
                examples[lab].append((crop, refm, (y0, y1, x0, x1),
                                      masks["norm_mean cellpose (pipeline)"][y0:y1, x0:x1],
                                      masks["raw_max cellpose (pipeline)"][y0:y1, x0:x1],
                                      dapi[y0:y1, x0:x1], m, (ry, rx)))
        print(crop, "done", flush=True)

    with open(OUT / "shape_cells.csv", "w", newline="") as fh:
        keys = sorted({k for row in cells for k in row}, key=lambda k: list(cells[-1]).index(k)
                      if k in cells[-1] else 99)
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(cells)

    scored = [x for x in cells if x["method"] != "_skipped"]
    skipped = defaultdict(int)
    for x in cells:
        if x["method"] == "_skipped":
            skipped[x["label"]] += 1

    def summarise(rows):
        iou = np.array([x["iou"] for x in rows])
        rec = np.array([x["recall"] for x in rows])
        pre = np.array([x["precision"] for x in rows], float)
        cap = np.array([x["capture"] for x in rows])
        area = np.array([x["cell_area_um2"] for x in rows])
        refa = np.array([x["ref_area_um2"] for x in rows])
        found = area > 0
        return {"n": len(rows), "found_pct": round(100 * float(found.mean()), 1),
                "mean_iou": round(float(iou.mean()), 3),
                "mean_iou_found": round(float(iou[found].mean()), 3) if found.any() else 0,
                "median_recall_found": round(float(np.median(rec[found])), 3) if found.any() else 0, "median_recall": round(float(np.median(rec)), 3),
                "median_precision": round(float(np.nanmedian(pre)), 3),
                "mean_capture": round(float(cap.mean()), 3),
                "median_cell_area_um2": round(float(np.median(area[area > 0])), 1) if (area > 0).any() else 0,
                "median_ref_area_um2": round(float(np.median(refa)), 1)}

    by_method, by_type = defaultdict(list), defaultdict(list)
    for x in scored:
        by_method[x["method"]].append(x)
        by_type[(x["method"], x["label"])].append(x)
    summary = sorted(({"method": k, **summarise(v)} for k, v in by_method.items()),
                     key=lambda d: -d["mean_iou"])
    for name, rows in (("shape_summary.csv", summary),
                       ("shape_by_type.csv", [{"method": m, "label": lab, **summarise(v)}
                                              for (m, lab), v in sorted(by_type.items())])):
        with open(OUT / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    print("skipped (marker too weak):", dict(skipped))
    kept = defaultdict(int)
    for x in by_method["raw_max cellpose (pipeline)"]:
        kept[x["label"]] += 1
    print("scored per type:", dict(kept))
    print(f"{'method':42s} {'n':>4s} {'found':>6s} {'IoU':>6s} {'IoUfnd':>6s} {'recF':>6s} {'prec':>6s} {'capt':>6s} {'area':>6s} {'ref':>6s}")
    for s in summary:
        print(f"{s['method']:42s} {s['n']:4d} {s['found_pct']:6.1f} {s['mean_iou']:6.3f} {s['mean_iou_found']:6.3f} {s['median_recall_found']:6.3f} "
              f"{s['median_precision']:6.3f} {s['mean_capture']:6.3f} {s['median_cell_area_um2']:6.1f} {s['median_ref_area_um2']:6.1f}")

    for lab, exs in examples.items():
        fig, axes = plt.subplots(3, 8, figsize=(16, 6.5))
        for ax in axes.ravel():
            ax.axis("off")
        for i, (crop, refm, _, nm, raw, d, m, (ry, rx)) in enumerate(exs):
            for j, (img, title) in enumerate(((m, MARKERS[lab].split("_")[-1]), (d, "DAPI"))):
                ax = axes.ravel()[2 * i + j] if 2 * i + j < axes.size else None
                if ax is None:
                    continue
                g = np.clip(img / max(np.percentile(img, 99.5), 1), 0, 1)
                rgb = np.stack([g, g, g], -1)
                rgb[find_boundaries(refm, mode="inner")] = [0, 1, 1]
                if nm[ry, rx]:
                    rgb[find_boundaries(nm == nm[ry, rx], mode="inner")] = [1, 0, 0]
                if raw[ry, rx]:
                    rgb[find_boundaries(raw == raw[ry, rx], mode="inner")] = [1, 1, 0]
                ax.imshow(rgb)
                ax.set_title(f"{crop} {title}", fontsize=7)
        fig.suptitle(f"{lab}: reference outline (cyan), norm_mean cell (red), raw_max cell (yellow); "
                     f"window 24 um", fontsize=10)
        fig.tight_layout()
        fig.savefig(OUT / f"ref_gallery_{lab}.png", dpi=80)
        plt.close(fig)


if __name__ == "__main__":
    main()
