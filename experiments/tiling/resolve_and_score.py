#!/usr/bin/env python3
"""
Stitch the per-patch Cellpose polygons two ways and score them.

Stitching:
  sopa   sopa.segmentation.solve_conflicts(threshold=0.5), as the pipeline does:
         overlapping cells from different patches are merged into their union
         when the overlap is >= 0.5 x the smaller cell, otherwise both are kept.
  owner  each cell is taken from one patch: a patch keeps the cells whose
         centroid lies in its core (the patch minus half the overlap on each
         inner side) and that do not touch an inner patch edge. Remaining
         overlaps go to the cell with the higher mean cell probability; the
         other cell loses the overlap and is dropped if it keeps < 50% of its area.

Scores, per normalisation and stitching:
  cells, overlapping cell pairs left, straight edges (polygon segments >= 15 px
  that are exactly horizontal or vertical, not on the region border), how many
  of them lie on a patch edge, cells with a straight edge, and labelled cells
  inside a cell / alone in a cell, for labels near a seam (<= 150 px) and elsewhere.
Also the agreement of the two patches' outlines for cells in the shared strips.

Run inside the sopa-cellpose container:
  python resolve_and_score.py <patches dir> <out csv>
"""
import csv, io, json, sys
from pathlib import Path
import numpy as np, geopandas as gpd, shapely
from shapely.strtree import STRtree
from sopa.segmentation import solve_conflicts

D = Path(sys.argv[1]); OUT = Path(sys.argv[2])
meta = json.loads((D / "meta.json").read_text())
patches = np.array(meta["patches"], float)
X0, Y0, X1, Y1 = meta["x0"], meta["y0"], meta["x1"], meta["y1"]
half = meta["overlap"] / 2
IMAGE = "20260826_163033_3_k1ghXo_CSL_batch6_cohort_25_1592_1A_26_948_1_Y_Z_33ch_AS_1.ome.tif"
GT_FILES = ["/vast/projects/CSL_COMET/csl_artherosclerosis/ground_truth_8_10_25.csv",
            "/vast/projects/CSL_COMET/csl_artherosclerosis/ground_truth_6_10_26.csv"]

inner_x = sorted({v for p in patches for v in (p[0], p[2]) if X0 < v < X1})
inner_y = sorted({v for p in patches for v in (p[1], p[3]) if Y0 < v < Y1})
seam_x = sorted({(a + b) / 2 for a in inner_x for b in inner_x if 0 < b - a <= meta["overlap"]})
seam_y = sorted({(a + b) / 2 for a in inner_y for b in inner_y if 0 < b - a <= meta["overlap"]})


def load(norm):
    parts = [gpd.read_parquet(f) for f in sorted((D / norm).glob("patch_*.parquet"))]
    return gpd.GeoDataFrame(pd_concat(parts), geometry="geometry")


def pd_concat(parts):
    import pandas as pd
    return pd.concat(parts, ignore_index=True)


def core_box(p):
    x0, y0, x1, y1 = p
    return (x0 + (half if x0 > X0 else 0), y0 + (half if y0 > Y0 else 0),
            x1 - (half if x1 < X1 else 0), y1 - (half if y1 < Y1 else 0))


def owner(cells):
    keep = []
    for k, p in enumerate(patches):
        cx0, cy0, cx1, cy1 = core_box(p)
        sub = cells[cells.patch == k]
        c = sub.geometry.centroid
        b = sub.geometry.bounds
        in_core = (c.x >= cx0) & (c.x < cx1) & (c.y >= cy0) & (c.y < cy1)
        touches = (((b.minx <= p[0] + 1) & (p[0] > X0)) | ((b.maxx >= p[2] - 1) & (p[2] < X1)) |
                   ((b.miny <= p[1] + 1) & (p[1] > Y0)) | ((b.maxy >= p[3] - 1) & (p[3] < Y1)))
        keep.append(sub[in_core & ~touches])
    kept = gpd.GeoDataFrame(pd_concat(keep), geometry="geometry")
    kept = kept.sort_values("mean_cellprob", ascending=False).reset_index(drop=True)
    geoms = list(kept.geometry)
    tree = STRtree(geoms)
    accepted, out = np.zeros(len(geoms), bool), []
    for i, g in enumerate(geoms):
        others = [j for j in tree.query(g, predicate="intersects") if accepted[j] and j != i]
        cell = g
        if others:
            cell = g.difference(shapely.union_all([geoms[j] for j in others]))
            if cell.area < 0.5 * g.area:
                continue
            if cell.geom_type != "Polygon":
                cell = max(cell.geoms, key=lambda q: q.area)
        accepted[i] = True
        geoms[i] = cell
        out.append(cell)
    return gpd.GeoDataFrame(geometry=out)


def edge_info(cells):
    """distance from each cell centroid to the nearest inner edge of its own patch, and whether it touches one"""
    dist, touch = np.zeros(len(cells)), np.zeros(len(cells), bool)
    c = cells.geometry.centroid; b = cells.geometry.bounds
    for k, p in enumerate(patches):
        sel = (cells.patch == k).values
        d = np.full(sel.sum(), 1e9)
        t = np.zeros(sel.sum(), bool)
        if p[0] > X0: d = np.minimum(d, c.x[sel] - p[0]); t |= (b.minx[sel] <= p[0] + 1).values
        if p[2] < X1: d = np.minimum(d, p[2] - c.x[sel]); t |= (b.maxx[sel] >= p[2] - 1).values
        if p[1] > Y0: d = np.minimum(d, c.y[sel] - p[1]); t |= (b.miny[sel] <= p[1] + 1).values
        if p[3] < Y1: d = np.minimum(d, p[3] - c.y[sel]); t |= (b.maxy[sel] >= p[3] - 1).values
        dist[sel], touch[sel] = d, t
    return dist, touch


def match(cells, by):
    """group versions of a cell from different patches (overlap >= 0.5 x smaller), keep one per group"""
    cells = cells.reset_index(drop=True)
    g = list(cells.geometry); pt = cells.patch.values
    dist, touch = edge_info(cells)
    parent = np.arange(len(g))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    tree = STRtree(g)
    for i, j in tree.query(g, predicate="intersects").T:
        if i < j and pt[i] != pt[j] and g[i].intersection(g[j]).area >= 0.5 * min(g[i].area, g[j].area):
            parent[find(i)] = find(j)
    groups = {}
    for i in range(len(g)):
        groups.setdefault(find(i), []).append(i)
    keep = []
    for members in groups.values():
        if by == "context":   # not cut by a patch edge first, then most context
            best = max(members, key=lambda i: (not touch[i], dist[i], cells.mean_cellprob[i]))
        else:                 # highest mean cell probability
            best = max(members, key=lambda i: (cells.mean_cellprob[i], not touch[i]))
        keep.append(best)
    kept = cells.iloc[keep].sort_values("mean_cellprob", ascending=False).reset_index(drop=True)
    geoms = list(kept.geometry); tree = STRtree(geoms)
    accepted, out = np.zeros(len(geoms), bool), []
    for i, gg in enumerate(geoms):      # leftover overlaps: the higher cell probability keeps the pixels
        others = [j for j in tree.query(gg, predicate="intersects") if accepted[j] and j != i]
        cell = gg
        if others:
            cell = gg.difference(shapely.union_all([geoms[j] for j in others]))
            if cell.is_empty or cell.area < 0.5 * gg.area:
                continue
            if cell.geom_type != "Polygon":
                cell = max(cell.geoms, key=lambda q: q.area)
        accepted[i] = True; geoms[i] = cell; out.append(cell)
    return gpd.GeoDataFrame(geometry=out)


def straight_stats(geoms):
    total = on_edge = 0
    cells_with = 0
    for poly in geoms:
        if poly.geom_type != "Polygon":
            poly = max(poly.geoms, key=lambda q: q.area)
        r = np.asarray(poly.exterior.coords)
        s = np.diff(r, axis=0)
        L = np.hypot(*s.T)
        has = False
        for k in np.where(L >= 15)[0]:
            if abs(s[k, 1]) < 1e-6:
                y = r[k, 1]
                if Y0 + 2 < y < Y1 - 2:
                    total += 1; has = True
                    on_edge += min((abs(y - e) for e in inner_y), default=99) <= 1.5
            elif abs(s[k, 0]) < 1e-6:
                x = r[k, 0]
                if X0 + 2 < x < X1 - 2:
                    total += 1; has = True
                    on_edge += min((abs(x - e) for e in inner_x), default=99) <= 1.5
        cells_with += has
    return total, int(on_edge), cells_with


def overlaps(geoms):
    tree = STRtree(geoms)
    pairs = tree.query(geoms, predicate="intersects")
    n = 0
    for i, j in pairs.T:
        if i < j and geoms[i].intersection(geoms[j]).area > 1:
            n += 1
    return n


def gt_points():
    pts = []
    for f in GT_FILES:
        lines = [l for l in open(f) if not l.startswith("#")]
        for r in csv.DictReader(io.StringIO("".join(lines))):
            if r["Image"] == IMAGE:
                x, y = float(r["CentroidX"]), float(r["CentroidY"])
                if X0 <= x < X1 and Y0 <= y < Y1:
                    pts.append((x, y))
    return pts


def gt_scores(geoms, pts):
    tree = STRtree(geoms)
    res = {"seam": [], "interior": []}
    hits = []
    for x, y in pts:
        p = shapely.Point(x, y)
        idx = [j for j in tree.query(p) if geoms[j].contains(p)]
        hits.append(idx[0] if idx else -1)
    counts = np.bincount([h for h in hits if h >= 0], minlength=len(geoms)) if geoms else []
    for (x, y), h in zip(pts, hits):
        near = min([abs(x - s) for s in seam_x] + [abs(y - s) for s in seam_y] + [1e9]) <= 150
        res["seam" if near else "interior"].append((h >= 0, h >= 0 and counts[h] == 1))
    out = {}
    for k, v in res.items():
        v = np.array(v) if v else np.zeros((0, 2), bool)
        out[f"gt_{k}_n"] = len(v)
        out[f"gt_{k}_in_cell"] = round(100 * v[:, 0].mean(), 1) if len(v) else ""
        out[f"gt_{k}_alone"] = round(100 * v[:, 1].mean(), 1) if len(v) else ""
    return out


def strip_agreement(cells):
    stats = []
    for i in range(len(patches)):
        for j in range(i + 1, len(patches)):
            a, b = shapely.box(*patches[i]), shapely.box(*patches[j])
            ov = a.intersection(b)
            if ov.area == 0:
                continue
            core = ov.buffer(-40)
            if core.is_empty:
                continue
            ci = [g for g in cells[cells.patch == i].geometry if core.contains(g.centroid)]
            cj = [g for g in cells[cells.patch == j].geometry if ov.intersects(g)]
            if not ci or not cj:
                continue
            tree = STRtree(cj)
            for g in ci:
                best = max([g.intersection(cj[k]).area / g.union(cj[k]).area for k in tree.query(g)], default=0)
                stats.append(best)
    s = np.array(stats)
    return {"strip_cells": len(s), "strip_same_pct": round(100 * (s >= 0.8).mean(), 1),
            "strip_differs_pct": round(100 * ((s >= 0.5) & (s < 0.8)).mean(), 1),
            "strip_conflict_pct": round(100 * (s < 0.5).mean(), 1)}


rows = []
pts = gt_points()
for norm in ["patch", "global"]:
    if not (D / norm).exists():
        continue
    cells = load(norm)
    agree = strip_agreement(cells)
    for method in ["sopa", "owner", "match_context", "match_cellprob"]:
        if method == "sopa":
            res = solve_conflicts(cells, threshold=0.5, patch_indices=cells.patch.values)
        elif method == "owner":
            res = owner(cells)
        else:
            res = match(cells, method.split("_")[1])
        geoms = [g for g in res.geometry if g is not None and not g.is_empty]
        total, on_edge, with_ = straight_stats(geoms)
        row = {"norm": norm, "stitch": method, "cells": len(geoms), "overlapping_pairs": overlaps(geoms),
               "straight_edges": total, "straight_on_patch_edge": on_edge,
               "cells_with_straight_edge_pct": round(100 * with_ / len(geoms), 2), **gt_scores(geoms, pts), **agree}
        res.to_parquet(D / f"resolved_{norm}_{method}.parquet")
        rows.append(row)
        print(row, flush=True)
if (D / "whole" / "whole_region.parquet").exists():   # reference: one Cellpose call, no patches
    w = gpd.read_parquet(D / "whole" / "whole_region.parquet")
    geoms = [g for g in w.geometry if g is not None and not g.is_empty]
    total, on_edge, with_ = straight_stats(geoms)
    row = {"norm": "whole region", "stitch": "none", "cells": len(geoms), "overlapping_pairs": overlaps(geoms),
           "straight_edges": total, "straight_on_patch_edge": on_edge,
           "cells_with_straight_edge_pct": round(100 * with_ / len(geoms), 2), **gt_scores(geoms, pts),
           **{k: "" for k in ["strip_cells", "strip_same_pct", "strip_differs_pct", "strip_conflict_pct"]}}
    rows.append(row)
    print(row, flush=True)
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
