"""Raw-max cells that norm_mean / norm_max lack (best IoU < 0.3), split by DAPI and by what the
other input has at that place.
  clear nucleus: cell's DAPI p90 >= 50% of the median production cell's (per crop)
  empty: < 20% of the cell is covered by the other input's cells
  merged: one other cell covers >= 50% of it and is >= 1.5x its area
  other: anything else (partial overlap / shifted)
python missed_cells.py  (writes out/missed_vs_raw_max.csv)"""
import csv
import numpy as np, tifffile
from scipy import ndimage as ndi
T = '/vast/scratch/users/mckay.m/sp_segment_tests/'
R = T + 'membrane_combine_all/results/cellmeasurement/'
rows = []
for crop in ['b6_immune', 'b6_macro', 'b2_endo', 'b2_vsmc', 'b2_tcell']:
    dapi = tifffile.imread(T + f'membrane_combine_csl/inputs/{crop}_source.ome.tif')[0].astype(np.float32)
    prod = np.load(T + f'membrane_artefacts/out/{crop}_production_mask.npy')
    def p90(mask, labs):
        return np.array(ndi.labeled_comprehension(dapi, mask, labs, lambda v: np.percentile(v, 90), float, 0))
    pl = np.unique(prod); pl = pl[pl > 0]; ref = np.median(p90(prod, pl))
    raw = tifffile.imread(R + f'{crop}_raw_max_cellpose_mask.tiff').astype(np.int64)
    rl = np.unique(raw); rl = rl[rl > 0]
    rarea = np.bincount(raw.ravel())
    rdapi = dict(zip(rl, p90(raw, rl)))
    for other in ['norm_mean', 'norm_max', 'single']:
        o = tifffile.imread(R + f'{crop}_{other}_cellpose_mask.tiff').astype(np.int64)
        oarea = np.bincount(o.ravel())
        both = (raw > 0) & (o > 0)
        pairs, inter = np.unique(np.stack([raw[both], o[both]]), axis=1, return_counts=True)
        best, cover, top, topb = {}, {}, {}, {}
        for (a, b), k in zip(pairs.T, inter):
            iou = k / (rarea[a] + oarea[b] - k)
            best[a] = max(best.get(a, 0), iou); cover[a] = cover.get(a, 0) + k
            if k > top.get(a, 0): top[a], topb[a] = k, b
        cnt = {'missing': 0, 'clear_nucleus': 0, 'empty': 0, 'merged': 0, 'other': 0}
        for a in rl:
            if best.get(a, 0) >= 0.3:
                continue
            cnt['missing'] += 1
            if rdapi[a] < 0.5 * ref:
                continue
            cnt['clear_nucleus'] += 1
            if cover.get(a, 0) < 0.2 * rarea[a]:
                cnt['empty'] += 1
            elif top[a] >= 0.5 * rarea[a] and oarea[topb[a]] >= 1.5 * rarea[a]:
                cnt['merged'] += 1
            else:
                cnt['other'] += 1
        rows.append({'crop': crop, 'vs': other, 'raw_max_cells': len(rl), 'other_cells': int((oarea[1:] > 0).sum()),
                     'raw_max_not_in_other': cnt['missing'], 'with_clear_nucleus': cnt['clear_nucleus'],
                     'other_empty_there': cnt['empty'], 'other_merged_into_bigger_cell': cnt['merged'],
                     'other_partial_or_shifted': cnt['other']})
with open(T + 'membrane_artefacts/out/missed_vs_raw_max.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
