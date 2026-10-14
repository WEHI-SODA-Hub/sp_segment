"""Cells whose 90th-percentile DAPI is < 20% of the median production cell's (per crop).
python low_dapi.py  (writes out/low_dapi.csv)"""
import csv
import numpy as np, tifffile
from scipy import ndimage as ndi
T = '/vast/scratch/users/mckay.m/sp_segment_tests/'
STROM = {'b6_immune': 'fap', 'b6_macro': 'fap', 'b2_vsmc': 'asma'}
rows = []
for crop in ['b6_immune', 'b6_macro', 'b2_endo', 'b2_vsmc', 'b2_tcell']:
    dapi = tifffile.imread(T + f'membrane_combine_csl/inputs/{crop}_source.ome.tif')[0].astype(np.float32)
    masks = {m: T + f'membrane_combine_all/results/cellmeasurement/{crop}_{m}_cellpose_mask.tiff'
             for m in ['raw_max', 'norm_max', 'norm_mean', 'single']}
    if crop in STROM:
        t = STROM[crop]
        for m, lab in [('raw_max', 'FAP+aSMA max'), ('norm_max', 'FAP+aSMA rescaled max'),
                       ('single', 'FAP alone' if t == 'fap' else 'aSMA alone')]:
            masks[lab] = T + f'membrane_combine_stromal/results/cellmeasurement/{crop}_{t}_{m}_cellpose_mask.tiff'
    def p90(mask):
        labs = np.unique(mask); labs = labs[labs > 0]
        return np.array(ndi.labeled_comprehension(dapi, mask, labs, lambda v: np.percentile(v, 90), float, 0))
    prod = np.load(T + f'membrane_artefacts/out/{crop}_production_mask.npy')
    pv = p90(prod); ref = np.median(pv)
    for name, path in list(masks.items()) + [('production', None)]:
        v = pv if path is None else p90(tifffile.imread(path))
        rows.append({'crop': crop, 'input': name, 'cells': len(v), 'low_dapi': int((v < 0.2 * ref).sum()),
                     'low_dapi_pct': round(100 * float((v < 0.2 * ref).mean()), 1), 'ref_p90': round(float(ref))})
with open(T + 'membrane_artefacts/out/low_dapi.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
tot = {}
for r in rows:
    t = tot.setdefault(r['input'], [0, 0, 0]); t[0] += r['cells']; t[1] += r['low_dapi']; t[2] += 1
for k, (c, l, n) in tot.items():
    print(f'{k:24s} crops={n} cells={c:5d} low_dapi={l:4d} ({100*l/c:.1f}%)')
print(); print(' '.join(f"{r['crop']}/{r['input']}={r['low_dapi']}" for r in rows))
