"""Cellpose-SAM setting variants on the 5 CSL crops, for the cell-shape test.
For norm_mean and raw_max inputs: base (diameter 30, cellprob 0), diameter 40 and 50,
cellprob -1 and -2. Plus a DAPI-only Cellpose run (nuclei) per crop. Cells < 100 px removed
(the pipeline's cellpose_min_area). Writes masks/<crop>_<input>_<variant>.tiff
python segment_variants.py"""
import numpy as np, tifffile
from cellpose import models
T = '/vast/scratch/users/mckay.m/sp_segment_tests/'
I = T + 'membrane_combine_csl/inputs/'
O = T + 'cell_shape/masks/'
model = models.CellposeModel(gpu=True, pretrained_model="cpsam_v2")
VARIANTS = {"base": dict(diameter=30, cellprob_threshold=0.0),
            "d40": dict(diameter=40, cellprob_threshold=0.0),
            "d50": dict(diameter=50, cellprob_threshold=0.0),
            "cp-1": dict(diameter=30, cellprob_threshold=-1.0),
            "cp-2": dict(diameter=30, cellprob_threshold=-2.0)}

def clean(mask, min_px=100):
    area = np.bincount(mask.ravel())
    small = np.flatnonzero(area < min_px)
    small = small[small > 0]
    mask = np.where(np.isin(mask, small), 0, mask)
    _, inv = np.unique(mask, return_inverse=True)
    return inv.reshape(mask.shape).astype(np.int32)

def save(name, mask):
    mask = clean(mask)
    tifffile.imwrite(O + name + ".tiff", mask)
    print(name, int(mask.max()), flush=True)

for crop in ["b6_immune", "b6_macro", "b2_endo", "b2_vsmc", "b2_tcell"]:
    for inp in ["norm_mean", "raw_max"]:
        img = tifffile.imread(I + f"{crop}_{inp}.ome.tif")
        x = np.stack([img[1], img[0]], axis=-1).astype(np.float32)   # membrane, DAPI
        for v, kw in VARIANTS.items():
            m = model.eval(x.copy(), channels=[1, 2], flow_threshold=0.4, **kw)[0]
            save(f"{crop}_{inp}_{v}", m)
    img = tifffile.imread(I + f"{crop}_norm_mean.ome.tif")
    x = np.stack([img[0], np.zeros_like(img[0])], axis=-1).astype(np.float32)  # DAPI only
    m = model.eval(x.copy(), channels=[1, 0], diameter=30, flow_threshold=0.4, cellprob_threshold=0.0)[0]
    save(f"{crop}_nuclei", m)
