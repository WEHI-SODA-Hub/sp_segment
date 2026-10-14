"""Second round of Cellpose-SAM variants on the 5 CSL crops (sopa defaults, CLAHE, 3-channel input).

Variants (raw-max membrane input unless named otherwise):
  f2_cp0, f2_cp-2, f2_cp-4, f2_cp-6   flow threshold 2 with cell-probability threshold 0/-2/-4/-6
  clahe                                sopa preprocessing: Gaussian sigma 1, then CLAHE clip 0.2 per channel
  clahe_f2_cp-2                        the same + flow 2, cellprob -2
  sopa_default                         the same + flow 2, cellprob -6 (sopa's own defaults)
  3ch, 3ch_f2_cp-2, 3ch_f2_cp-4        3 channels: max of immune markers, max of stromal/vascular markers, DAPI
  3ch_clahe_f2_cp-2                    3 channels with sopa preprocessing, flow 2, cellprob -2
  norm_mean_f2_cp-2                    norm_mean input, flow 2, cellprob -2
Diameter 30. Cells < 100 px removed (the pipeline's cellpose_min_area).
Writes masks/<crop>_<variant>.tiff
python segment_variants2.py
"""
import json

import numpy as np
import tifffile
from cellpose import models
from scipy.ndimage import gaussian_filter
from skimage import exposure

T = "/vast/scratch/users/mckay.m/sp_segment_tests/"
I = T + "membrane_combine_csl/inputs/"
O = T + "cell_shape/masks/"
CFG = json.load(open("/vast/scratch/users/mckay.m/sp_segment/experiments/membrane_combine/config_csl.json"))
MARKERS = CFG["markers"]
IMMUNE = ["R2_CSL_CD11c", "R1_CSL_S100A9", "R1_CSL_MMP9", "R1_CSL_CD45", "R1_CSL_MPO",
          "R1_CSL_HistoneH3", "R1_CSL_CD3", "R2_CSL_CD68", "R2_CSL_Adipophilin"]
STROMAL = ["R2_CSL_VECadherin", "R1_CSL_aSMA", "R1_CSL_FAP"]
assert sorted(IMMUNE + STROMAL) == sorted(MARKERS)

model = models.CellposeModel(gpu=True, pretrained_model="cpsam_v2")


def sopa_preprocess(chans):
    """sopa StainingSegmentation: per channel gaussian_filter(sigma=1) on the integer image,
    then skimage equalize_adapthist(clip_limit=0.2, kernel_size=None)."""
    out = []
    for c in chans:
        c = gaussian_filter(c, sigma=1)
        out.append(exposure.equalize_adapthist(c, clip_limit=0.2).astype(np.float32))
    return out


def clean(mask, min_px=100):
    area = np.bincount(mask.ravel())
    small = np.flatnonzero(area < min_px)
    small = small[small > 0]
    mask = np.where(np.isin(mask, small), 0, mask)
    _, inv = np.unique(mask, return_inverse=True)
    return inv.reshape(mask.shape).astype(np.int32)


def run(name, chans, flow=0.4, cellprob=0.0):
    x = np.stack([np.asarray(c, np.float32) for c in chans], axis=-1)
    m = model.eval(x.copy(), channel_axis=2, diameter=30, flow_threshold=flow,
                   cellprob_threshold=cellprob)[0]
    m = clean(m)
    tifffile.imwrite(O + name + ".tiff", m)
    print(name, int(m.max()), flush=True)


for crop in ["b6_immune", "b6_macro", "b2_endo", "b2_vsmc", "b2_tcell"]:
    src = tifffile.imread(I + f"{crop}_source.ome.tif")          # DAPI + 12 markers, uint16
    dapi = src[0]
    raw = tifffile.imread(I + f"{crop}_raw_max.ome.tif")[1]
    nm = tifffile.imread(I + f"{crop}_norm_mean.ome.tif")[1]
    immune = src[[1 + MARKERS.index(m) for m in IMMUNE]].max(0)
    stromal = src[[1 + MARKERS.index(m) for m in STROMAL]].max(0)
    run(f"{crop}_raw_max_base_check", [raw, dapi])
    for cp in (0, -2, -4, -6):
        run(f"{crop}_raw_max_f2_cp{cp}", [raw, dapi], flow=2, cellprob=cp)
    pre = sopa_preprocess([raw, dapi])
    run(f"{crop}_raw_max_clahe", pre)
    run(f"{crop}_raw_max_clahe_f2_cp-2", pre, flow=2, cellprob=-2)
    run(f"{crop}_raw_max_sopa_default", pre, flow=2, cellprob=-6)
    run(f"{crop}_3ch", [immune, stromal, dapi])
    run(f"{crop}_3ch_f2_cp-2", [immune, stromal, dapi], flow=2, cellprob=-2)
    run(f"{crop}_3ch_f2_cp-4", [immune, stromal, dapi], flow=2, cellprob=-4)
    run(f"{crop}_3ch_clahe_f2_cp-2", sopa_preprocess([immune, stromal, dapi]), flow=2, cellprob=-2)
    run(f"{crop}_norm_mean_f2_cp-2", [nm, dapi], flow=2, cellprob=-2)
