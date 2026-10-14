#!/usr/bin/env python3
"""
Rotation consensus for Cellpose: run at 0 and 45 degrees, rotate the 45-degree
flows and cell probability back, average both runs, and build masks once with
cellpose.dynamics.compute_masks. Writes <out>_base.tiff, <out>_consensus.tiff.
Masks smaller than 100 px are removed (pipeline cellpose_min_area).
Run inside the sopa-cellpose container on a GPU.
"""
import sys, numpy as np, tifffile, torch
from scipy import ndimage as ndi
from cellpose import models, dynamics

inp, out = sys.argv[1], sys.argv[2]
img = tifffile.imread(inp)
x = np.stack([img[1], img[0]], -1).astype(np.float32)          # membrane, DAPI
H, W = x.shape[:2]
model = models.CellposeModel(gpu=True, pretrained_model="cpsam_v2")
kw = dict(diameter=30, flow_threshold=0.4, cellprob_threshold=0.0)
dev = torch.device("cuda")


def min_area(mask, n=100):
    ids, counts = np.unique(mask, return_counts=True)
    small = ids[(counts < n) & (ids > 0)]
    mask = mask.copy(); mask[np.isin(mask, small)] = 0
    from skimage.segmentation import relabel_sequential
    return relabel_sequential(mask)[0]


mask0, flows0, _ = model.eval(x, **kw)
dP0, cp0 = flows0[1].astype(np.float32), flows0[2].astype(np.float32)

pad = int(0.21 * max(H, W)) + 8
xp = np.pad(x, ((pad, pad), (pad, pad), (0, 0)))
xr = np.stack([ndi.rotate(xp[..., c], 45, reshape=False, order=1) for c in range(2)], -1)
_, flowsr, _ = model.eval(xr, **kw)
dPr, cpr = flowsr[1].astype(np.float32), flowsr[2].astype(np.float32)


def back(a, order=1):
    return ndi.rotate(a, -45, reshape=False, order=order, cval=0)[pad:pad + H, pad:pad + W]


cpb = back(cpr)
dy, dx = back(dPr[0]), back(dPr[1])
c, s = np.cos(np.pi / 4), np.sin(np.pi / 4)
inside = cp0 > 0
best = None
for sign in (1, -1):                     # pick the vector rotation that matches the 0-degree flows
    vy = c * dy - sign * s * dx
    vx = sign * s * dy + c * dx
    corr = np.corrcoef(np.concatenate([vy[inside], vx[inside]]),
                       np.concatenate([dP0[0][inside], dP0[1][inside]]))[0, 1]
    if best is None or corr > best[0]:
        best = (corr, sign, np.stack([vy, vx]))
print(f"{inp}: vector rotation sign {best[1]}, flow correlation with 0-degree run {best[0]:.3f}", flush=True)
dPb = best[2]

def masks(dP, cp):
    res = dynamics.compute_masks(dP, cp, niter=200, cellprob_threshold=0.0, flow_threshold=0.4, device=dev)
    return res[0] if isinstance(res, tuple) else res


base = masks(dP0, cp0)                    # same maps as the 0-degree eval, rebuilt with the same call
cons = masks((dP0 + dPb) / 2, (cp0 + cpb) / 2)
base, cons = min_area(np.asarray(base)), min_area(np.asarray(cons))
print(f"{inp}: eval mask {int(mask0.max())} cells, recomputed base {int(base.max())}, consensus {int(cons.max())}", flush=True)
tifffile.imwrite(f"{out}_base.tiff", base.astype(np.int32))
tifffile.imwrite(f"{out}_consensus.tiff", cons.astype(np.int32))
