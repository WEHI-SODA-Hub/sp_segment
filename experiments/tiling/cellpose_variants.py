"""Cellpose variants on one crop: baseline, niter 1000/2000, augment, rotate 45 deg."""
import sys, numpy as np, tifffile
from cellpose import models
from scipy import ndimage as ndi
inp, out = sys.argv[1], sys.argv[2]
img = tifffile.imread(inp)
x = np.stack([img[1], img[0]], axis=-1).astype(np.float32)          # membrane, DAPI
model = models.CellposeModel(gpu=True, pretrained_model="cpsam_v2")
kw = dict(diameter=30, channels=[1, 2], flow_threshold=0.4, cellprob_threshold=0.0)
def save(name, mask): tifffile.imwrite(f"{out}_{name}.tiff", mask.astype(np.int32)); print(out, name, int(mask.max()), flush=True)
mask, flows, _ = model.eval(x, **kw)
save("base", mask)
np.savez_compressed(f"{out}_base_flows.npz", dP=flows[1].astype(np.float16), cellprob=flows[2].astype(np.float16))
for n in [1000, 2000]:
    save(f"niter{n}", model.eval(x, niter=n, **kw)[0])
save("augment", model.eval(x, augment=True, **kw)[0])
# rotate 45 degrees (pad so nothing is cut), segment, rotate the mask back
H, W = x.shape[:2]; pad = int(0.21 * max(H, W)) + 8
xp = np.pad(x, ((pad, pad), (pad, pad), (0, 0)))
xr = np.stack([ndi.rotate(xp[..., c], 45, reshape=False, order=1) for c in range(2)], -1)
mr = model.eval(xr, **kw)[0]
mb = ndi.rotate(mr, -45, reshape=False, order=0)[pad:pad + H, pad:pad + W]
save("rot45", mb)
