#!/usr/bin/env python3
"""
Segment a 4 x 4 block of production-sized patches (3000 px, 250 px overlap) with
Cellpose as the pipeline does, normalising either per patch (current) or once
for the whole region ("global": 1st/99th percentile per channel over the region).

Writes <out>/<norm>/patch_<k>.parquet: polygons in image coordinates with
columns patch, mean_cellprob. Run inside the sopa-cellpose container.
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np, tifffile, zarr, geopandas as gpd
from scipy import ndimage as ndi
from cellpose import models
from sopa.shapes import vectorize

ap = argparse.ArgumentParser()
ap.add_argument("--tiff"); ap.add_argument("--x0", type=int); ap.add_argument("--y0", type=int)
ap.add_argument("--n", type=int, default=4); ap.add_argument("--width", type=int, default=3000)
ap.add_argument("--overlap", type=int, default=250); ap.add_argument("--out")
ap.add_argument("--norm", nargs="+", default=["patch", "global"])
args = ap.parse_args()
sys.path.insert(0, "/vast/scratch/users/mckay.m/sp_segment/bin")
from tiff_channels import get_channel_names
MEMBRANE = ["R2_CSL_CD11c", "R2_CSL_VECadherin", "R1_CSL_aSMA", "R2_CSL_Adipophilin", "R1_CSL_S100A9",
            "R1_CSL_FAP", "R1_CSL_MMP9", "R1_CSL_CD45", "R1_CSL_MPO", "R1_CSL_HistoneH3", "R1_CSL_CD3", "R2_CSL_CD68"]
names = get_channel_names(args.tiff)
step = args.width - args.overlap
store = tifffile.imread(args.tiff, aszarr=True)
g = zarr.open(store, mode="r")
arr = g["0"] if not hasattr(g, "shape") else g
H, W = arr.shape[-2:]
X1 = min(args.x0 + step * (args.n - 1) + args.width, W); Y1 = min(args.y0 + step * (args.n - 1) + args.width, H)
t = time.time()
dapi = np.asarray(arr[names.index("R1_DAPI"), args.y0:Y1, args.x0:X1])
mem = np.zeros_like(dapi)
for m in MEMBRANE:                                   # raw per-pixel max, as COMBINECHANNELS
    np.maximum(mem, np.asarray(arr[names.index(m), args.y0:Y1, args.x0:X1]), out=mem)
print(f"read region {dapi.shape} in {time.time()-t:.0f}s", flush=True)
img = np.stack([mem, dapi], -1).astype(np.float32)   # membrane, DAPI (pipeline channel order)
lowhigh = [np.percentile(img[..., c][::4, ::4], [1, 99]).tolist() for c in range(2)]
meta = {"x0": args.x0, "y0": args.y0, "x1": X1, "y1": Y1, "width": args.width, "overlap": args.overlap,
        "n": args.n, "global_lowhigh": lowhigh, "patches": []}
model = models.CellposeModel(gpu=True, pretrained_model="cpsam_v2")
c0 = img[1000:1512, 1000:1512]                       # sanity check before the full run
# model.eval normalises its input in place, so every call gets a copy
n_p = model.eval(c0.copy(), diameter=30, normalize=True)[0].max()
n_g = model.eval(c0.copy(), diameter=30, normalize={"normalize": True, "lowhigh": lowhigh})[0].max()
print(f"sanity check 512 px: per-patch {n_p} cells, global {n_g} cells", flush=True)
if n_p > 0 and n_g == 0:
    sys.exit("global normalisation gives no cells; stopping")
for norm in args.norm:
    out = Path(args.out) / norm; out.mkdir(parents=True, exist_ok=True)
    k = 0
    for j in range(args.n):
        for i in range(args.n):
            px0, py0 = i * step, j * step
            px1, py1 = min(px0 + args.width, img.shape[1]), min(py0 + args.width, img.shape[0])
            if norm == "patch":
                meta["patches"].append([args.x0 + px0, args.y0 + py0, args.x0 + px1, args.y0 + py1])
            patch = img[py0:py1, px0:px1].copy()        # eval normalises in place
            normalize = True if norm == "patch" else {"normalize": True, "lowhigh": lowhigh}
            t = time.time()
            mask, flows, _ = model.eval(patch, diameter=30, flow_threshold=0.4, cellprob_threshold=0.0,
                                        normalize=normalize)
            ids = np.arange(1, mask.max() + 1)
            cp = ndi.mean(flows[2], mask, ids) if len(ids) else []
            cells = vectorize(mask) if len(ids) else gpd.GeoDataFrame(geometry=[])
            # score each polygon by the mean cell probability of the mask label under it
            labs = [int(mask[min(int(pt.y), mask.shape[0] - 1), min(int(pt.x), mask.shape[1] - 1)])
                    for pt in cells.geometry.representative_point()] if len(ids) else []
            cells["mean_cellprob"] = [float(cp[l - 1]) if l > 0 else np.nan for l in labs]
            cells = cells[cells.area >= 100].copy()
            cells["geometry"] = cells.geometry.translate(args.x0 + px0, args.y0 + py0)
            cells["patch"] = k
            cells.reset_index(drop=True).to_parquet(out / f"patch_{k}.parquet")
            print(f"{norm} patch {k}: {len(cells)} cells, {time.time()-t:.0f}s", flush=True)
            k += 1
(Path(args.out) / "meta.json").write_text(json.dumps(meta, indent=1))
print("DONE")
