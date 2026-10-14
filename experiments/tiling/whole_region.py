"""Reference: Cellpose on the whole 4x4-patch region in one call (no sopa patches)."""
import sys, json, time, numpy as np, tifffile, zarr, geopandas as gpd
from pathlib import Path
from cellpose import models
from sopa.shapes import vectorize
sys.path.insert(0, "/vast/scratch/users/mckay.m/sp_segment/bin")
from tiff_channels import get_channel_names
D = Path("/vast/scratch/users/mckay.m/sp_segment_tests/tiling_test")
meta = json.loads((D / "patches/meta.json").read_text())
tif = sys.argv[1]
MEMBRANE = ["R2_CSL_CD11c", "R2_CSL_VECadherin", "R1_CSL_aSMA", "R2_CSL_Adipophilin", "R1_CSL_S100A9",
            "R1_CSL_FAP", "R1_CSL_MMP9", "R1_CSL_CD45", "R1_CSL_MPO", "R1_CSL_HistoneH3", "R1_CSL_CD3", "R2_CSL_CD68"]
names = get_channel_names(tif); g = zarr.open(tifffile.imread(tif, aszarr=True), mode="r"); arr = g["0"] if not hasattr(g, "shape") else g
y0, y1, x0, x1 = meta["y0"], meta["y1"], meta["x0"], meta["x1"]
dapi = np.asarray(arr[names.index("R1_DAPI"), y0:y1, x0:x1]); mem = np.zeros_like(dapi)
for m in MEMBRANE: np.maximum(mem, np.asarray(arr[names.index(m), y0:y1, x0:x1]), out=mem)
img = np.stack([mem, dapi], -1).astype(np.float32)
t = time.time()
mask, flows, _ = models.CellposeModel(gpu=True, pretrained_model="cpsam_v2").eval(img, diameter=30, flow_threshold=0.4, cellprob_threshold=0.0)
print(f"whole region {img.shape[:2]}: {mask.max()} cells in {time.time()-t:.0f}s", flush=True)
cells = vectorize(mask); cells = cells[cells.area >= 100].copy()
cells["geometry"] = cells.geometry.translate(x0, y0); cells["patch"] = 0; cells["mean_cellprob"] = 0.0
out = D / "patches/whole"; out.mkdir(exist_ok=True); cells.reset_index(drop=True).to_parquet(out / "whole_region.parquet")
print("DONE", len(cells))
