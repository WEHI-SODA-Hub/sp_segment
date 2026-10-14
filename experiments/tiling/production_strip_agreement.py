"""For one production image, compare the two patches' Cellpose outlines of each cell in the shared strips (IoU)."""
import pyarrow.parquet as pq, shapely, numpy as np, glob, sys
from shapely.strtree import STRtree
z=sys.argv[1]
def geoms(p):
    t=pq.read_table(p); col=[c for c in t.column_names if "geom" in c.lower()][0]
    g=shapely.from_wkb(t.column(col).to_pylist()); return list(g)
patches=geoms(glob.glob(f"{z}/shapes/image_patches/*.parquet")[0])
cells={i:geoms(f"{z}/.sopa_cache/cellpose_boundaries/{i}.parquet") for i in range(len(patches))}
n_pairs=0; stats=[]
for i in range(len(patches)):
    for j in range(i+1,len(patches)):
        ov=patches[i].intersection(patches[j])
        if ov.area==0: continue
        n_pairs+=1
        core=ov.buffer(-40)                              # cells well inside the shared strip
        if core.is_empty: continue
        a=[c for c in cells[i] if core.contains(c.centroid)]
        bj=[c for c in cells[j] if ov.intersects(c)]
        if not a or not bj: continue
        tree=STRtree(bj)
        for c in a:
            cand=[bj[k] for k in tree.query(c)]
            best=max([c.intersection(d).area/c.union(d).area for d in cand], default=0)
            stats.append(best)
s=np.array(stats)
print(f"{len(patches)} patches, {n_pairs} overlapping pairs, {len(s)} cells in shared strips")
print(f"  same cell in both patches (IoU>=0.8): {100*(s>=0.8).mean():.1f}%")
print(f"  shape differs (0.5<=IoU<0.8):        {100*((s>=0.5)&(s<0.8)).mean():.1f}%")
print(f"  split/merged/missing (IoU<0.5):      {100*(s<0.5).mean():.1f}%")
