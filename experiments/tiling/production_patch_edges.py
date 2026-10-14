"""Count straight polygon edges (>= 15 px) in production Cellpose outlines that lie on a sopa patch edge, against chance."""
import pyarrow.parquet as pq, shapely, numpy as np, glob, sys
B="/vast/projects/CSL_COMET/csl_artherosclerosis"
zs=sorted(glob.glob(f"{B}/seg_batch*/sopa/*.zarr"))[::8][:8]
def geoms(p):
    t=pq.read_table(p); col=[c for c in t.column_names if "geom" in c.lower()][0]
    return shapely.from_wkb(t.column(col).to_pylist())
tot=[0,0,0]
for z in zs:
    patches=geoms(glob.glob(f"{z}/shapes/image_patches/*.parquet")[0])
    b=np.array([p.bounds for p in patches])
    xs=np.unique(np.concatenate([b[:,0],b[:,2]])); ys=np.unique(np.concatenate([b[:,1],b[:,3]]))
    W,H=b[:,2].max(),b[:,3].max()
    cells=geoms(glob.glob(f"{z}/shapes/cellpose_boundaries/*.parquet")[0])
    n=near=0
    for poly in cells:
        if poly.geom_type!="Polygon": poly=max(poly.geoms,key=lambda q:q.area)
        r=np.asarray(poly.exterior.coords); s=np.diff(r,axis=0); L=np.hypot(*s.T)
        for k in np.where(L>=15)[0]:
            if abs(s[k,1])<1e-6:
                y=r[k,1]; 
                if 2<y<H-2: n+=1; near+= np.min(np.abs(ys-y))<=1.5
            elif abs(s[k,0])<1e-6:
                x=r[k,0]
                if 2<x<W-2: n+=1; near+= np.min(np.abs(xs-x))<=1.5
    chance=n*(3*(len(ys)-2)/H+3*(len(xs)-2)/W)/2
    tot[0]+=n; tot[1]+=near; tot[2]+=chance
    print(f"{z.split('/')[-1][:60]}: {len(patches)} patches, {len(cells)} cells, straight segs {n}, on patch edges {near} (chance {chance:.1f})")
print("total", tot[0], "on patch edges", tot[1], "chance %.1f"%tot[2])
