"""Count straight (>= 20 px, horizontal or vertical) boundaries between two cells in each Cellpose and CellSAM mask."""
import numpy as np, tifffile
from pathlib import Path
T=Path("/vast/scratch/users/mckay.m/sp_segment_tests")
MIN=20
def cc_runs(mask):
    out=[]
    for axis in (0,1):
        m = mask if axis==0 else mask.T
        a,b=m[:-1],m[1:]
        e=(a!=b)&(a>0)&(b>0)
        for y in np.where(e.sum(1)>=MIN)[0]:
            d=np.diff(np.concatenate([[0],e[y].astype(np.int8),[0]]))
            out += [(axis,y+1,s) for s,t in zip(np.where(d==1)[0],np.where(d==-1)[0]) if t-s>=MIN]
    return out
tot={}
for d in [T/"membrane_combine_all/results/cellmeasurement", T/"membrane_combine_stromal/results/cellmeasurement"]:
    for f in sorted(d.glob("*_mask.tiff")):
        seg=f.name.split("_")[-2]
        if seg=="mesmer": continue
        r=cc_runs(tifffile.imread(f))
        on600=sum(1 for _,p,_ in r if min(abs(p-600),abs(p-1200),abs(p-1800))<=2)
        t=tot.setdefault(seg,[0,0,0]); t[0]+=1; t[1]+=len(r); t[2]+=on600
        if seg=="cellpose" and r: print(f.name, len(r), r[:4])
for seg,(n,k,s) in tot.items(): print(f"{seg}: {n} masks, straight cell|cell runs >= {MIN}px: {k} ({s} on the 600 px CellSAM grid)")
