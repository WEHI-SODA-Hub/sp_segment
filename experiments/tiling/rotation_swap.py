"""Prototype: replace base cells that have a straight edge with the matching cell from the 45-degree run."""
import numpy as np, tifffile
from pathlib import Path
from scipy import ndimage as ndi
D=Path("/vast/scratch/users/mckay.m/sp_segment_tests/cellpose_tile_overlap")
MIN=15
def straight_cells(mask):
    """labels of cells with an axis-aligned straight edge of >= MIN px"""
    bad=set()
    for m in (mask, mask.T):
        a,b=m[:-1],m[1:]; e=(a!=b)&((a>0)|(b>0)); e[:2]=False; e[-2:]=False
        for y in np.where(e.sum(1)>=MIN)[0]:
            d=np.diff(np.concatenate([[0],e[y].astype(np.int8),[0]]))
            for s,t in zip(np.where(d==1)[0],np.where(d==-1)[0]):
                if t-s>=MIN:
                    for lab in np.unique(np.concatenate([a[y,s:t],b[y,s:t]])):
                        if lab: bad.add(int(lab))
    return bad
for r in ["b6_immune","b6_macro","b2_vsmc"]:
    base=tifffile.imread(D/f"v_{r}_base.tiff"); rot=tifffile.imread(D/f"v_{r}_rot45.tiff")
    bad=straight_cells(base); rot_bad=straight_cells(rot)
    out=base.copy(); nxt=out.max()+1; replaced=kept=0
    objs=ndi.find_objects(base)
    for lab in sorted(bad):
        sl=objs[lab-1]
        if sl is None: continue
        sl=tuple(slice(max(s.start-20,0),s.stop+20) for s in sl)
        cell=base[sl]==lab; cand=rot[sl][cell]; cand=cand[cand>0]
        if not cand.size: kept+=1; continue
        j=np.bincount(cand).argmax(); other=rot[sl]==j
        iou=(cell&other).sum()/(cell|other).sum()
        if iou<0.5 or j in rot_bad: kept+=1; continue
        region=out[sl]
        region[cell]=0                                  # remove the straight-edged version
        region[other&(region==0)]=nxt; nxt+=1           # paste the 45-degree version where free
        replaced+=1
    tifffile.imwrite(D/f"v_{r}_consensus.tiff", out)
    print(f"{r}: {len(np.unique(base))-1} cells, {len(bad)} with a straight edge -> replaced {replaced}, kept {kept}; straight cells after: {len(straight_cells(out))}")
