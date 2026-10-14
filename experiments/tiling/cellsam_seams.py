"""Count CellSAM cells cut in two or left as half cells on its 600 px block seams (rows/cols 600, 1200, 1800)."""
import numpy as np, tifffile, json, sys
from pathlib import Path
sys.path.insert(0, "/vast/scratch/users/mckay.m/csl_composite_check/membrane_test")
R=Path("/vast/scratch/users/mckay.m/sp_segment_tests/membrane_combine_all/results/cellmeasurement")
MIN=15
def classify(mask, seams=(600,1200,1800)):
    out={"cut_in_two":set(),"half_cell":set()}
    for axis in (0,1):
        m = mask if axis==0 else mask.T
        for s in seams:
            a, b = m[s-1], m[s]          # last row of one block, first of the next
            for lab_side, other in ((a,b),(b,a)):
                for lab in np.unique(lab_side[lab_side>0]):
                    on = lab_side==lab
                    if on.sum() < MIN: continue
                    across = other[on]
                    if (across==lab).mean() > 0.5: continue      # continues across: joined
                    if (across==0).mean() > 0.8: out["half_cell"].add(int(lab))
                    elif ((across>0)&(across!=lab)).mean() > 0.5: out["cut_in_two"].add(int(lab))
    return {k:len(v) for k,v in out.items()}
for seg in ["cellsam","cellpose"]:
    for reg in ["b6_immune","b6_macro","b2_tcell"]:
        mask=tifffile.imread(R/f"{reg}_raw_max_{seg}_mask.tiff")
        print(seg, reg, "cells", len(np.unique(mask))-1, classify(mask))
