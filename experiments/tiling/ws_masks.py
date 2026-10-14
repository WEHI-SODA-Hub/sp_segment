#!/usr/bin/env python3
"""
Cellpose mask building with a watershed instead of the 11 x 11 seed box.

Cellpose's get_masks_torch grows each seed with a 3x3 max-pool, 5 times, inside
an 11 x 11 box in the space of final pixel positions, through bins with > 2
pixels. Pixels whose final position lies outside the box get no mask, and where
two boxes meet the higher label wins, so cells get straight edges.

get_masks_ws uses the same histogram and the same seeds, then labels the bins
with > 2 pixels by a watershed on the negative histogram from the seeds, with no
box limit. Everything else is cellpose.dynamics.compute_masks.

Usage (inside the sopa-cellpose container):
  python ws_masks.py <flows.npz with dP, cellprob> <out prefix>
writes <out>_cellpose.tiff (unchanged Cellpose) and <out>_ws.tiff.
"""
import sys, numpy as np, tifffile, torch, fastremap
from scipy import ndimage as ndi
from skimage.segmentation import watershed, relabel_sequential
from cellpose import dynamics, utils


def get_masks_ws(pt, inds, shape0, rpad=20, max_size_fraction=0.4):
    pt = pt.cpu().numpy().astype(np.int64) + rpad
    shape = tuple(np.array(shape0) + 2 * rpad)
    for i in range(2):
        pt[i] = np.clip(pt[i], 0, shape0[i] + rpad - 1)
    h = np.zeros(shape, np.int32)
    np.add.at(h, (pt[0], pt[1]), 1)
    hmax = ndi.maximum_filter(h, size=5, mode="constant")
    seeds = np.argwhere((h - hmax > -1e-6) & (h > 10))          # same seed rule as Cellpose
    if len(seeds) == 0:
        return np.zeros(shape0, "uint16")
    markers = np.zeros(shape, np.int32)
    markers[seeds[:, 0], seeds[:, 1]] = np.arange(1, len(seeds) + 1)
    support = (h > 2) | (markers > 0)                           # same growth condition, no 11 x 11 box
    L = watershed(-h.astype(np.float32), markers=markers, mask=support, connectivity=2)
    M0 = np.zeros(shape0, np.int32)
    M0[inds] = L[pt[0], pt[1]]
    uniq, counts = fastremap.unique(M0, return_counts=True)
    big = uniq[(counts > np.prod(shape0) * max_size_fraction) & (uniq > 0)]
    if len(big):
        M0[np.isin(M0, big)] = 0
    return relabel_sequential(M0)[0]


def compute(dP, cellprob, ws, niter=200, cellprob_threshold=0.0, flow_threshold=0.4, min_size=15,
            device=torch.device("cpu")):
    inds = np.nonzero(cellprob > cellprob_threshold)
    p = dynamics.follow_flows(dP * (cellprob > cellprob_threshold) / 5., inds=inds, niter=niter, device=device)
    p = p.int() if torch.is_tensor(p) else torch.from_numpy(p).int()
    if ws:
        mask = get_masks_ws(p, inds, dP.shape[1:])
    else:
        mask = dynamics.get_masks_torch(p, inds, dP.shape[1:], max_size_fraction=0.4)
    mask = mask.astype(np.int32)
    if mask.max() > 0:
        mask = dynamics.remove_bad_flow_masks(mask, dP, threshold=flow_threshold, device=device)
    mask = utils.fill_holes_and_remove_small_masks(mask, min_size=min_size)
    ids, cnt = np.unique(mask, return_counts=True)
    mask[np.isin(mask, ids[(cnt < 100) & (ids > 0)])] = 0       # pipeline cellpose_min_area
    return relabel_sequential(mask)[0]


if __name__ == "__main__":
    f, out = sys.argv[1], sys.argv[2]
    z = np.load(f)
    dP, cp = z["dP"].astype(np.float32), z["cellprob"].astype(np.float32)
    dev = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
    for name, ws in [("cellpose", False), ("ws", True)]:
        m = compute(dP, cp, ws, device=dev)
        tifffile.imwrite(f"{out}_{name}.tiff", m.astype(np.int32))
        print(out, name, int(m.max()), "cells", flush=True)
