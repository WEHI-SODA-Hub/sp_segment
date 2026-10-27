#!/usr/bin/env python3
"""
Expand the cells in a whole-cell label mask by a fixed distance.

Each background pixel within --distance-px of a cell is given the label of the
nearest cell, so growth stops where two cells meet. This is the same rule as
skimage.segmentation.expand_labels.

With --nuclear-mask, pixels that belong to a nucleus but to no whole-cell label
are left as background. Those nuclei have no membrane-based cell; cellmeasurement
grows them itself (--estimate-cell-boundary-dist) and needs their pixels free.

The mask is processed in strips of rows with a margin of --distance-px rows on
each side, so memory use does not scale with the image height. The output is an
uncompressed TIFF so that cellmeasurement can memory-map it.
"""

import argparse
import sys

import numpy as np
import tifffile
from scipy import ndimage as ndi


def open_mask(path):
    '''
    Return a single-plane label TIFF as a 2-D array.

    Uses a read-only memory map when the file is uncompressed and contiguous.
    The segmenters write compressed masks, which are decoded into memory.
    '''
    try:
        arr = tifffile.memmap(path, mode='r')
    except ValueError:
        arr = tifffile.imread(path)
    arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"{path}: expected a single-plane label mask, got shape {arr.shape}")
    return arr


def expand_block(labels, distance_px, blocked=None):
    '''
    Expand the labels in a 2-D block by up to distance_px pixels.

    Background pixels take the label of the nearest labelled pixel if it is at
    most distance_px away. Pixels where blocked is True stay background.
    '''
    if distance_px <= 0 or not labels.any():
        return labels.copy()
    background = labels == 0
    distances, (rows, cols) = ndi.distance_transform_edt(background, return_indices=True)
    grow = background & (distances <= distance_px)
    if blocked is not None:
        grow &= ~blocked
    out = labels.copy()
    out[grow] = labels[rows[grow], cols[grow]]
    return out


def expand_mask(whole_cell, distance_px, out, nuclear=None, strip_rows=2048):
    '''
    Expand a full mask strip by strip and write the result into out.

    Each strip is read with distance_px extra rows above and below. A pixel's
    nearest label within distance_px is always inside that margin, so the
    result matches expanding the whole mask in one pass.
    '''
    height = whole_cell.shape[0]
    for r0 in range(0, height, strip_rows):
        r1 = min(height, r0 + strip_rows)
        h0 = max(0, r0 - distance_px)
        h1 = min(height, r1 + distance_px)
        block = np.asarray(whole_cell[h0:h1])
        blocked = None
        if nuclear is not None:
            blocked = (np.asarray(nuclear[h0:h1]) > 0) & (block == 0)
        expanded = expand_block(block, distance_px, blocked)
        out[r0:r1] = expanded[r0 - h0:r1 - h0]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--whole-cell-mask", required=True,
                        help="Whole-cell label mask TIFF to expand.")
    parser.add_argument("--nuclear-mask", default=None,
                        help="Nuclear label mask TIFF. Nucleus pixels outside every "
                             "whole-cell label are not given to any cell.")
    parser.add_argument("--distance-px", type=int, required=True,
                        help="Maximum expansion distance in pixels.")
    parser.add_argument("--output", required=True,
                        help="Output label mask TIFF (uncompressed, same dtype as input).")
    parser.add_argument("--strip-rows", type=int, default=2048,
                        help="Rows processed at a time (default: 2048).")
    args = parser.parse_args(argv)

    if args.distance_px < 1:
        parser.error("--distance-px must be at least 1")
    if args.strip_rows < 1:
        parser.error("--strip-rows must be at least 1")

    whole_cell = open_mask(args.whole_cell_mask)
    nuclear = None
    if args.nuclear_mask:
        nuclear = open_mask(args.nuclear_mask)
        if nuclear.shape != whole_cell.shape:
            sys.exit(f"Nuclear mask shape {nuclear.shape} does not match "
                     f"whole-cell mask shape {whole_cell.shape}")

    print(f"[expand_cells] {args.whole_cell_mask}: shape={whole_cell.shape} "
          f"dtype={whole_cell.dtype} distance_px={args.distance_px} "
          f"nuclear_mask={'yes' if nuclear is not None else 'no'}",
          file=sys.stderr)

    out = tifffile.memmap(args.output, shape=whole_cell.shape, dtype=whole_cell.dtype)
    expand_mask(whole_cell, args.distance_px, out, nuclear, args.strip_rows)
    out.flush()
    del out


if __name__ == "__main__":
    main()
