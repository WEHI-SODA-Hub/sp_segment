#!/usr/bin/env python3
"""
Expand the cells in a whole-cell label mask by a fixed distance.

Each background pixel within --distance-px of a cell is given the label of the
nearest cell, so growth stops where two cells meet. This is the same rule as
skimage.segmentation.expand_labels.

Two modes:

- Expand every cell (default). With --nuclear-mask, pixels that belong to a
  nucleus but to no whole-cell label are left as background. Those nuclei have
  no membrane-based cell; cellmeasurement grows them itself
  (--estimate-cell-boundary-dist) and needs their pixels free.

- --match-iou T (needs --nuclear-mask): expand only the cells that the membrane
  channels did not change. The nuclear mask is a segmentation of the nuclear
  channel alone. A whole-cell label whose best-matching nucleus has an
  intersection over union of at least T is expanded. Each nucleus that touches
  no whole-cell label is added as a new cell and expanded. Every other
  whole-cell label keeps its outline. No label is removed. A pixel is given to
  an expanded cell only if that cell is the nearest cell of all, so expanded
  cells do not grow past cells that keep their outline.

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


def expand_block(labels, distance_px, blocked=None, grow_labels=None):
    '''
    Expand the labels in a 2-D block by up to distance_px pixels.

    Background pixels take the label of the nearest labelled pixel if it is at
    most distance_px away. Pixels where blocked is True stay background. With
    grow_labels (a boolean lookup indexed by label), a pixel is filled only if
    its nearest label is one marked True.
    '''
    if distance_px <= 0 or not labels.any():
        return labels.copy()
    background = labels == 0
    distances, (rows, cols) = ndi.distance_transform_edt(background, return_indices=True)
    grow = background & (distances <= distance_px)
    if blocked is not None:
        grow &= ~blocked
    nearest = labels[rows[grow], cols[grow]]
    if grow_labels is not None:
        keep = grow_labels[nearest]
        grow[grow] = keep
        nearest = nearest[keep]
    out = labels.copy()
    out[grow] = nearest
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


def _add_counts(total, counts):
    '''Add a bincount result to a running total, growing it as needed.'''
    if len(counts) > len(total):
        total = np.pad(total, (0, len(counts) - len(total)))
    total[:len(counts)] += counts
    return total


def match_to_nuclei(whole_cell, nuclear, match_iou, strip_rows=2048):
    '''
    Compare each whole-cell label with the nuclear labels it overlaps.

    Returns (expand_cell, new_nuclei, best_iou):
      expand_cell  boolean array indexed by whole-cell label: best IoU >= match_iou
      new_nuclei   boolean array indexed by nuclear label: the nucleus touches
                   no whole-cell label
      best_iou     float array indexed by whole-cell label (0 = no overlap)
    '''
    height = whole_cell.shape[0]
    cell_area = np.zeros(1, dtype=np.int64)
    nuc_area = np.zeros(1, dtype=np.int64)
    pair_keys, pair_counts = [], []
    for r0 in range(0, height, strip_rows):
        cell = np.asarray(whole_cell[r0:r0 + strip_rows]).ravel().astype(np.int64)
        nuc = np.asarray(nuclear[r0:r0 + strip_rows]).ravel().astype(np.int64)
        cell_area = _add_counts(cell_area, np.bincount(cell))
        nuc_area = _add_counts(nuc_area, np.bincount(nuc))
        both = (cell > 0) & (nuc > 0)
        pairs = (cell[both].astype(np.uint64) << np.uint64(32)) | nuc[both].astype(np.uint64)
        keys, counts = np.unique(pairs, return_counts=True)
        pair_keys.append(keys)
        pair_counts.append(counts)
    keys = np.concatenate(pair_keys) if pair_keys else np.zeros(0, np.uint64)
    counts = np.concatenate(pair_counts) if pair_counts else np.zeros(0, np.int64)
    keys, inverse = np.unique(keys, return_inverse=True)
    inter = np.bincount(inverse, weights=counts).astype(np.int64) if len(keys) else counts
    cell_ids = (keys >> np.uint64(32)).astype(np.int64)
    nuc_ids = (keys & np.uint64(0xFFFFFFFF)).astype(np.int64)

    iou = inter / (cell_area[cell_ids] + nuc_area[nuc_ids] - inter)
    best_iou = np.zeros(len(cell_area))
    np.maximum.at(best_iou, cell_ids, iou)
    expand_cell = best_iou >= match_iou
    expand_cell[0] = False

    new_nuclei = nuc_area > 0
    new_nuclei[0] = False
    new_nuclei[nuc_ids] = False
    return expand_cell, new_nuclei, best_iou


def expand_matched(whole_cell, nuclear, distance_px, out, match_iou, strip_rows=2048):
    '''
    Expand only the cells that match their nucleus, and add + expand the
    nuclei that have no whole-cell label. Writes into out; returns a summary.

    New cells from nuclei are labelled max(whole-cell label) + nuclear label.
    '''
    expand_cell, new_nuclei, _ = match_to_nuclei(whole_cell, nuclear, match_iou, strip_rows)
    offset = len(expand_cell) - 1
    max_label = offset + len(new_nuclei) - 1
    if max_label > np.iinfo(out.dtype).max:
        raise ValueError(f"output dtype {out.dtype} cannot hold label {max_label}")
    grow_labels = np.concatenate([expand_cell, new_nuclei[1:]])

    height = whole_cell.shape[0]
    for r0 in range(0, height, strip_rows):
        r1 = min(height, r0 + strip_rows)
        h0 = max(0, r0 - distance_px)
        h1 = min(height, r1 + distance_px)
        block = np.asarray(whole_cell[h0:h1]).astype(out.dtype)
        nuc = np.asarray(nuclear[h0:h1])
        add = (block == 0) & new_nuclei[nuc]
        block[add] = nuc[add].astype(out.dtype) + offset
        expanded = expand_block(block, distance_px, grow_labels=grow_labels)
        out[r0:r1] = expanded[r0 - h0:r1 - h0]
    return {
        "expanded_whole_cell": int(expand_cell.sum()),
        "added_nuclei": int(new_nuclei.sum()),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--whole-cell-mask", required=True,
                        help="Whole-cell label mask TIFF to expand.")
    parser.add_argument("--nuclear-mask", default=None,
                        help="Nuclear label mask TIFF. Without --match-iou, nucleus pixels "
                             "outside every whole-cell label are not given to any cell.")
    parser.add_argument("--distance-px", type=int, required=True,
                        help="Maximum expansion distance in pixels.")
    parser.add_argument("--match-iou", type=float, default=None,
                        help="Expand only whole-cell labels whose best-matching nucleus has "
                             "at least this intersection over union, plus nuclei that touch "
                             "no whole-cell label (added as new cells). Needs --nuclear-mask.")
    parser.add_argument("--output", required=True,
                        help="Output label mask TIFF (uncompressed).")
    parser.add_argument("--strip-rows", type=int, default=2048,
                        help="Rows processed at a time (default: 2048).")
    args = parser.parse_args(argv)

    if args.distance_px < 1:
        parser.error("--distance-px must be at least 1")
    if args.strip_rows < 1:
        parser.error("--strip-rows must be at least 1")
    if args.match_iou is not None:
        if not args.nuclear_mask:
            parser.error("--match-iou needs --nuclear-mask")
        if not 0 < args.match_iou <= 1:
            parser.error("--match-iou must be above 0 and at most 1")

    whole_cell = open_mask(args.whole_cell_mask)
    nuclear = None
    if args.nuclear_mask:
        nuclear = open_mask(args.nuclear_mask)
        if nuclear.shape != whole_cell.shape:
            sys.exit(f"Nuclear mask shape {nuclear.shape} does not match "
                     f"whole-cell mask shape {whole_cell.shape}")

    print(f"[expand_cells] {args.whole_cell_mask}: shape={whole_cell.shape} "
          f"dtype={whole_cell.dtype} distance_px={args.distance_px} "
          f"nuclear_mask={'yes' if nuclear is not None else 'no'} "
          f"match_iou={args.match_iou}",
          file=sys.stderr)

    if args.match_iou is None:
        out = tifffile.memmap(args.output, shape=whole_cell.shape, dtype=whole_cell.dtype)
        expand_mask(whole_cell, args.distance_px, out, nuclear, args.strip_rows)
    else:
        # Added nuclei take labels above the whole-cell labels: use 32 bits.
        out = tifffile.memmap(args.output, shape=whole_cell.shape, dtype=np.uint32)
        summary = expand_matched(whole_cell, nuclear, args.distance_px, out,
                                 args.match_iou, args.strip_rows)
        print(f"[expand_cells] expanded {summary['expanded_whole_cell']} whole-cell labels "
              f"that match their nucleus (IoU >= {args.match_iou}); added and expanded "
              f"{summary['added_nuclei']} nuclei with no whole-cell label",
              file=sys.stderr)
    out.flush()
    del out


if __name__ == "__main__":
    main()
