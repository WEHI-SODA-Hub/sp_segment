# Expand cells

Use `cell_expansion_microns` when the cell masks are smaller than the cells,
for example when the membrane markers mark only some cells and the nuclear
channel is used as the membrane channel. The masks are then nucleus-sized, and
membrane markers are measured over the nucleus only.

## How expansion works

Each background pixel within `cell_expansion_microns` of a cell is added to the
nearest cell. Growth stops where two cells meet, so cells do not overlap. The
distance is converted to pixels with `pixel_size_microns` and rounded, with a
minimum of 1 pixel.

Expansion runs after segmentation and mask smoothing, and before cell
measurement. The published GeoJSON and `*_mask.tiff` include the expansion.
The KRONOS2 input is the segmented whole-cell mask, without the expansion.

## Grow every cell

```bash
--membrane_channels DAPI     # in the samplesheet
--cell_expansion_microns 1.5
```

With `use_whole_cell_only true` (the default), every cell is grown by 1.5 µm.

## Grow only cells the membrane channels did not change

Use this when some cells have membrane signal and others do not: cells with a
membrane outline keep it, and nucleus-sized cells are grown.

```bash
--cell_expansion_microns 1.5
--expand_only_without_membrane true
```

The pipeline segments each sample twice: from the nuclear channel alone, and
from the nuclear and membrane channels (the whole-cell segmentation). Both
use the same segmenter and settings. Each whole-cell mask is compared with
the nuclear-channel mask it overlaps best, by intersection over union (IoU):

| Case                                                         | Result                        |
| ------------------------------------------------------------ | ----------------------------- |
| Whole-cell mask with IoU of at least `expand_match_iou`      | Grown by the distance.        |
| Whole-cell mask with a lower IoU, or no nuclear-channel mask | Kept as segmented, not grown. |
| Nuclear-channel mask that touches no whole-cell mask         | Added as a cell and grown.    |

No cell is removed. A pixel is given to a grown cell only if that cell is
the nearest cell, so grown cells do not take space next to cells that keep
their outline.

`expand_match_iou` (default `0.8`, hidden in `--help`) sets how close the two
masks must be. If one mask lies inside the other, an IoU of 0.8 means the
larger mask has 1.25 times the area of the smaller one. Raise it to grow fewer
cells. The two segmentations do not outline a nucleus in exactly the same
way, so values above about 0.85 also miss cells that have no membrane signal.

!!! warning "Cells that are almost all nucleus are also grown"

    The option cannot tell a cell with no membrane signal from a cell that has
    almost no cytoplasm, for example tumour cells in dense areas. Both have a
    whole-cell mask that matches the nucleus, so both are grown, although the
    outline of the second is already correct. Inside densely packed areas
    growth is small, because cells stop where they meet. Cells at the edge of
    a dense area can grow by the full distance into the space next to them.
    Check the outlines in such areas, and use a smaller
    `cell_expansion_microns` or leave the option off if the cells there should
    not grow.

Both segmentations run, so segmentation takes up to about twice as long as in
whole-cell-only mode. If `use_whole_cell_only` is `true`, the run logs a
warning about this at startup, and no nuclear measurements are produced. With
`use_whole_cell_only false`, cell measurement then pairs the cells with the
nuclear-channel masks as usual.

The KRONOS2 input is the whole-cell mask, so it does not include the cells
added from the nuclear-channel segmentation.

## Grow every cell in paired mode

With `use_whole_cell_only false` and `expand_only_without_membrane false`,
whole-cell masks are grown first, but not into nuclei that have no whole-cell
mask. Those nuclei are then grown by the same distance into the pixels that
are left.

## Choosing the distance

Measure how far a membrane marker extends past the nucleus on cells that are
positive for it, and pick a distance that covers the bright ring next to the
nucleus. Larger distances mostly add space between cells, which belongs to no
cell. As a reference, a 1.5 µm expansion of a nucleus 7 µm across doubles its
area.

`cell_expansion_microns` replaces `estimate_cell_boundary_dist`, which sets the
growth of unmatched nuclei in pixels, when both are set.
