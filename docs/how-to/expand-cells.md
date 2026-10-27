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

## Grow only cells with no membrane signal

This follows 10x Xenium multimodal cell segmentation: a cell found from the
membrane channels keeps its outline, and a nucleus with no such cell is grown
to make one.

```bash
--cell_expansion_microns 1.5
--expand_only_without_membrane true
```

The pipeline segments nuclei from the nuclear channel and whole cells from the
nuclear and membrane channels. Nuclei are segmented for this whatever
`use_whole_cell_only` is set to. Cell measurement then pairs each nucleus with
the whole-cell mask it overlaps most, or else with the nearest unpaired
whole-cell mask within `dist_threshold` pixels:

| Case                                    | Result                                           |
| --------------------------------------- | ------------------------------------------------ |
| Nucleus with a whole-cell mask          | The cell is the whole-cell mask, not grown.      |
| Nucleus with no whole-cell mask         | The cell is the nucleus grown by the distance.   |
| Whole-cell mask with no nucleus         | Dropped.                                         |

Both segmentations run, so segmentation takes up to about twice as long as in
whole-cell-only mode. If `use_whole_cell_only` is `true`, the run logs a
warning about this at startup. Nuclear measurements are also produced.

Membrane cells cannot be told apart from nucleus-sized cells using the
whole-cell mask and the images alone, so the separate nuclear segmentation is
needed.

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
