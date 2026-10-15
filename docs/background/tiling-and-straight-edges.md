# Tiling and straight cell edges

Some cell outlines have straight horizontal or vertical edges that do not
follow the image. This page describes where they come from, which ones are
caused by tiling, and what was found when fixes were tested.

## Summary

| Source                          | Segmenter | Caused by tiling | Effect                                             |
| ------------------------------- | --------- | ---------------- | -------------------------------------------------- |
| CellSAM block stitching         | CellSAM   | Yes              | some cells are cut or halved at block seams        |
| sopa patch stitching            | Cellpose  | Yes              | straight edges on patch edges were rare in testing |
| The Cellpose-SAM network output | Cellpose  | No               | straight edges appear with or without patches      |
| Image border                    | all       | No               | cells cut by the edge of the image                 |

A straight edge here is a part of a cell boundary at least 15 pixels long that
is exactly horizontal or exactly vertical. Natural cell boundaries are rarely
aligned with the pixel grid over that length. Edges on the image border are
not counted.

## CellSAM block seams

CellSAM (`cellsam_use_wsi = true`) segments the image in blocks of
`cellsam_block_size` (600) pixels with `cellsam_overlap` (250) pixels of
padding on each side. It then joins cells across a seam when the two blocks'
versions overlap with an IoU of at least `cellsam_iou_threshold` (0.5), and
trims every block back to its 600 × 600 core. Each cell that crosses a seam is
therefore made from two predictions, one per side. This gives three faults:

- The two predictions are joined but differ in shape, so the outline steps at
  the seam.
- The two predictions overlap with an IoU below 0.5, so the cell stays as two
  cells split along the seam.
- One block did not detect the cell, so only half of it is kept.

In testing, straight boundaries between CellSAM cells often lay on the 600 px
block grid. Cellpose masks of the same images had none on that grid.

## Cellpose patches (sopa)

The Cellpose path segments `patch_width_pixel` (3000) pixel patches with
`patch_overlap_pixel` (250) pixels of overlap. Each patch is normalised to its
own 1st and 99th percentiles. sopa then merges two overlapping cells from
different patches into their union when the overlap is at least half the area
of the smaller cell (`solve_conflicts`, threshold 0.5), and otherwise keeps
both.

Because each patch is normalised on its own, two neighbouring patches can see
the same cells at different contrast, and their outlines of a cell in the
shared strip can differ.

Normalisations and stitching rules compared:

- Normalisation per patch (the current behaviour), or with the 1st and 99th
  percentiles of each channel over the whole region passed to Cellpose as fixed
  values (global).
- Stitching by sopa (the current behaviour); keeping one version per cell,
  chosen by context (a version that does not touch a patch edge, then the one
  furthest from its patch edge); keeping one version per cell, chosen by mean
  cell probability; or centre ownership (a patch keeps the cells whose centre
  lies in its core).
- One run of the whole region in a single Cellpose call, without patches.

Findings:

- Global normalisation made the two patches' outlines agree more often but
  found fewer labelled cells away from seams. Buglakova et al. (2025) describe
  the same trade-off.
- Segmenting the whole region in one call scored lowest.
- Straight edges on patch edges were at about the level expected by chance for
  every method.
- Keeping one version per cell, chosen by context, scored the same as sopa and
  left no overlapping outlines. Choosing by cell probability or by centre
  ownership scored lower; centre ownership loses cells whose two versions have
  their centres on opposite sides of the midline.

## Straight edges from Cellpose itself

Cellpose produces straight edges even when it runs on a whole image without
patches, at the same places as in the pipeline.

Changes tested, with Cellpose run directly on the image:

| Change                             | Effect on straight edges |
| ---------------------------------- | ------------------------ |
| `niter` 1000 or 2000               | none                     |
| Flip augmentation (`augment`)      | none                     |
| `tile_overlap` 0.5 (default 0.1)   | none                     |
| Watershed instead of the seed box  | none (identical masks)   |
| Image rotated 45° and rotated back | most removed             |

- The straight edges do not lie on Cellpose's internal 256 px tiles or on the
  8 or 16 px patches of the network.
- On a rotated image, Cellpose makes about as many straight edges along its own
  pixel axes, in different places. The straight edges follow the pixel grid
  that Cellpose runs on, not the tissue.
- The straight edges are already in the cell probability map, so they come
  from the network output, not from mask construction.
- Averaging the flows and cell probability of the 0° and 45° runs before
  building masks did not remove them, because the straight edges of the 0° map
  remain in the average.

Taking whole cells from a second run at 45° does remove most of them. Where a
cell from the default run has a straight edge, and the 45° run has a matching
cell (IoU ≥ 0.5) without one, the 45° cell replaces it. Cell counts, median
areas and labelled cells found did not change. The rotated image is padded so
that no tissue is cut, so the second run takes about twice the Cellpose GPU
time of the first. This is not part of the pipeline.

## Not tested or not done

- The straight-edge measure does not separate wrong straight edges from
  correct ones (for example, the flat contact between two pressed nuclei).
- No hand-drawn outlines were available, so outline accuracy was not measured.
- None of the fixes on this page is in the pipeline. It uses sopa stitching,
  per-patch normalisation and a single Cellpose run.

## Related work

- Buglakova et al. (2025), _Tiling artifacts and trade-offs of feature
  normalization in the segmentation of large biological images_,
  [arXiv:2503.19545](https://arxiv.org/abs/2503.19545).
- Squidpy `calculate_tiling_qc` and `assign_stitch_groups` (experimental) flag
  cells with long straight, axis-aligned edges and group pieces cut at tile
  seams.
- StarDist `predict_instances_big` discards a context margin from each block
  before stitching.

!!! note

    Cellpose's `model.eval` normalises its input array in place. Pass a copy
    when calling it more than once on the same array or on overlapping views of
    one array.
