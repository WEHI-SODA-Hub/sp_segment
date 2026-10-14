# Tiling and straight cell edges

Some cell outlines have straight horizontal or vertical edges that do not
follow the image. This page records where they come from, which ones are
caused by tiling, and what was tested to remove them. The tests used the CSL
atherosclerosis COMET images (0.28 µm/pixel) on 2026-10-14.

## Summary

| Source                                    | Segmenter | Caused by tiling | Size of the problem                               |
| ----------------------------------------- | --------- | ---------------- | ------------------------------------------------- |
| CellSAM block stitching                   | CellSAM   | Yes              | about 1% of cells cut at 600 px block seams       |
| sopa patch stitching                      | Cellpose  | Yes              | at or near chance level in the test region        |
| The Cellpose-SAM network output           | Cellpose  | No               | 1.6% of cells in the outlines; 8–23% in the masks |
| Image border (cells cut by the crop edge) | all       | No               | only at the edge of the image                     |

## How straight edges were counted

Two counts are used. Natural cell boundaries are rarely exactly aligned with
the pixel grid over 15 pixels (4.2 µm), so both count edges of that length
that are exactly horizontal or exactly vertical. Edges on the image border are
not counted.

- Outlines: a polygon segment of at least 15 pixels in the cell outlines, as
  written by sopa. sopa smooths outlines, so this count is lower.
- Masks: a run of at least 15 pixels along which the cell label changes between
  two neighbouring pixel rows or columns, in a label mask.

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

On the raw-maximum membrane input, CellSAM had cells cut in two or left as
half cells on the seams in every crop tested:

| Crop      | Cells | Labels cut in two | Half cells |
| --------- | ----- | ----------------- | ---------- |
| b6_immune | 1088  | 6                 | 6          |
| b6_macro  | 507   | 5                 | 1          |
| b2_tcell  | 957   | 0                 | 5          |

Across the 33 CellSAM masks of the membrane-combination tests, 204 straight
boundaries of at least 20 px between two cells were found, 37 of them on the
600 px grid. The 33 Cellpose masks of the same samples had 40, none on that
grid.

## Cellpose patches (sopa)

The Cellpose path segments `patch_width_pixel` (3000) pixel patches with
`patch_overlap_pixel` (250) pixels of overlap. Each patch is normalised to its
own 1st and 99th percentiles. sopa then merges two overlapping cells from
different patches into their union when the overlap is at least half the area
of the smaller cell (`solve_conflicts`, threshold 0.5), and otherwise keeps
both.

In 8 production images (129,000 cells), 79 of 1651 straight edges lay on a
patch edge, against about 3 expected by chance.

### Per-patch normalisation

In one production image, the 99th percentile used to normalise each patch
varied 5.3-fold for DAPI (112 to 589) and 7.1-fold for the combined membrane
channel (350 to 2493) across 38 patches with tissue. Two neighbouring patches
therefore see the same cells at different contrast. For cells in the shared
strips of that image, the two patches' outlines agreed (IoU 0.8 or higher) for
75.7% of cells, differed in shape (IoU 0.5 to 0.8) for 10.6%, and disagreed
(split, merged or missing; IoU below 0.5) for 13.7%.

### Test region

A 4 × 4 block of production-sized patches (16 patches, x 11000 to 22250,
y 13750 to 24360) of the batch 6 image 25_1592_1A AS_1 was segmented with
Cellpose (`cpsam_v2`, diameter 30, flow threshold 0.4, cell probability
threshold 0, minimum area 100 px). The region holds 756 manually labelled
cells: 137 within 150 px of a seam and 619 elsewhere. Two normalisations and
four stitching rules were compared, and one run of the whole region in a single
Cellpose call (no patches) was used as a reference.

Normalisation:

- Per patch: the current behaviour.
- Global: the 1st and 99th percentiles of each channel over the whole region,
  passed to Cellpose as fixed values.

Stitching:

- sopa: the current behaviour.
- One version per cell, by context: versions of a cell from different patches
  are grouped when they overlap by at least half the smaller one. One version
  is kept: one that does not touch a patch edge, then the one furthest from
  its patch edge. Remaining overlaps go to the cell with the higher mean cell
  probability.
- One version per cell, by cell probability: the same grouping, keeping the
  version with the higher mean cell probability.
- Centre ownership: a patch keeps the cells whose centre lies in its core (the
  patch minus half the overlap on each inner side).

| Normalisation | Stitching                   | Cells  | Overlapping pairs | Straight edges | On a patch edge | Seam labels in a cell / alone | Other labels in a cell / alone |
| ------------- | --------------------------- | ------ | ----------------- | -------------- | --------------- | ----------------------------- | ------------------------------ |
| Per patch     | sopa                        | 12,390 | 943               | 208            | 1               | 100.0% / 79.6%                | 99.5% / 94.7%                  |
| Per patch     | One version, by context     | 12,390 | 0                 | 214            | 0               | 100.0% / 79.6%                | 99.5% / 94.7%                  |
| Per patch     | One version, by probability | 12,390 | 0                 | 220            | 4               | 99.3% / 78.8%                 | 99.5% / 94.7%                  |
| Per patch     | Centre ownership            | 12,209 | 0                 | 214            | 0               | 97.1% / 76.6%                 | 99.5% / 94.7%                  |
| Global        | sopa                        | 12,330 | 952               | 196            | 2               | 100.0% / 79.6%                | 97.6% / 92.7%                  |
| Global        | One version, by context     | 12,330 | 0                 | 204            | 1               | 100.0% / 79.6%                | 97.6% / 92.7%                  |
| Global        | One version, by probability | 12,330 | 0                 | 217            | 7               | 99.3% / 78.8%                 | 97.6% / 92.7%                  |
| Global        | Centre ownership            | 12,186 | 0                 | 201            | 0               | 92.7% / 73.7%                 | 97.6% / 92.7%                  |
| Whole region  | none (no patches)           | 12,188 | 746               | 236            | 2               | 90.5% / 70.1%                 | 95.6% / 90.8%                  |

"In a cell" is the share of labelled cells whose centroid is inside a cell.
"Alone" also requires that no other labelled cell is in the same cell.
Overlapping pairs count pairs of outlines that share more than 1 px²; with
sopa, 25 of the 943 pairs overlap by more than 10% of the smaller cell, and
the rest are slivers from sopa's polygon smoothing.

Agreement of the two patches' outlines in the shared strips:

| Normalisation | Cells in strips | Same (IoU ≥ 0.8) | Shape differs | Disagree (IoU < 0.5) |
| ------------- | --------------- | ---------------- | ------------- | -------------------- |
| Per patch     | 1472            | 84.0%            | 7.4%          | 8.6%                 |
| Global        | 1425            | 89.5%            | 5.0%          | 5.5%                 |

Findings:

- Global normalisation made the two patches agree more often (84.0% to 89.5%)
  but found fewer labelled cells away from seams (99.5% to 97.6% in a cell).
  Buglakova et al. (2025) describe the same trade-off.
- Segmenting the whole region in one call scored lowest. Labelled cells near
  seams scored lower than the rest in every run, including the run without
  patches, so the difference comes from the tissue near those lines, not from
  the seams.
- In this region, straight edges on patch edges were at the chance level (0 to 7) for every method.
- Keeping one version per cell, chosen by context, scored the same as sopa and
  left no overlapping outlines. Choosing by cell probability or by centre
  ownership scored lower; centre ownership loses cells whose two versions have
  their centres on opposite sides of the midline.

## Straight edges from Cellpose itself

Cellpose produces straight edges even when it runs on a whole image without
patches. Running Cellpose directly on three 2048 × 2048 crops, outside sopa,
gave straight lines at the same places as the pipeline, for example the
40 px line between two stacked nuclei at row 517 of b6_immune. In those crops,
8.4% (b6_immune), 14.5% (b6_macro) and 22.6% (b2_vsmc) of cells have a straight
edge in the mask.

Tests on the crops b6_immune, b6_macro and b2_vsmc
(`cellpose_variant_scores.csv`):

| Change                             | Straight edges between a cell and background, b6_immune |
| ---------------------------------- | ------------------------------------------------------- |
| Default                            | 90                                                      |
| `niter` 1000 or 2000               | 89                                                      |
| Flip augmentation (`augment`)      | 97                                                      |
| `tile_overlap` 0.5 (default 0.1)   | 99                                                      |
| Image rotated 45° and rotated back | 12                                                      |

- The straight edges do not lie on Cellpose's internal 256 px tiles or on the
  8 or 16 px patches of the network.
- When Cellpose runs on the rotated image, its output has 55 straight edges
  along its own pixel axes in b6_immune. The default run keeps 56 after the
  same rotation and back, so both runs make about the same number, in
  different places. The straight edges follow the pixel grid that Cellpose runs
  on, not the tissue.
- The straight edges are already in the cell probability map: the boundary of
  cell probability > 0 has 112 straight edges in b6_immune, against 90 in the
  final mask. They therefore come from the network output, not from mask
  construction. Building masks with a watershed instead of Cellpose's 11 × 11
  seed box gave identical masks (`ws_scores.csv`).
- Averaging the flows and cell probability of the 0° and 45° runs before
  building masks did not remove them (`flow_average_scores.csv`), because the
  straight edges of the 0° map remain in the average.

Taking whole cells from a second run at 45° does remove them. Where a cell
from the default run has a straight edge, and the 45° run has a matching cell
(IoU ≥ 0.5) without one, the 45° cell replaces it (`swap_scores.csv`):

| Crop      | Cells with a straight edge, default | After replacement | Labelled cells alone, default | After replacement |
| --------- | ----------------------------------- | ----------------- | ----------------------------- | ----------------- |
| b6_immune | 108                                 | 28                | 87.0%                         | 87.0%             |
| b6_macro  | 77                                  | 20                | 99.3%                         | 99.3%             |
| b2_vsmc   | 56                                  | 9                 | 100.0%                        | 100.0%            |

Cell counts and median areas did not change. The rotated image is padded so
that no tissue is cut (2924 × 2924 px for a 2048 px crop), so the second run
takes about twice the Cellpose GPU time of the first.
This is not part of the pipeline.

## Not tested or not done

- The straight-edge measure does not separate wrong straight edges from
  correct ones (for example, the flat contact between two pressed nuclei).
- No hand-drawn outlines were available, so outline accuracy was not measured.
- The fixes on this page are experiments in `experiments/tiling/`. The pipeline
  still uses sopa stitching, per-patch normalisation and a single Cellpose run.

## Related work

- Buglakova et al. (2025), _Tiling artifacts and trade-offs of feature
  normalization in the segmentation of large biological images_,
  [arXiv:2503.19545](https://arxiv.org/abs/2503.19545).
- Squidpy `calculate_tiling_qc` and `assign_stitch_groups` (experimental) flag
  cells with long straight, axis-aligned edges and group pieces cut at tile
  seams.
- StarDist `predict_instances_big` discards a context margin from each block
  before stitching.

## Repeating the tests

The scripts and result tables are in `experiments/tiling/`; see its
`README.md`. They read WEHI paths and run in the pipeline's sopa-Cellpose
container.

!!! note

    Cellpose's `model.eval` normalises its input array in place. Pass a copy
    when calling it more than once on the same array or on overlapping views of
    one array.
