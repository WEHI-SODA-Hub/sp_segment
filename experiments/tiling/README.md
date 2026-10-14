# Tiling and straight cell edges

Tests of where straight horizontal and vertical cell edges come from, and of
fixes for them. The method, results and conclusions are in
[docs/background/tiling-and-straight-edges.md](../../docs/background/tiling-and-straight-edges.md).

The scripts read WEHI paths (the CSL atherosclerosis project and
`/vast/scratch/users/mckay.m/sp_segment_tests`) and were run on 2026-10-14.
Scripts that call Cellpose or sopa run inside the pipeline's sopa-Cellpose
container (`community.wave.seqera.io/library/python_pip_sopacellpose_cellpose:2bb51160896b005b`)
on a GPU. The others need numpy, tifffile, scipy, scikit-image, shapely and
pyarrow.

## Patches (sopa Cellpose path)

- `segment_patches.py`: segments a 4 × 4 block of 3000 px patches with 250 px
  overlap, with per-patch or global normalisation, and writes each patch's
  polygons with their mean cell probability.
- `whole_region.py`: segments the same region in one Cellpose call, as a
  reference without patches.
- `resolve_and_score.py`: stitches the patches with sopa's `solve_conflicts`,
  centre ownership, and one version per cell (by context or by cell
  probability), and scores each result.
- `tiling_scores.csv`: results of `resolve_and_score.py`.
- `production_patch_edges.py`: straight outline edges on sopa patch edges in
  production images.
- `production_strip_agreement.py`: agreement of two patches' outlines in the
  shared strips of one production image.

## CellSAM

- `cellsam_seams.py`: cells cut in two or left as half cells on CellSAM's
  600 px block seams.
- `straight_splits.py`: straight boundaries between two cells in Cellpose and
  CellSAM masks.

## Cellpose's own straight edges

- `cellpose_variants.py`: Cellpose on one crop with default settings, `niter`
  1000 and 2000, flip augmentation, and the image rotated 45°.
- `cellpose_variant_scores.csv`: straight edges for those runs and for
  `tile_overlap` 0.1 and 0.5.
- `ws_masks.py`: Cellpose mask construction with a watershed instead of the
  11 × 11 seed box. `ws_scores.csv`: results.
- `rot_consensus.py`: averages the flows of the 0° and 45° runs before building
  masks. `flow_average_scores.csv`: results.
- `rotation_swap.py`: replaces cells that have a straight edge with the matching
  cell from the 45° run. `swap_scores.csv`: results.
