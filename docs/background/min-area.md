# Minimum cell area

Cellpose and CellSAM remove cells smaller than a minimum area. Both filters
default to 200 pixels², so the same image is filtered the same way whichever
method segments it.

| Parameter           | Default     | Measures                                  | Applied                                    |
| ------------------- | ----------- | ----------------------------------------- | ------------------------------------------ |
| `cellpose_min_area` | 200 pixels² | area of the cell polygon, after smoothing | by sopa, on each patch during segmentation |
| `cellsam_min_area`  | 200 pixels² | number of pixels in the cell's label      | once, after the tiles are merged           |

- Both keep cells whose area is greater than or equal to the threshold.
- `0` keeps every cell.
- The unit is pixels², not µm², so the value does not follow
  `pixel_size_microns`. At 0.28 µm/pixel (COMET), 200 pixels² is 15.7 µm², about
  a 4.5 µm diameter disc. For another pixel size, recalculate the value.

## Polygon area against pixel count

The two filters measure area differently. On round and ragged test cells, the
polygon area was 0.93 to 1.00 times the pixel count (median 0.99). So 200 on the
CellSAM side is about the same as 202 on the Cellpose side.

Cellpose applies the filter per patch. This does not remove cells cut by a patch
edge: patches overlap by `patch_overlap_pixel`, so a cut cell is whole in the
neighbouring patch, and the whole copy is the one kept.

## Mesmer

`mesmer_min_nuclei_area` (default 40 pixels²) is a separate filter. It removes
small nuclei, not cells, and is applied inside the Mesmer tool.
