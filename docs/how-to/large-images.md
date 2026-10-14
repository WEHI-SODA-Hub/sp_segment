# Large and dense images

Whole-slide images can need hundreds of GB of memory and many hours per task.
Memory use depends more on the number of cells than on the file size.

## Choose a resource profile

| Image                                                | Profile                    |
| ---------------------------------------------------- | -------------------------- |
| up to 150 GB                                         | `small`                    |
| up to 300 GB                                         | `medium`                   |
| up to 600 GB                                         | `large`                    |
| millions of cells (for example large tumour regions) | `wehi_med` or `wehi_large` |

The full requests are in [Resource profiles](../reference/profiles.md). If a
task fails for lack of memory, it is retried once with twice the memory. If it
fails again, move to the next profile, or raise that one process with a `-c`
config file. See [Changing resources](../usage.md#changing-resources).

## Reduce memory

- List only the membrane channels you need in `membrane_channels`.
- Keep `gzip_geojson = true` (the default). Uncompressed GeoJSON for a whole
  slide can exceed 100 GB.
- Leave `neighbors = 0` (the default) unless you need neighbourhood features.
- For a quick look, set `downsample_factor 2` (or higher) to measure on a
  downsampled image. Area and length measurements are scaled back to full
  resolution.
- Set `enable_measurements false` to write cell outlines only.

## Speed up cell measurement

Cell measurement runs on tiles in parallel, one tile per CPU.

- `tile_size` (default 2048): aim for 2 to 4 times more tiles than CPUs.
- `tile_overlap` (default 200): at least the diameter of the largest cell.

Polygon extraction is the slowest part of cell measurement on images with
millions of cells. If the task can hit its time limit, set
`geometry_checkpoint_dir`. See
[Resume a killed measurement task](resume-measurement.md).

## Cellpose patches

Cellpose runs on square patches of `patch_width_pixel` pixels (default 3000)
that overlap by `patch_overlap_pixel` (default 250). Each patch is a separate
GPU task.

- Patch size changes the number of tasks, not the result: with the default
  preprocessing, about 82-87% of cells keep the same outline (IoU ≥ 0.90) when
  the patch width changes from 1500 to 3000 pixels. See
  [Cellpose settings](../background/cellpose.md#patch-size).
- Larger patches give fewer, longer tasks. This helps when the GPU queue is
  busy, because each task waits in the queue.
- Set `patch_overlap_pixel` to at least the diameter of the largest cell.

## QuPath cannot load the GeoJSON

Cells with very complex outlines can make QuPath fail with
`StackOverflowError`. Set `--smooth_masks true` to simplify the outlines before
measurement. `--smooth_method shapely` reduces the number of vertices directly.
See [Troubleshooting](../troubleshooting.md).
