# Cellpose settings

This page explains the Cellpose defaults and the measurements behind them. For
the parameter list, see
[Parameters](../reference/parameters.md#cellpose-options).

All models the pipeline offers (`cpsam_v2`, `cpsam`, `cpdino`, `cpdino-vitb`)
are Cellpose 4 transformer models. The defaults are the values the Cellpose
authors use at test time. For most data, choose a model and change nothing
else.

## What each setting does

| Parameter                      | Default    | Change it?                                                                                                       |
| ------------------------------ | ---------- | ---------------------------------------------------------------------------------------------------------------- |
| `cellpose_pretrained_model`    | `cpsam_v2` | Yes. This is the main choice. See [Choose a segmentation method](../how-to/choose-segmenter.md#cellpose-models). |
| `cellpose_diameter`            | 30         | Only if cells are outside 7.5-120 pixels.                                                                        |
| `cellpose_flow_threshold`      | 0.4        | Only if masks are visibly missing or ragged.                                                                     |
| `cellpose_cellprob_threshold`  | 0.0        | Only if faint cells are missed or background is segmented.                                                       |
| `cellpose_min_area`            | 100        | Set it for your cell type and pixel size. See [Minimum cell area](min-area.md).                                  |
| `cellpose_gaussian_sigma`      | 0          | No. See [Preprocessing](#preprocessing).                                                                         |
| `cellpose_clip_limit`          | 0          | No. See [Preprocessing](#preprocessing).                                                                         |
| `cellpose_clahe_kernel_size`   | not set    | No effect while `cellpose_clip_limit` is 0.                                                                      |
| `cellpose_tile_norm_blocksize` | 0          | No. See [Tile-local normalisation](#tile-local-normalisation).                                                   |
| `patch_width_pixel`            | 3000       | Only for memory or throughput. See [Patch size](#patch-size).                                                    |

### Diameter

Cellpose resizes the image by `30 / cellpose_diameter` before segmentation and
resizes the masks back afterwards. The default of 30 means no resizing.

Cellpose-SAM was trained on images resized between 0.25 and 4 times around a
mean cell diameter of 30 pixels, so it handles cells from about 7.5 to 120
pixels without a diameter setting. Set the real diameter only for cells outside
that range. A larger value, for example 90, shrinks the image 3 times and makes
the run faster, but small objects are resolved less well.

Cellpose 4 has no automatic diameter estimate. A diameter of 0 means "do not
resize", not "estimate". (The QuPath Cellpose extension's `.diameter(0.0)`
estimate comes from Cellpose 2 and 3.)

The diameter does not remove any cells. Use `cellpose_min_area` for that.

### Flow threshold

Cellpose recomputes the flows that each finished mask implies and compares them
with the flows the network predicted. Masks with a mean squared error above
`cellpose_flow_threshold` are removed. This filter acts on shape only.

- Raise it (0.6, 0.8) if visible cells are missing.
- Lower it (0.2, 0.3) if you get ragged masks.
- `0` turns the check off.

### Cell probability threshold

Pixels with a cell probability above `cellpose_cellprob_threshold` are used to
build cells. This changes which cells are found and where each boundary falls,
so it changes cell areas and every per-cell intensity measured over them.

The value is a logit, not a probability between 0 and 1. The useful range is
about -6 to +6, and 0 is the midpoint.

- Lower it (-1, -2) to find faint cells and make cells larger.
- Raise it (1, 2) to stop segmenting dim background and make cells smaller.

### Model type

`cellpose_model_type` no longer exists. Cellpose 4 ignores `--model-type`, so
older runs that set `cyto3` actually used `cpsam`. Setting it now stops the
pipeline with an error. Use `cellpose_pretrained_model`.

## Preprocessing

Cellpose normalises each image by rescaling its 1st percentile to 0 and its
99th percentile to 1. It does nothing else: its `sharpen_radius` and
`smooth_radius` default to 0, and it has no CLAHE step. The Cellpose-SAM paper
states that its training images were normalised this way.

sopa, which runs Cellpose on each patch, adds two steps of its own before
Cellpose:

| Step                                       | Added by | Setting in this pipeline           |
| ------------------------------------------ | -------- | ---------------------------------- |
| gaussian filter (sopa default `sigma = 1`) | sopa     | off: `cellpose_gaussian_sigma = 0` |
| CLAHE (sopa default `clip_limit = 0.2`)    | sopa     | off: `cellpose_clip_limit = 0`     |
| rescale to 1st and 99th percentiles        | Cellpose | on                                 |

To use sopa's steps, set `--cellpose_clip_limit 0.2 --cellpose_gaussian_sigma 1`.

### Why CLAHE is off

Every normalisation step is computed from the patch, so the input a cell gets
depends on what else is in its patch. Patch overlap does not change this:
overlap handles cells that cross a patch edge, not differences in
normalisation.

CLAHE caused most of this dependence. Its kernel defaults to 1/8 of the patch
size, so the contrast-enhancement scale changes with `patch_width_pixel`: at
0.28 µm/pixel, 52 µm for 1500-pixel patches, 105 µm for 3000 and 210 µm for 6000.

Cellpose-SAM was trained with random changes to global brightness and contrast,
blur and noise, so it tolerates those. CLAHE is a local, position-dependent
remapping, which was not part of its training.

Measured on the DAPI channel of a whole-slide COMET image: a
1024 × 1024 region was held fixed while the patch layout changed, and the change
in the values Cellpose received for the same pixels was recorded (mean per-pixel
spread, worst of patch size changed and patch origin shifted). A gaussian-only
control measured 0.0000.

| Nuclei per mm² | DAPI-positive area | sopa preprocessing | Cellpose preprocessing |
| -------------- | ------------------ | ------------------ | ---------------------- |
| 2,250          | 10.1%              | 0.1211             | 0.0265                 |
| 2,932          | 14.6%              | 0.0651             | 0.0052                 |
| 4,647          | 24.2%              | 0.0817             | 0.0075                 |
| 12,006         | 66.2%              | 0.0758             | 0.0553                 |

### Patch size

The same 1024 × 1024 region was segmented as part of a 1500-pixel patch and as
part of a 3000-pixel patch, and the cells were matched by IoU. This measures how
many cells keep the same outline when only `patch_width_pixel` changes:

| Model      | sopa: IoU ≥ 0.90 | Cellpose: IoU ≥ 0.90 | sopa mean IoU | Cellpose mean IoU |
| ---------- | ---------------- | -------------------- | ------------- | ----------------- |
| `cpsam_v2` | 103/267 (38.6%)  | 267/326 (81.9%)      | 0.869         | 0.928             |
| `cpsam`    | 155/359 (43.2%)  | 271/323 (83.9%)      | 0.860         | 0.929             |
| `cpdino`   | 74/201 (36.8%)   | 317/366 (86.6%)      | 0.855         | 0.936             |

With the default (Cellpose) preprocessing, `patch_width_pixel` mainly changes
the number and length of tasks.

### Effect on results

Turning sopa's preprocessing off changes results. On dense patches:

| Model      | Cell count   | Median cell area |
| ---------- | ------------ | ---------------- |
| `cpsam_v2` | +8% to +35%  | −23% to −43%     |
| `cpsam`    | −3% to −10%  | −25% to −37%     |
| `cpdino`   | +28% to +54% | −25% to −42%     |

Median cell area falls for every model, so per-cell intensities change. Do not
pool results from pipeline versions before 0.5.0 with later results in one
analysis.

The cells are now about 24-27 pixels in equivalent diameter against the
configured 30. This is inside the range Cellpose-SAM was trained on, so
`cellpose_diameter` does not need changing.

CLAHE also stretches an all-zero patch to the full intensity range, which skips
Cellpose's check for empty patches (a 1st-99th percentile range below 0.001).
On a test image, 13 of 84 patches at `patch_width_pixel = 3000` were all zero.
`cpsam` found 0 cells in such patches with and without CLAHE, so this cost time
but did not create cells.

### Visual check

The defaults were run on a COMET image region (1683 × 2349 pixels at
0.28 µm/pixel, DAPI plus four membrane markers combined with `max`) and
inspected in QuPath:

| Method              | Cells | Median area            | Median equivalent diameter |
| ------------------- | ----- | ---------------------- | -------------------------- |
| Cellpose `cpsam_v2` | 1,742 | 864 pixels² (68 µm²)   | 33.2 pixels (9.3 µm)       |
| Cellpose `cpdino`   | 1,751 | 916 pixels² (72 µm²)   | 34.2 pixels (9.6 µm)       |
| CellSAM             | 1,648 | 1,129 pixels² (89 µm²) | 37.9 pixels (10.6 µm)      |

The three methods agree within 6% on cell count. CellSAM draws larger outlines
around membrane signal.

### Limits of these measurements

- They measure reproducibility, not accuracy. No result was scored against
  annotations.
- The normalisation measurements used the DAPI channel only. Whole-cell
  segmentation uses membrane channels, which are dimmer, and was not measured.
  CLAHE may help dim membrane channels. If your panel has weak membrane markers,
  compare both settings on your data.
- CellSAM also applies CLAHE to each block, but its model was trained with
  CLAHE and uses a fixed 128-pixel kernel. Cellpose was not trained with CLAHE.

## Tile-local normalisation

`cellpose_tile_norm_blocksize` makes Cellpose compute the 1st and 99th
percentiles per block inside the patch instead of once per patch. The default,
0, is Cellpose's own default (whole patch).

If you set it, use a block much larger than a cell. At 0.28 µm/pixel with a
30-pixel diameter, 512 pixels is 143 µm (about 17 nuclei across). 128 pixels is
36 µm (about 4 nuclei), and a block that size in sparse tissue can contain no
nuclei. In testing at 128 pixels, the smallest per-block percentile range was
5.0 against 1670 for the whole patch, noise was stretched to the full range,
and patch dependence became 3 to 24 times worse. Using it together with CLAHE
also made results worse.

## GPU

Cellpose runs on a GPU (`SOPA_SEGMENTATIONCELLPOSE` has the `process_gpu`
label). On a machine without a GPU, override that label's settings in a `-c`
config file. Cellpose then runs on CPU, more slowly.
