# WEHI-SODA-Hub/sp_segment: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Changed

- **:warning: `use_whole_cell_only` now defaults to `true`.** Nuclei are not
  segmented or measured unless `--use_whole_cell_only false` is set. Runs that
  need nuclear masks or nuclear measurements must set it.
- **:warning: `cellpose_min_area` now defaults to `100` (was `200`).** Cellpose
  keeps cells from 100 to 199 pixels² that were removed before. Set
  `--cellpose_min_area 200` to match earlier results. `cellsam_min_area` still
  defaults to `200`.
- The **Mask smoothing options** and **Cellpose preprocessing options**
  parameters are hidden in `--help` and the Seqera Platform launch form. Show
  them with `--show_hidden`. They can still be set, and their defaults are
  unchanged.

### Added

- `cellsam_remove_border_cells` (default `true`). Set it to `false` to keep
  CellSAM cells that touch the image edge. Border cells were previously always
  removed.

### Fixed

- **:warning: Background subtraction could use another sample's markers
  file.** With two or more `run_backsub` samples in one run, an image could be
  corrected with another slide's exposure times and background channels.
  **Re-run background subtraction for runs with more than one backsub
  sample.**
- **:warning: `--combine_method prod` saturated the combined membrane
  channel.** Every pixel whose product exceeded 65535 was set to 65535. The
  product is now rescaled so its maximum is 65535. The default, `max`, was not
  affected.
- `extract_markers` could give a marker another channel's background channel
  or exposure time when the OME-XML lists ChannelPriv elements in a different
  order from the channels, or stores one channel per Image. Output for
  standard COMET files is unchanged.
- The samplesheet is checked before any task runs. Sample names must be
  unique; a row can set at most one of `run_mesmer`, `run_cellpose` and
  `run_cellsam`; a row must run something; segmentation rows need
  `nuclear_channel` and `membrane_channels`. These cases previously overwrote
  outputs, dropped samples or crashed.
- Channel names can contain any character except `"`, `'`, `:` and `,`, and
  names with spaces are passed correctly to every tool.
- `combine_channels` reads images that store one channel per OME Image (COMET
  `_fixed` exports), panels with a repeated channel name, and TIFFs with a
  plain JSON description. It reads only the nuclear and membrane channels and
  writes valid OME-XML.
- CellSAM fails when it finds no cells instead of writing an empty mask. The
  whole-cell and nuclear tasks no longer race on the shared model download.
- An invalid `kronos_marker_mapping` (not an existing file and not valid JSON)
  fails the run. A mapping file is staged into the container, so relative
  paths work. Mapping entries that match no channel are listed in the marker
  report. `enable_kronos` without `kronos_model_path` fails at startup.
- `software_versions.yml` lists every tool. CellSAM and Mesmer versions are
  read from the installed packages instead of being hardcoded.
- Cellpose resolve and measurement for a sample start when that sample's
  patches finish, instead of after every sample's patches.
- GPU tasks get `--nv` under Apptainer or Singularity with any profile. sopa
  Cellpose previously ran on CPU unless a `small`, `medium`, `large` or
  `wehi_*` profile was used.
- The test profiles read their samplesheets from the pipeline directory.
- KRONOS2, Mesmer and CellSAM no longer need a writable home directory inside
  the container. They failed with `Read-only file system` when the home
  directory was not mounted. Without `deepcell_cache_dir`, Mesmer and CellSAM
  now keep their model files in the task directory in that case.
- `deepcell_cache_dir` is created if it does not exist. Previously every
  Mesmer and CellSAM task failed to start until it was created by hand.

### Documentation

- Added a documentation site built with MkDocs Material, deployed to GitHub
  Pages with versions per release.
- The parameter and samplesheet references are generated from the JSON schemas
  by `docs/scripts/generate_reference.py`.
- `nextflow_schema.json` now lists the default of every parameter that has one
  in `nextflow.config`, and its help text is corrected.
- Removed the KRONOS2 migration working notes and the unused
  `docs/examples/additional_markers.csv`.

## v0.5.0 - 2026-09-07

### Changed

- **:warning: Label masks are always `uint32` (were `uint16`).** Masks with
  more than 65,535 cells were silently truncated: every cell past 65,535 was
  merged into low label IDs. **Regenerate any mask with more than 65,535 cells
  made by an earlier version.** Masks with fewer cells are unchanged.
  `cellsam_segment.py` fails if the cell count exceeds `uint32`.
- Intermediate label masks (`PARQUETTOTIFF`, `CELLSAMSEGMENT`) are written with
  deflate compression. A dense whole-slide mask fell from 2.87 GB to 31 MB. The
  smoothed mask is still uncompressed.
- `neighbors` defaults to `0` (disabled). It was `5`.
- Resource profiles: `small`, `medium` and `large` are unchanged. New WEHI
  Milton profiles `wehi_small`, `wehi_med` and `wehi_large` request most of a
  node, for samples with many cells.
- cellmeasurement-py updated to 0.2.2 (faster GeoJSON export overlap step) via
  0.2.1 (fixes channel detection on OPAL OME-TIFFs).
- **:warning: Cellpose preprocessing now matches Cellpose: sopa's gaussian
  filter and CLAHE are off.** New parameters `cellpose_clip_limit` (`0`),
  `cellpose_gaussian_sigma` (`0`), `cellpose_clahe_kernel_size` (unset) and
  `cellpose_tile_norm_blocksize` (`0`). Restore the old behaviour with
  `--cellpose_clip_limit 0.2 --cellpose_gaussian_sigma 1`. Results no longer
  depend strongly on `patch_width_pixel`: cells keeping the same outline
  (IoU ≥ 0.90) when the patch width changes rose from 37-43% to 82-87%.
  **Results change:** median cell area falls 23-43% and cell counts change by
  -10% to +54% depending on the model. Do not pool Cellpose results from before
  and after this version.
- The Singularity and Apptainer image pull timeout is 1 hour (was 20 minutes).
  CI caches the converted images.
- **:warning: `cellsam_min_area` defaults to `200` (was `0`)**, the same as
  `cellpose_min_area`. CellSAM objects smaller than 200 pixels² are removed.
  Set `--cellsam_min_area 0` for the old behaviour. Do not pool CellSAM results
  from before and after this version.
- The Cellpose parameter descriptions state which direction finds more cells.
- New parameter `container_cache_dir` puts converted container images on a
  path of your choice. If it is set while `APPTAINER_CACHEDIR` is under
  `$HOME`, the pipeline stops at startup. If the layer cache is under `$HOME`,
  it warns.
- `deepcell_cache_dir` is mounted into the Mesmer and CellSAM containers.
  Before, the weights were downloaded again in every task.
- Cellpose updated to 4.2.1.1 (sopa 2.2.6) and runs on the GPU. With the model
  fixed, the update does not change results.
- `cellpose_pretrained_model` selects the model and defaults to `cpsam_v2`
  (was unset). It accepts `cpsam_v2`, `cpsam`, `cpdino`, `cpdino-vitb` or a path
  to a custom model.
- `cellpose_model_type` is removed and setting it is an error. Cellpose 4
  ignored it, so runs that set `cyto3` used `cpsam`.
- Cellpose weights are downloaded once per run by the new `CELLPOSEMODEL`
  process, not once per patch task. This fixes HTTP 429 errors.
- A missing custom Cellpose model path or `cellpose_models_dir` stops the
  pipeline at startup. Before, Cellpose used its built-in weights without
  warning.
- The nf-test snapshot suite was repaired, including `sopa_patchifyimage - zarr`
  and the `backgroundsubtract` snapshot (which no longer records the output
  md5, because the file contains a new UUID on every run).
- **KRONOS1 is replaced by KRONOS2.** Embeddings have 768 values (were 384).
  **KRONOS1 and KRONOS2 embeddings cannot be compared.**
- KRONOS2 runs once at the top level, not inside each segmentation subworkflow.
- KRONOS2 patches contain only the target cell: pixels outside the cell are
  set to 0. Disable with `--kronos_isolate_cell false`.
- The KRONOS2 intensity divisor is derived from the image data type (uint8 255,
  uint16 65535, float 400). Float images were divided by 65535 before.
- The samplesheet `nuclear_channel` is passed to KRONOS2 as its DAPI channel, so
  DRAQ5 and Hoechst slides get DAPI statistics.
- `kronos_batch_size` defaults to `16` (was `32`), which reproduces the
  published values.
- The per-cell embeddings CSV is no longer written. Embeddings are in the
  cell measurement GeoJSON.

### Added

- `geometry_checkpoint_dir` and `geometry_batch_size` (default `2000`): a
  `CELLMEASUREMENT` task killed during polygon extraction resumes from the last
  completed batch. The directory must be outside the work directory.
- `kronos_exclude_markers`: channels KRONOS2 does not read, for channels that
  are not markers (for example autofluorescence). Other steps still use them.
  Excluding a channel changes every cell's embedding.
- The `cpdino` and `cpdino-vitb` Cellpose models work. The container includes
  Meta's `dinov3` code. See `CITATIONS.md` for the licence.
- `cellpose_models_dir`: a directory of pre-downloaded Cellpose weights.
- `CELLPOSE_MODELS_DIR` environment variable for nf-test runs
  (`tests/nextflow.config`).
- `bin/kronos2_common.py`, shared code for per-cell foundation-model
  embeddings.
- KRONOS2 container pinned by digest
  (`community.wave.seqera.io/library/kronos2embeddings:70590359503eec08`, torch
  2.6.0+cu124).
- `kronos_nuclear_marker`, `kronos_isolate_cell` and
  `kronos_allow_novel_defaults` parameters, and the `test_mesmer_kronos2`
  profile.
- KRONOS2 fails by default when a channel name is not in its vocabulary. Fix
  with `kronos_marker_mapping`, or accept default statistics with
  `kronos_allow_novel_defaults`.
- KRONOS2 cell footprints are drawn with Pillow instead of rasterio. Embeddings
  change by a cosine similarity of 0.99974.

### Removed

- `kronos_marker_metadata`: KRONOS2 has its own 288-marker vocabulary.
- `kronos_config_path`: KRONOS2 has no config file.
- `kronos_distance_threshold`: cells are joined to embeddings by GeoJSON feature
  order.

### Fixed

- **:warning: KRONOS2 read only the first channel of OME-TIFFs that store each
  channel as its own series** (for example COMET `_fixed` exports), and used it
  for every marker. **Regenerate embeddings made from such images by earlier
  versions.** The pipeline now reads each channel from its own series and fails
  if the channel count and the name count differ.
- A channel with neither an OME `Name` nor an `ID` shifted every later marker
  onto the wrong channel in KRONOS2. Names are now assigned by position.
- Stub runs of the Cellpose path failed with `For input string: ""`.
- The KRONOS2 module passed `--nv` under Docker. It now passes `--nv` or
  `--gpus all` according to the container engine.

## v0.4.0 - 2026-04-20

### Changed

- cellmeasurement: the Python implementation
  ([cellmeasurement-py](https://github.com/WEHI-SODA-Hub/cellmeasurement-py))
  replaces the Groovy one.
- cellmeasurement: measurements are computed on tiles in parallel. Tile size
  and overlap are configurable.
- cellmeasurement: erosion measurements use 5 equal-area rings named
  `<channel>: <compartment>: ErosionBin_<N>: Mean` and `Median`. The Groovy
  names (`Eroded_<N>px`) are no longer produced.
- cellmeasurement: `--erosion-steps` is a boolean and no longer takes pixel
  depths.
- The pipeline can run with whole-cell segmentation only.

### Added

- KRONOS embeddings (`kronosembeddings` module, 384 values per cell), disabled
  by default. Replaced by KRONOS2 in v0.5.0.
- CellSAM segmentation (`cellsam_segment`), with `cellsam_bbox_threshold`,
  `cellsam_block_size`, `cellsam_overlap`, `cellsam_iou_threshold` and
  `cellsam_use_wsi`.
- Mask smoothing (`smooth_masks`, `smooth_method`, `smooth_kernel_size`) with
  the `morphological` and `shapely` methods.
- `dist_threshold` (default `10.0`) and `downsample_factor` (default `1.0`)
  cell measurement parameters.

### Removed

- `gradle_cache_dir`.

## v0.3.0 - 2025-11-13

- Renamed to sp_segment.
- Fixed OME-XML processing.
- Added a background image preview and embedded report support to the
  segmentation report.
- Added percentile measurements.
- Updated the workflow diagram.

## v0.2.0 - 2025-09-12

- Added the segmentation report.

## v0.1.0 - 2025-08-29

Initial release of WEHI-SODA-Hub/sp_segment, created with the
[nf-core](https://nf-co.re/) template.
