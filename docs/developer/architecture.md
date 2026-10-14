# Architecture

sp_segment follows the nf-core template (version 3.3.2). This page maps the
code.

## Layout

| Path                       | Contents                                                            |
| -------------------------- | ------------------------------------------------------------------- |
| `main.nf`                  | Entry point. Runs initialisation, `SP_SEGMENT`, completion.         |
| `workflows/sp_segment.nf`  | The main workflow. Splits samples by method and step.               |
| `subworkflows/local/`      | One subworkflow per segmentation path, plus utilities.              |
| `modules/local/`           | Pipeline processes.                                                 |
| `modules/nf-core/backsub/` | Background subtraction module from nf-core/modules.                 |
| `bin/`                     | Python scripts called by the modules.                               |
| `conf/base.config`         | Default resources per label.                                        |
| `conf/modules.config`      | Tool arguments (`ext.args`), publishing, GPU options.               |
| `conf/test*.config`        | Test profiles.                                                      |
| `nextflow.config`          | Parameters, container and resource profiles.                        |
| `nextflow_schema.json`     | Parameter schema. Drives validation, `--help`, Seqera and the docs. |
| `assets/schema_input.json` | Samplesheet schema.                                                 |
| `tests/`                   | Test samplesheets, synthetic data, pytest tests.                    |

## Startup checks

`PIPELINE_INITIALISATION`
(`subworkflows/local/utils_nfcore_sp_segment_pipeline/main.nf`) runs before any
task. It validates parameters against the schema, reads the samplesheet, and
stops with an error for:

- an unknown Cellpose model name or a missing custom model path
- a missing `cellpose_models_dir`, or one under `/opt`
- `cellpose_model_type` being set
- `container_cache_dir` set while the Apptainer layer cache is under `$HOME`

## Main workflow

`SP_SEGMENT` reads each samplesheet row and sends it to one path:

| Row                                   | Subworkflow                |
| ------------------------------------- | -------------------------- |
| `run_backsub` only                    | `BACKGROUNDSUBTRACT`       |
| Mesmer, with background subtraction   | `MESMER_SEGMENT_WBACKSUB`  |
| Mesmer                                | `MESMER_SEGMENT`           |
| Cellpose, with background subtraction | `SOPA_SEGMENT_WBACKSUB`    |
| Cellpose                              | `SOPA_SEGMENT`             |
| CellSAM, with background subtraction  | `CELLSAM_SEGMENT_WBACKSUB` |
| CellSAM                               | `CELLSAM_SEGMENT`          |

Each `*_WBACKSUB` subworkflow runs `BACKGROUNDSUBTRACT` and passes the
background-subtracted image to the matching segmentation subworkflow.

Two steps run once at the top level instead of inside each path:

- `CELLPOSEMODEL` downloads the Cellpose weights once per run. Every Cellpose
  patch task receives the same copy. It is skipped when `cellpose_models_dir`
  is set or `cellpose_pretrained_model` is a file.
- `KRONOS2EMBEDDINGS` runs on the cell measurement output of every path.

## Segmentation subworkflows

All three run `COMBINECHANNELS` first. It writes a two-channel image: the
nuclear channel and the membrane channels combined into one.

- `SOPA_SEGMENT` runs `SOPA_SEGMENT_COMPARTMENT` for the whole-cell mask, and
  for the nuclear mask when `use_whole_cell_only = false`. Each run does
  `SOPA_CONVERT` (image to SpatialData zarr), `SOPA_PATCHIFYIMAGE`,
  `SOPA_CELLPOSE` (one `SOPA_SEGMENTATIONCELLPOSE` task per patch, then
  `SOPA_RESOLVECELLPOSE`) and `PARQUETTOTIFF`.
- `MESMER_SEGMENT` runs `MESMERSEGMENT` as `MESMERWC` and `MESMERNUC`.
- `CELLSAM_SEGMENT` runs `CELLSAMSEGMENT` as `CELLSAMWC` and `CELLSAMNUC`.

Each then runs `SMOOTHMASKS` (when `smooth_masks = true`), `CELLMEASUREMENT`
and `SEGMENTATIONREPORT` (when `generate_report = true`).

## Modules and scripts

| Module                    | Script or tool                                       | Container                                                                           |
| ------------------------- | ---------------------------------------------------- | ----------------------------------------------------------------------------------- |
| `EXTRACTMARKERS`          | `bin/extract_markers.py`                             | `community.wave.seqera.io/library/tifffile_pandas_typer`                            |
| `BACKSUB`                 | Background_subtraction (Schapiro lab)                | `ghcr.io/schapirolabor/background_subtraction:v0.4.1`                               |
| `COMBINECHANNELS`         | `bin/combine_channels.py`                            | `community.wave.seqera.io/library/tifffile_xarray_numpy_typer`                      |
| `SOPA_*`, `CELLPOSEMODEL` | sopa 2.2.6 and cellpose 4.2.1.1                      | `community.wave.seqera.io/library/python_pip_sopacellpose_cellpose`                 |
| `PARQUETTOTIFF`           | `bin/parquet_to_tiff.py`                             | `community.wave.seqera.io/library/tifffile_pyarrow_rasterio_shapely_pruned`         |
| `MESMERSEGMENT`           | MesmerSegmentation CLI                               | `ghcr.io/wehi-soda-hub/mesmersegmentation:0.3.1`                                    |
| `CELLSAMSEGMENT`          | `bin/cellsam_segment.py`                             | `community.wave.seqera.io/library/python_tifffile_scikit-image_scikit-learn_pruned` |
| `SMOOTHMASKS`             | `bin/smooth_masks.py`                                | `community.wave.seqera.io/library/python_tifffile_numpy_scipy_pruned`               |
| `CELLMEASUREMENT`         | cellmeasurement-py                                   | `ghcr.io/wehi-soda-hub/cellmeasurement-py:0.2.2`                                    |
| `KRONOS2EMBEDDINGS`       | `bin/kronos2_embeddings.py`, `bin/kronos2_common.py` | `community.wave.seqera.io/library/kronos2embeddings`                                |
| `SEGMENTATIONREPORT`      | Quarto template with spatialVis                      | `ghcr.io/wehi-soda-hub/spatialvis:0.2.0`                                            |

The Wave images are pinned by tag. The exact image for each process is in its `main.nf`.

## Publishing

`conf/modules.config` publishes every process's outputs to
`<outdir>/<name>/`, where `<name>` is the last part of the process name before
any `_`, in lower case (for example `SOPA_CONVERT` to `sopa/`, `MESMERWC` to
`mesmerwc/`). `CELLPOSEMODEL` is not published. `KRONOS2EMBEDDINGS` writes its
GeoJSON to `cellmeasurement/` and its report to `kronos2embeddings/`.
