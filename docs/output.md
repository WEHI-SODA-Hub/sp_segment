# Outputs

All paths are relative to `--outdir`. `<sample>` is the `sample` value from the
samplesheet. Which directories appear depends on the steps each sample runs.

## Summary

| Directory                   | Written when             | Main files                                          |
| --------------------------- | ------------------------ | --------------------------------------------------- |
| `extractmarkers/`           | `run_backsub = true`     | `<sample>_markers.csv`                              |
| `backsub/`                  | `run_backsub = true`     | `<sample>.ome.tif`, markers CSV                     |
| `combinechannels/`          | any segmentation         | `<sample>_combined_channels.tiff`                   |
| `sopa/`                     | Cellpose                 | `<sample>.zarr`                                     |
| `parquettotiff/`            | Cellpose                 | `<sample>_whole-cell.tiff`, `<sample>_nuclear.tiff` |
| `mesmerwc/`, `mesmernuc/`   | Mesmer                   | `<sample>_whole-cell.tiff`, `<sample>_nuclear.tiff` |
| `cellsamwc/`, `cellsamnuc/` | CellSAM                  | `<sample>_whole-cell.tiff`, `<sample>_nuclear.tiff` |
| `smoothmasks/`              | `smooth_masks = true`    | `*_smoothed.tiff`                                   |
| `cellmeasurement/`          | any segmentation         | `<sample>.geojson.gz`, `<sample>_mask.tiff`         |
| `kronos2embeddings/`        | `enable_kronos = true`   | `<sample>_marker_report.txt`                        |
| `segmentationreport/`       | `generate_report = true` | `<sample>/<sample>.html`                            |
| `pipeline_info/`            | always                   | Nextflow reports, software versions                 |

Nuclear masks (`*_nuclear.tiff`, `mesmernuc/`, `cellsamnuc/`) are written only
when `use_whole_cell_only = false`.

## Background subtraction

COMET only.

- `extractmarkers/<sample>_markers.csv`: one row per channel with the marker
  name, its background channel and its exposure time, read from the OME-XML.
- `backsub/<sample>.ome.tif`: the background-subtracted image. Channels in
  `remove_markers` (default `TRITC,Cy5`) are removed.

When a sample also runs segmentation, the segmentation uses this image.

## Combined channels

- `combinechannels/<sample>_combined_channels.tiff`: a two-channel image with
  the nuclear channel and the membrane channels combined into one with
  `combine_method`. This is the image each segmentation method reads.

## Segmentation masks

Each mask is a TIFF in which every cell has its own integer label and the
background is 0. Masks are `uint32`, so images with more than 65,535 cells are
supported.

- Cellpose: `parquettotiff/<sample>_whole-cell.tiff` and
  `parquettotiff/<sample>_nuclear.tiff`. `sopa/<sample>.zarr` holds the
  SpatialData store with the image patches and the Cellpose cell boundaries.
- Mesmer: `mesmerwc/<sample>_whole-cell.tiff` and
  `mesmernuc/<sample>_nuclear.tiff`.
- CellSAM: `cellsamwc/<sample>_whole-cell.tiff` and
  `cellsamnuc/<sample>_nuclear.tiff`. CellSAM fails the task if the mask
  contains no cells.

## Smoothed masks

When `smooth_masks = true`, `smoothmasks/` holds a smoothed copy of each mask,
named after the input mask with `_smoothed.tiff` added, for example
`<sample>_whole-cell_smoothed.tiff`. Cell measurement then uses the smoothed
masks. See `smooth_method` in [Parameters](reference/parameters.md#mask-smoothing-options).

## Cell measurement

- `cellmeasurement/<sample>.geojson.gz`: one GeoJSON feature per cell, with
  the cell outline, the nucleus outline when nuclei were segmented, and the
  measurements. With `gzip_geojson = false` the file is `<sample>.geojson`.
- `cellmeasurement/<sample>_mask.tiff`: the final cell label mask, rasterised
  from the cells in the GeoJSON.

### Measurement names

Measurements use QuPath-style names: `<channel>: <compartment>: <statistic>`.

| Measurement                               | Example                                  | Written when                   |
| ----------------------------------------- | ---------------------------------------- | ------------------------------ |
| Intensity per compartment                 | `CD45: Cell: Mean`                       | `enable_measurements = true`   |
| Intensity percentile                      | `CD45: Cell: Percentile: 95`             | `percentiles` is set           |
| Erosion ring (1-5, boundary inward)       | `CD45: Cell: ErosionBin_1: Mean`         | `erosion_steps = true`         |
| Expansion ring (1-5, up to 20 µm outward) | `CD45: Cell: ExpansionBin_3: Median`     | `expansion_steps = true`       |
| 20 µm environment zone                    | `CD45: Cell: Environment_20um: Std.Dev.` | `environment_expansion = true` |
| Shape                                     | `Cell: Area µm^2`, `Cell: Solidity`      | `enable_measurements = true`   |
| Neighbourhood mean                        | `Neighbours: Mean: CD45: Cell: Mean`     | `neighbors` > 0                |
| KRONOS2 embedding (768 values)            | `kronos_emb_0` to `kronos_emb_767`       | `enable_kronos = true`         |

- Intensity statistics are `Mean`, `Median`, `Min`, `Max` and `Std.Dev.`.
  Erosion and expansion rings have `Mean` and `Median` only.
- Compartments are `Cell`, `Nucleus`, `Cytoplasm` and `Membrane`. `Nucleus`
  measurements are written only when `use_whole_cell_only = false`.
- Shape measurements are `Area µm^2`, `Circularity`, `Length µm`,
  `Max diameter µm`, `Min diameter µm` and `Solidity`, for `Cell` and, with
  nuclei, `Nucleus`. With nuclei, `Nucleus/Cell area ratio` is also written.
- Distances and areas in µm use `pixel_size_microns` (default 0.28).

## KRONOS2 embeddings

When `enable_kronos = true`:

- `cellmeasurement/<sample>.geojson.gz` is replaced by a copy with 768
  `kronos_emb_*` measurements added to every cell. Embedding row _i_ belongs to
  cell feature _i_. Running KRONOS2 again on the same file replaces the
  `kronos_emb_*` values.
- `kronos2embeddings/<sample>_marker_report.txt` lists, for each channel, the
  name given to KRONOS2, any mapping applied, any channel left out with
  `kronos_exclude_markers`, any name not in the KRONOS2 vocabulary, and any
  `kronos_marker_mapping` entry that matches no channel. The
  first line gives the number of channels KRONOS2 used. The report is written
  even when the run fails on an unmatched marker.

768 numbers per cell make the GeoJSON much larger. Keep `gzip_geojson = true`
(the default).

## Segmentation report

When `generate_report = true`:

- `segmentationreport/<sample>/<sample>.html`: a QC report for the sample. It is
  most complete with `enable_measurements = true`.

## Pipeline information

`pipeline_info/` contains:

- Nextflow reports: `execution_report_<timestamp>.html`,
  `execution_timeline_<timestamp>.html`, `execution_trace_<timestamp>.txt` and
  `pipeline_dag_<timestamp>.html`.
- `nf_core_pipeline_software_versions.yml`: the version of every tool used.
- `params_<timestamp>.json`: the parameters of the run.
- `pipeline_report.html` and `pipeline_report.txt`, only when `--email` or
  `--email_on_fail` is set.

The Nextflow reports show the run time, CPU and memory use of every task, which
helps when choosing resources.
