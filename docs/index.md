# sp_segment

![sp_segment](images/wehi-soda-hub-sp_segment_logo_light.png#only-light)
![sp_segment](images/wehi-soda-hub-sp_segment_logo_dark.png#only-dark)

sp_segment is a Nextflow pipeline that segments cells in multiplex spatial
proteomics images from COMET, MIBI and OPAL. It writes one GeoJSON file per
image that contains every cell outline and its measurements. You can open the
GeoJSON in QuPath.

For each image, the pipeline can:

- subtract background (COMET only)
- segment whole cells, and optionally nuclei, with Cellpose, Mesmer or CellSAM
- smooth the label masks
- measure shape and per-channel intensity for each cell compartment
- add a 768-value KRONOS2 embedding to every cell
- write an HTML QC report

## Workflow

Dashed boxes are optional steps. Nuclear segmentation is optional for every
method and is off by default (`use_whole_cell_only = true`).

```mermaid
flowchart TD
  input["Multiplex TIFF"] --> extract["Extract markers"]
  extract --> bgsub["Background subtraction"]
  bgsub --> combine["Combine membrane channels"]
  input --> combine
  combine --> choose{"Segmentation method"}

  choose -- Cellpose --> convert["sopa convert"]
  convert --> patchify["sopa patchify"]
  patchify --> cellpose["Cellpose on each patch"]
  cellpose --> resolve["sopa resolve"]
  resolve --> totiff["Parquet to TIFF"]
  totiff --> smooth["Smooth masks"]

  choose -- Mesmer --> mesmer["Mesmer"]
  mesmer --> smooth

  choose -- CellSAM --> cellsam["CellSAM"]
  cellsam --> smooth

  smooth --> measure["Cell measurement"]
  measure --> geojson["GeoJSON"]
  measure --> kronos["KRONOS2 embeddings"]
  kronos --> kgeojson["GeoJSON with embeddings"]
  kronos --> markers["Marker report"]
  measure --> report["Segmentation report"]

  style extract stroke-dasharray: 5 5
  style bgsub stroke-dasharray: 5 5
  style smooth stroke-dasharray: 5 5
  style kronos stroke-dasharray: 5 5
  style report stroke-dasharray: 5 5
```

## Supported inputs

| Input                | Background subtraction | Cellpose | Mesmer | CellSAM |
| -------------------- | ---------------------- | -------- | ------ | ------- |
| COMET OME-TIFF       | yes                    | yes      | yes    | yes     |
| OPAL OME-TIFF        | no                     | yes      | yes    | yes     |
| MIBI multi-page TIFF | no                     | no       | yes    | yes     |

Cellpose requires an OME-TIFF, because sopa converts the image with
`technology = ome_tif`. See
[Choose a segmentation method](how-to/choose-segmenter.md).

## Where to start

- New to the pipeline: [Install and test](getting-started.md).
- At WEHI: [Running at WEHI](wehi.md).
- Writing a samplesheet: [Samplesheet](reference/samplesheet.md).
- Looking up a parameter: [Parameters](reference/parameters.md).
- A run failed: [Troubleshooting](troubleshooting.md).

## Citation

If you use sp_segment, cite
[10.5281/zenodo.17103183](https://doi.org/10.5281/zenodo.17103183) and the
tools listed in [Citations](citations.md).
