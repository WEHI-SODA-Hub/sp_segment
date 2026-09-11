# WEHI-SODA-Hub/sp_segment: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

Depending on the mode(s) run, the pipeline can produce background subtracted TIFFs, TIFF masks containing nuclear and whole-cell segmentations, and resolved whole-cell and nuclear segmentations in GeoJSON format.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- Background subtraction (COMET only) - generates a background subtracted TIFF
- Segmentation via Cellpose (COMET), Mesmer (COMET/MIBI), or CellSAM
  (COMET/MIBI) for nuclear and whole cell masks
- Cell measurement
  - generates a GeoJSON file with consolidated whole-cell/nuclear segmentations
  - calculates cell compartment measurements and channel intensities
- Optional KRONOS2 embeddings merged into the cellmeasurement GeoJSON
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Pipeline outputs

<details markdown="1">
<summary>Output files</summary>

#### Background subtraction (COMET only)

- `extractmarkers/sample.csv` -- marker names, background per channel and
  exposure time
- `backsub/sample.tiff` -- bakground subtracted tiff image.

#### Cellpose segmentation

- `sopa/sample.zarr` -- SpatialData converted image containing segmentations.
- `parquettotiffwholecell/sample_whole-cell.tiff` -- tiff label masks from cellpose segmentation for
  whole cell segmentation.
- `parquettotiffnuclear/sample_nuclear.tiff` -- tiff label masks from cellpose segmentation for
  nuclear segmentation.

#### AVITI24 cytoprofiling segmentation

AVITI produces tile-level viewer outputs plus stitched per-well outputs:

- `avitisegmentation/sample/CellSegmentation/WellA1/<tile>_Cell.tif` -- per-tile whole-cell mask (uint16, instance-labeled) published in the Cytocanvas layout. Produced by the Cellpose v4 SAM path, or by the Cellpose 3.x membrane-model path when the samplesheet row sets `membrane_model`; both write to this same path.
- `avitisegmentation/sample/CellSegmentation/WellA1/<tile>_Nuclear.tif` -- per-tile nuclear mask (uint8, binary 0/1 presence) published in the Cytocanvas layout, matching Elembio's own `cells2stats`/Cytocanvas convention.
- `avitistitched/sample/WellA1/sample__WellA1_cell_stitched.tif` -- stitched whole-cell mask for one well.
- `avitistitched/sample/WellA1/sample__WellA1_nuclear_stitched.tif` -- stitched **instance-labeled** nuclear mask for one well (not the binary Nuclear.tif convention above -- see note below).
- `avitistitched/sample/WellA1/sample__WellA1_image_stitched.tif` -- stitched multi-channel image for downstream cellmeasurement/KRONOS.

The tile names preserve the source AVITI tile basename, and the stitched
outputs are the ones passed into the existing `CELLMEASUREMENT`,
`SEGMENTATIONREPORT`, and `KRONOS2EMBEDDINGS` wiring.

Nuclear.tif's binary presence-mask convention (required by Elembio's
`cells2stats`/Cytocanvas) can't carry per-nucleus instance IDs, but
`CELLMEASUREMENT` needs them to match nuclei to whole cells by centroid. To
support both, nuclear segmentation also writes an internal, unpublished
`<tile>_Nuclear_label.tif` (uint16, instance-labeled) alongside the published
binary `Nuclear.tif`; it is this label mask -- not the binary one -- that
gets stitched into `sample__WellA1_nuclear_stitched.tif` above.

##### Plate-level assembly

When `aviti_plate_assembly` is enabled (the default), every well's stitched
output for a sample is additionally combined into one plate-level view,
alongside -- not instead of -- the per-well outputs above:

- `avitiplate/sample/sample_plate.ome.tif` -- a single pyramidal, tiled
  BigTIFF OME-TIFF laying out every well on the plate grid (letter = column,
  number = row: A1 top-left, A2 directly below A1, B1 the next column
  across). Opens directly in QuPath with smooth zoom from the whole plate
  down to individual cells; the OME `PhysicalSizeX/Y` reflects the run's own
  pixel size (see below), so QuPath's scale bar is correct.
- `avitiplate/sample/sample_plate.geojson.gz` -- every cell from every well's
  cellmeasurement annotations, translated into plate pixel coordinates, plus
  one labelled rectangle annotation per well. Import this into QuPath
  alongside the plate image to overlay cell objects and measurements across
  the whole plate. (`.gz` unless `--gzip_geojson false`.)
- `avitiplate/sample/sample_plate_layout.csv` -- the well placements
  (`well,col,row,x0,y0,width,height`) used to build both artefacts above.
- `avitiplate/sample/sample_plate_overview.tif` -- a small, single-resolution
  overview image (the smallest pyramid level) for a quick-look outside QuPath.

The plate image can be tens of GB uncompressed; `avitiplate/` is published by
symlink rather than copy for that reason. See
[docs/usage.md](usage.md#aviti24-cytoprofiling-segmentation) for the
plate-related parameters and the rule for when a per-well vs. plate-level
`SEGMENTATIONREPORT` is generated.

##### AVITI pixel size

AVITI tile discovery reads the run's own `ImageInfo.PixelSizeUm` from
`RunParameters.json` and uses it -- rather than the pipeline-wide
`pixel_size_microns` default -- for well/plate stitching gaps, the plate
image's OME `PhysicalSizeX/Y`, and `CELLMEASUREMENT`'s µm-based measurements.
`pixel_size_microns` is only a fallback for a run directory whose
`RunParameters.json` predates this field.

#### Mesmer segmentation

- `mesmerwc/sample_whole-cell.tiff` -- whole-cell label mask generated by Mesmer.
- `mesmernuc/sample_nuclear.tiff` -- nuclear label mask generated by Mesmer.

#### CellSAM segmentation

- `cellsamwc/sample_whole-cell.tiff` -- whole-cell label mask generated by
  CellSAM.
- `cellsamnuc/sample_nuclear.tiff` -- nuclear label mask generated by CellSAM
  (not produced when `use_whole_cell_only=true`).

#### Cell measurement

- `cellmeasurement/sample.geojson` -- resolved whole-cell and nuclear
  segmentations, optionally containing measurements and intensity values per
  cell, compatible with QuPath.

#### KRONOS2 embeddings

When `enable_kronos=true`:

- `cellmeasurement/sample.geojson` (or `.geojson.gz`) -- same output path as the
  cellmeasurement annotations, with 768 `kronos_emb_*` measurements added to
  every cell. This is the **only** embedding artefact: the separate CSV that
  earlier versions wrote held the same vectors twice, so it has been dropped.
- `kronos2embeddings/sample_marker_report.txt` -- per-channel record of the name
  given to KRONOS2, any mappings applied, any channels withheld with
  `kronos_exclude_markers`, and any markers outside the model's vocabulary. The
  header gives the embedded/not-shown split, so the number of channels the model
  actually saw is on the first line. It is written even when the run then fails
  on an unmatched marker, because that is when it is most useful.

Embedding row _i_ corresponds to cell feature _i_ in the input GeoJSON, so the
join is exact by construction. Re-running clears any previous `kronos_emb_*`
keys first, so embeddings are replaced rather than accumulated.

Note that 768 float measurements per cell is a substantial addition to the
GeoJSON. `gzip_geojson` (on by default) keeps this manageable.

#### Segmentation report

- `segmentationreport/sample/sample.html` -- html file for visualising report.

#### Data provenance

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
