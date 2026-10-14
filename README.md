<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/wehi-soda-hub-sp_segment_logo_dark.png">
    <img alt="WEHI-SODA-Hub/sp_segment" src="assets/wehi-soda-hub-sp_segment_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/WEHI-SODA-Hub/sp_segment/actions/workflows/ci.yml/badge.svg)](https://github.com/WEHI-SODA-Hub/sp_segment/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/WEHI-SODA-Hub/sp_segment/actions/workflows/linting.yml/badge.svg)](https://github.com/WEHI-SODA-Hub/sp_segment/actions/workflows/linting.yml)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](https://wehi-soda-hub.github.io/sp_segment/)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.17103183-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.17103183)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.04.2-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with docker](http://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/WEHI-SODA-Hub/sp_segment)

## Introduction

**WEHI-SODA-Hub/sp_segment** is a Nextflow pipeline that segments cells in
multiplex spatial proteomics images from COMET, MIBI and OPAL. For each image it
can subtract background (COMET only), segment whole cells and optionally nuclei
with Cellpose, Mesmer or CellSAM, measure shape and per-channel intensity for
each cell, add KRONOS2 embeddings, and write a QC report. The cells and their
measurements are written to a GeoJSON file that opens in QuPath.

**Documentation: <https://wehi-soda-hub.github.io/sp_segment/>** (also readable
in [`docs/`](docs/README.md)).

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

Dashed boxes are optional steps. Nuclear segmentation is optional for every
method and is off by default.

## Quick start

Requirements: Nextflow 24.04.2 or later, and Apptainer, Singularity or Docker.
Cellpose, CellSAM and KRONOS2 need a GPU.

Test your setup:

```bash
nextflow run WEHI-SODA-Hub/sp_segment -profile test,singularity --outdir results
```

Write a samplesheet with one row per image:

```csv
sample,run_backsub,run_cellpose,tiff,nuclear_channel,membrane_channels
slide1,true,true,/data/slide1.ome.tiff,DAPI,CD45:Pan-Cytokeratin
```

Run the pipeline:

```bash
nextflow run WEHI-SODA-Hub/sp_segment \
    -r 0.5.0 \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results
```

Set parameters on the command line or with `-params-file`, not in a `-c`
config file. Mesmer needs a DeepCell token:
`nextflow secrets set DEEPCELL_ACCESS_TOKEN <token>`.

Next steps:

- [Install and test](docs/getting-started.md)
- [Running at WEHI](docs/wehi.md) (Milton and Seqera Platform)
- [Samplesheet](docs/reference/samplesheet.md) and
  [parameters](docs/reference/parameters.md)
- [Outputs](docs/output.md)
- [Troubleshooting](docs/troubleshooting.md)

## Tools

- [Background_subtraction](https://github.com/SchapiroLabor/Background_subtraction): COMET background subtraction
- [Cellpose](https://www.cellpose.org) through [sopa](https://github.com/gustaveroussy/sopa): patch-based segmentation
- [MesmerSegmentation](https://github.com/WEHI-SODA-Hub/mesmersegmentation): Mesmer segmentation
- [CellSAM](https://github.com/vanvalenlab/cellSAM): foundation-model segmentation
- [cellmeasurement-py](https://github.com/WEHI-SODA-Hub/cellmeasurement-py): nucleus-to-cell matching and measurements
- [KRONOS2](https://huggingface.co/MahmoodLab/KRONOS2): per-cell embeddings
- [spatialVis](https://github.com/WEHI-SODA-Hub/spatialVis): segmentation report plots

## Credits

WEHI-SODA-Hub/sp_segment was originally written by the WEHI SODA Hub.

We thank the following people for their extensive assistance in the development
of this pipeline:

- Michael McKay (@mikemcka)
- Emma Watson

## Contributions and support

See the [contributing guide](docs/developer/contributing.md). Report problems in
the [issue tracker](https://github.com/WEHI-SODA-Hub/sp_segment/issues).

## Citations

If you use WEHI-SODA-Hub/sp_segment for your analysis, please cite it using the
following DOI:
[10.5281/zenodo.17103183](https://doi.org/10.5281/zenodo.17103183)

References for the tools used by the pipeline are in
[`CITATIONS.md`](CITATIONS.md).

This pipeline was created using the `nf-core` template. You can cite the
`nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
