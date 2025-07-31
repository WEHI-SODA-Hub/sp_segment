<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/wehi-soda-hub-spatialproteomics_logo_dark.png">
    <img alt="WEHI-SODA-Hub/spatialproteomics" src="docs/images/wehi-soda-hub-spatialproteomics_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/WEHI-SODA-Hub/spatialproteomics/actions/workflows/ci.yml/badge.svg)](https://github.com/WEHI-SODA-Hub/spatialproteomics/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/WEHI-SODA-Hub/spatialproteomics/actions/workflows/linting.yml/badge.svg)](https://github.com/WEHI-SODA-Hub/spatialproteomics/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/WEHI-SODA-Hub/spatialproteomics)

## Introduction

**WEHI-SODA-Hub/spatialproteomics** is a pipeline for running cell segmentation
on COMET and MIBI data. For COMET, background subtraction can be performed
followed by patched cellpose segmentation, or non-patched mesmer segmentation.
For MIBI, mesmer segmentation can be run. Whole-cell and nuclear segmentations
are run separately, and then consolidated into whole cells with nuclei with full
shape and intensity measurements per compartment. The output GeoJSON files can
be viewed in QuPath.

The pipeline performs these steps:

- Background subtraction (COMET only) -- generates a background subtracted TIFF
- Segmentation via Cellpose (COMET) or Mesmer (COMET/MIBI) for nuclear and whole
  cell
- Cell measurement
  - generates a GeoJSON file with consolidated whole-cell/nuclear segmentations
  - calculates cell compartment measurements and channel intensities

![spatialproteomics workflow](assets/spatialproteomics_workflow.png)

The pipeline uses the following tools:

- [Background_subtraction](https://github.com/SchapiroLabor/Background_subtraction)
  -- background subtraction tool for COMET.
- [MesmerSegmentation](https://github.com/WEHI-SODA-Hub/mesmersegmentation) -- a
  CLI for running Mesmer segmentation of MIBI and OME-XML TIFFs.
- [cellmeasurement](https://github.com/WEHI-SODA-Hub/cellmeasurement) -- a
  Groovy app that matches whole-cell segmentations with nuclei, and uses the
  QuPath API to calculate compartment measurements and intensities.
- [sopa](https://github.com/gustaveroussy/sopa) -- we use the sopa CLI tool to
  patchify images and perform cellpose segmentation.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` (to test cellpose segmentation) or `-profile test_mesmer` to test mesmer segmentation before running the workflow on actual data.

If you are running this pipeline from WEHI, it has been set up to run on [Seqera Platform](https://seqera.services.biocommons.org.au/).

> [!NOTE]
> If you don't have a .gradle directory in your home, make sure you create it with `mkdir $HOME/.gradle` before runnning the pipeline. You con't need to do this if you are running via WEHI's Seqera Platform mentioned above.

Usage will depend on your desired steps.

### Background subtraction

> [!NOTE]
> This step will only work with COMET OME-TIF files.

Prepare a sample sheet as follows:

`samplesheet.csv`:

```csv
sample,run_backsub,tiff
sample1,true,/path/to/sample1.tiff
sample2,true,/path/to/sample2.tiff
```

You may also prefer to use YAML for your samplesheet, either is supported:

`samplesheet.yml`:

```yaml
- sample: sample1
  run_backsub: backsub
  tiff: /path/to/sample1.tiff
- sample: sample2
  run_backsub: backsub
  tiff: /path/to/sample2.tiff
```

> [!WARNING]
> Please ensure that your image name and all directories in your file path do not contain spaces.

If you don't specify any segmentation algorithm to run (mesmer or cellpose), the
pipeline will run a background subtraction only.

Now, you can run the pipeline using:

```bash
nextflow run WEHI-SODA-Hub/spatialproteomics \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

### Mesmer segmentation

If you want to run Mesmer as your segmentation algorithm, you can specify a
config file like so:

```csv
sample,run_backsub,run_mesmer,tiff,nuclear_channel,membrane_channels
sample1,true,true,/path/to/sample1.tiff,DAPI,CD45:CD8
sample2,false,true,/path/to/sample2.tiff,DAPI,CD45
```

Nuclear channels only support one entry; membrane channels may have multiple
values separated by `:` characters. You can also set the following parameters,
either via CLI (e.g., `--mesmer_combine_method prod` or in a config file pass to
the workflow via `-c`:

| Parameter Name             | Description                                                                                                                             |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| mesmer_combine_method      | Method used to combine membrane channels (product or max).                                                                              |
| mesmer_segmentation_level  | Segmentation level (legacy parameter).                                                                                                  |
| mesmer_maxima_threshold    | Controls segmentation level directly in mesmer, (lower values = more cells, higher values = fewer cells).                               |
| mesmer_interior_threshold  | Controls how conservative model is in distinguishing cell from background (lower values = larger cells, higher values = smaller cells). |
| mesmer_maxima_smooth       | Controls what is considered a unique cell (lower values = more separate cells, higher values = fewer cells).                            |
| mesmer_min_nuclei_area     | Minimum area of nuclei to keep in square pixels.                                                                                        |
| mesmer_remove_border_cells | Remove cells that touch the image border.                                                                                               |
| mesmer_pixel_expansion     | Manual pixel expansion after segmentation.                                                                                              |
| mesmer_padding             | Number of pixels to crop the image by on each side before segmentation.                                                                 |

> [!WARNING]
> You cannot run both Mesmer and Cellpose segmentation on the same sample (with
> the same name). If you want to run both on a sample, put it on a different
> line and give it a different sample name.

### Cellpose segmentation

If you want to run Cellpose as your segmentation algorithm, you can specify a
config file like so:

```csv
sample,run_backsub,run_cellpose,tiff,nuclear_channel,membrane_channels
sample1,true,true,/path/to/sample1.tiff,DAPI,CD45:CD8
sample2,false,true,/path/to/sample2.tiff,DAPI,CD45
```

As with Mesmer, nuclear channels only support one entry; membrane channels may
have multiple values separated by `:` characters. You can also set the following
parameters, either via CLI (e.g., `--mesmer_combine_method prod` or in a config
file pass to the workflow via `-c`:

| Parameter Name              | Description                                          |
| --------------------------- | ---------------------------------------------------- |
| cellpose_diameter           | Diameter of cells in pixels for cellpose.            |
| cellpose_min_area           | Minimum area of cells in square pixels for cellpose. |
| cellpose_flow_threshold     | Flow threshold for cellpose.                         |
| cellpose_cellprob_threshold | Cell probability threshold for cellpose.             |
| cellpose_model_type         | Cellpose model type to use for segmentation.         |
| cellpose_pretrained_model   | Path to a pre-trained Cellpose model.                |

Cellpose will run in a parallelised patched workflow using sopa. To control the
patching process, you can use the `patch_width_pixel` and `patch_overlap_pixel`
parameters.

If you want to skip measurements (this may take some time for large images), you
can use set the parameter `skip_measurements` to `true`.

## Dealing with large images

You can run the pipeline with different profiles for different size images:

- `small`: for images <150GB
- `medium`: for images <300GB
- `large`: for images <600GB

## Pipeline output

The pipeline will create the following outputs with background subtraction
(COMET only):

- `extractmarkers/sample.csv` -- marker names, background per channel and
  exposure time
- `backsub/sample.ome.tif` -- bakground subtracted tiff image.

For sopa cellpose segmentation:

- `sopa/sample.zarr` -- SpatialData converted image containing segmentations.
- `parquettotiffwholecell/sample_whole-cell.tiff` -- tiff label masks from cellpose segmentation for
  whole cell segmentation.
- `parquettotiffnuclear/sample_nuclear.tiff` -- tiff label masks from cellpose segmentation for
  nuclear segmentation.

For mesmer segmentation:

- `mesmerwc/sample_whole-cell.tiff` -- tiff label masks for cellpose whole-cell segmentation.
- `mesmernuc/sample_nuclear.tiff` -- tiff label masks for cellpose nuclear segmentation.

And for either method:

- `cellmeasurement/sample.geojson` -- resolved whole-cell and nuclear
  segmentations, optionally containing measurements and intensity values per
  cell, compatible with QuPath.

## Credits

WEHI-SODA-Hub/spatialproteomics was originally written by the WEHI SODA-Hub.

We thank the following people for their extensive assistance in the development of this pipeline:

- Michael McKay (@mikemcka)
- Emma Watson

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use WEHI-SODA-Hub/spatialproteomics for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline was created using the `nf-core` template. You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
