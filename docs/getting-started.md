# Install and test

## Requirements

- Nextflow 24.04.2 or later.
- A container engine: Apptainer, Singularity or Docker. Conda is not supported.
- A GPU for Cellpose, CellSAM and KRONOS2. These steps request one GPU per task
  (the `process_gpu` label). Mesmer, background subtraction and cell
  measurement run on CPU.
- A DeepCell access token for Mesmer. CellSAM uses it if present. See
  [Secrets](#secrets).

If you are new to Nextflow, see the
[nf-core installation guide](https://nf-co.re/docs/usage/installation).

## Get the pipeline

Nextflow downloads the pipeline the first time you run it. To download or
update it without running:

```bash
nextflow pull WEHI-SODA-Hub/sp_segment
```

Use `-r` to run a fixed release, for example `-r 0.5.0`. Releases are listed on
the [releases page](https://github.com/WEHI-SODA-Hub/sp_segment/releases).
Record the release you used so the run can be repeated.

## Secrets

Mesmer downloads its model from DeepCell and needs an access token. Create one
at [users.deepcell.org](https://users.deepcell.org/login/), then store it as a
Nextflow secret:

```bash
nextflow secrets set DEEPCELL_ACCESS_TOKEN <your-token>
```

CellSAM uses the same secret to download the CellSAM v1.2 weights. Without it,
CellSAM uses its default weights.

## Test profiles

Each test profile runs a small synthetic COMET image
(`tests/data/comet/synthetic_multichannel.tif`). Combine a test profile with a
container profile, for example `-profile test,singularity`.

| Profile                   | Method              | Background subtraction | Notes                                                             |
| ------------------------- | ------------------- | ---------------------- | ----------------------------------------------------------------- |
| `test`                    | Cellpose            | yes                    | This is the profile CI runs.                                      |
| `test_no_backsub`         | Cellpose            | no                     |                                                                   |
| `test_multisample`        | Cellpose and Mesmer | Mesmer sample only     | Three samples: two Cellpose, one Mesmer.                          |
| `test_mesmer`             | Mesmer              | yes                    | Needs `DEEPCELL_ACCESS_TOKEN`.                                    |
| `test_mesmer_no_backsub`  | Mesmer              | no                     | Needs `DEEPCELL_ACCESS_TOKEN`. Currently fails: see below.        |
| `test_mesmer_kronos2`     | Mesmer              | yes                    | Also runs KRONOS2. Set `KRONOS2_MODEL_DIR` to a KRONOS2 snapshot. |
| `test_cellsam`            | CellSAM             | yes                    |                                                                   |
| `test_cellsam_no_backsub` | CellSAM             | no                     |                                                                   |

The test samplesheets give image paths relative to the repository root, so run
the test profiles from the root of a clone:

```bash
git clone https://github.com/WEHI-SODA-Hub/sp_segment.git
cd sp_segment
nextflow run . -profile test_cellsam,singularity --outdir results
```

The test profiles limit each task to 4 CPUs, 15 GB of memory and 1 hour.

`test_mesmer_no_backsub` fails with `PhysicalSizeX or PhysicalSizeY not found
in the XML`. Mesmer needs the pixel size in the OME-XML, and the synthetic test
image has none. Background subtraction adds it, so `test_mesmer` passes.

## First run on your data

Write a samplesheet with one row per image. See
[Samplesheet](reference/samplesheet.md) for every column.

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

The cells are in `results/cellmeasurement/slide1.geojson.gz`. See
[View results in QuPath](how-to/qupath.md).

To run on a Slurm cluster, add a resource profile such as `medium`. See
[Resource profiles](reference/profiles.md). For all run options, see
[Running the pipeline](usage.md).
