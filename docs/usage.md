# Running the pipeline

## Basic command

```bash
nextflow run WEHI-SODA-Hub/sp_segment \
    -r 0.5.0 \
    -profile <container>[,<resources>] \
    --input samplesheet.csv \
    --outdir <OUTDIR>
```

- `--input`: the samplesheet. See [Samplesheet](reference/samplesheet.md).
- `--outdir`: the results directory. See [Outputs](output.md).
- `-profile`: see [Profiles](#profiles).
- `-r`: the release to run. Always set it for real data.

Nextflow creates these in the launch directory:

```text
work/        # task working directories (unless NXF_WORK or -work-dir is set)
<OUTDIR>/    # results
.nextflow/   # cache used by -resume
.nextflow.log
```

## Samplesheet

The samplesheet is a CSV or YAML file with one row per image. The minimum is a
sample name, an image path and one step to run:

```csv
sample,run_backsub,run_cellpose,tiff,nuclear_channel,membrane_channels
slide1,true,true,/data/slide1.ome.tiff,DAPI,CD45:CD8
slide2,false,true,/data/slide2.ome.tiff,DAPI,CD45
```

Rules:

- `sample` must be unique.
- At most one of `run_mesmer`, `run_cellpose` and `run_cellsam` can be `true`
  in a row. To run two methods on one image, see
  [Run several methods on one image](how-to/multiple-methods.md).
- A row must have at least one `run_*` column set to `true`. A row with only
  `run_backsub` set runs background subtraction and nothing else.
- A row that runs a segmentation method needs `nuclear_channel` and
  `membrane_channels`.
- Choose membrane channels that outline most cells. See
  [Combining membrane channels](background/membrane-channels.md).
- Channel names are written without quotes. Separate membrane channels with
  `:`. Channel names can contain spaces but not `"`, `'`, `:` or `,`.
- Paths must not contain spaces.

The full column list is in [Samplesheet](reference/samplesheet.md). An example
is in
[`assets/samplesheet.csv`](https://github.com/WEHI-SODA-Hub/sp_segment/blob/main/assets/samplesheet.csv).

## Parameters

Set parameters on the command line (`--cellpose_diameter 40`) or in a params
file:

```yaml title="params.yaml"
input: samplesheet.csv
outdir: results
cellpose_diameter: 40
generate_report: true
```

```bash
nextflow run WEHI-SODA-Hub/sp_segment -r 0.5.0 -profile singularity -params-file params.yaml
```

> [!WARNING]
> Do not set parameters in a config file passed with `-c`. Use `-c` only for
> resources, executor settings and module arguments. See the
> [nf-core configuration guide](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

Every parameter is listed in [Parameters](reference/parameters.md). Run
`nextflow run WEHI-SODA-Hub/sp_segment --help` to list them in the terminal, and
add `--show_hidden` to include hidden parameters.

## Profiles

Combine profiles with commas, for example `-profile apptainer,medium`. Later
profiles override earlier ones.

- Container profiles: `apptainer`, `singularity`, `docker`, `podman`,
  `shifter`, `charliecloud`. Use one of these. `wave` adds Seqera Wave
  containers and is used together with one of them. `conda` and `mamba` exist
  but are not supported.
- Resource profiles for Slurm: `small`, `medium`, `large`, `wehi_small`,
  `wehi_med`, `wehi_large`. See [Resource profiles](reference/profiles.md).
- Test profiles: see [Test profiles](getting-started.md#test-profiles).
- `debug`: prints each task's host name and keeps work directories.

Without a container profile, every tool must be installed on your `PATH`. This
is not supported.

The pipeline also loads institutional profiles from
[nf-core/configs](https://github.com/nf-core/configs) at startup. Set
`NXF_OFFLINE=true` to skip this when there is no internet access.

## Resuming a run

Add `-resume` to reuse the results of tasks that already finished with the same
inputs and parameters:

```bash
nextflow run WEHI-SODA-Hub/sp_segment -r 0.5.0 -profile singularity --input samplesheet.csv --outdir results -resume
```

To resume a specific earlier run, give its name: `-resume <run-name>`.
`nextflow log` lists run names.

`-resume` restarts a failed task from the beginning. For a cell measurement task
on a very large image, see
[Resume a killed measurement task](how-to/resume-measurement.md).

## Caches on a cluster

The pipeline downloads data that is neither input nor output:

- container images: tens of GB in total
- Cellpose weights: about 1.2 GB
- CellSAM and DeepCell weights: about 1.7 GB

By default all of this goes to your home directory. On a cluster with a home
quota, a full home directory makes image pulls fail with:

```text
FATAL: ... write /home/.../.apptainer/cache/oci-tmp/...: disk quota exceeded
```

Put every cache on scratch:

```bash
SCRATCH=/path/to/scratch/$USER

export APPTAINER_CACHEDIR=$SCRATCH/apptainer_cache
export SINGULARITY_CACHEDIR=$APPTAINER_CACHEDIR
mkdir -p "$APPTAINER_CACHEDIR"

nextflow run WEHI-SODA-Hub/sp_segment \
    -profile apptainer,medium \
    -work-dir $SCRATCH/nextflow/work \
    --container_cache_dir $SCRATCH/apptainer_images \
    --cellpose_models_dir $SCRATCH/cellpose_models \
    --deepcell_cache_dir $SCRATCH/deepcell \
    --input samplesheet.csv \
    --outdir results
```

| Cache                          | Set with                     | Holds                                         |
| ------------------------------ | ---------------------------- | --------------------------------------------- |
| Apptainer OCI layers (largest) | `APPTAINER_CACHEDIR` env var | Image layers unpacked while converting to SIF |
| Converted container images     | `--container_cache_dir`      | The `.img` files Nextflow reuses across runs  |
| Cellpose weights               | `--cellpose_models_dir`      | Pre-downloaded Cellpose models                |
| DeepCell and CellSAM weights   | `--deepcell_cache_dir`       | `<dir>/<username>/.deepcell/models/`          |

Notes:

- `APPTAINER_CACHEDIR` must be exported in the shell that starts Nextflow.
  Nextflow pulls images itself, not inside a task, so a config file cannot set
  it.
- If `--container_cache_dir` is set and `APPTAINER_CACHEDIR` is unset or under
  `$HOME`, the pipeline stops at startup with an error that shows the commands
  to run.
- With Apptainer or Singularity and the layer cache under `$HOME`, the pipeline
  prints a warning at startup.
- The pipeline creates `--deepcell_cache_dir` if it does not exist and mounts
  it into the Mesmer and CellSAM containers.

### Cellpose model weights

Without `--cellpose_models_dir`, the pipeline downloads the Cellpose model once
per run and gives the same copy to every patch task.

If compute nodes have no internet access, download the model once on a login
node, in a Python environment with cellpose 4.2.1.1 installed:

```bash
CELLPOSE_LOCAL_MODELS_PATH=/path/to/scratch/cellpose_models \
    python -c "from cellpose import models; models.CellposeModel(gpu=False, pretrained_model='cpsam_v2')"
```

Then run with `--cellpose_models_dir /path/to/scratch/cellpose_models`.

- Name the model you will run (`cpsam_v2` above).
- The directory must exist. The pipeline checks this at startup.
- The directory must not be under `/opt`. The pipeline rejects this at startup,
  because mounting the host `/opt` hides the container's software.
- cellpose.org serves the current weights for each model name. For results that
  can be repeated exactly, keep using the same `--cellpose_models_dir`.

## Changing resources

Each process has a label (`process_low`, `process_medium`, `process_high`,
`process_multi`, `process_single`, `process_gpu`) that sets its CPUs, memory and
time. The resource profiles set these per label. See
[Resource profiles](reference/profiles.md).

To change one process, write a config file and pass it with `-c`:

```groovy title="custom.config"
process {
    withName: 'CELLMEASUREMENT' {
        cpus   = 32
        memory = 400.GB
        time   = 72.h
    }
}
```

A task that fails with an exit code from 130 to 145, or 104, is retried once
with twice the memory. Other failures stop the pipeline after running tasks
finish.

See the nf-core guides on
[tuning resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources)
and
[tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments).

## Running in the background

The Nextflow process must keep running until the pipeline finishes. Use one of:

- a Slurm job (see [Running at WEHI](wehi.md#submit-the-run))
- `nextflow run ... -bg`, which detaches from the terminal and writes to
  `.nextflow.log`
- `screen` or `tmux`

To limit the memory of the Nextflow process itself, set:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

## Updating the pipeline

Nextflow caches the pipeline code after the first run and keeps using that copy.
To update the cached copy:

```bash
nextflow pull WEHI-SODA-Hub/sp_segment
```

With `-r`, Nextflow uses the named release regardless of the cached default
branch.
