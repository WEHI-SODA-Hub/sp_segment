# Running at WEHI

This page covers Milton (the WEHI HPC) and the Seqera Platform.

## Milton

### Load Nextflow

```bash
module load nextflow/25.10.0
```

Name the version. The default `nextflow` module is 23.10.0, which is older than
the 24.04.2 the pipeline requires. The module also loads Apptainer and sets:

| Variable                   | Value                                                  |
| -------------------------- | ------------------------------------------------------ |
| `NXF_WORK`                 | `/vast/scratch/users/$USER/nextflow/work`              |
| `NXF_SINGULARITY_CACHEDIR` | `/vast/scratch/users/$USER/nextflow/singularity_cache` |

So the work directory and the converted container images are already on
scratch.

### Move the remaining caches off your home directory

The pipeline downloads tens of GB of container images and about 3 GB of model
weights. Your home directory has a quota. Put every cache on scratch:

```bash
export APPTAINER_CACHEDIR=/vast/scratch/users/$USER/apptainer_cache
export SINGULARITY_CACHEDIR=$APPTAINER_CACHEDIR
mkdir -p "$APPTAINER_CACHEDIR"
```

and pass these parameters:

```bash
--cellpose_models_dir /vast/scratch/users/$USER/cellpose_models \
--deepcell_cache_dir /vast/scratch/users/$USER/deepcell
```

`cellpose_models_dir` must exist and contain the model before the run. See
[Cellpose model weights](usage.md#cellpose-model-weights). If you leave it
unset, the pipeline downloads the weights once per run instead.

`APPTAINER_CACHEDIR` must be exported in the shell that starts Nextflow. It
cannot be set in a config file. See [Caches on a cluster](usage.md#caches-on-a-cluster).

### Submit the run

Do not run Nextflow on a login node. Submit the Nextflow process as a Slurm job.
It submits each pipeline task as its own Slurm job and must keep running until
the pipeline finishes.

```bash title="run_sp_segment.sh"
#!/bin/bash
#SBATCH --job-name=sp_segment
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

module load nextflow/25.10.0
export APPTAINER_CACHEDIR=/vast/scratch/users/$USER/apptainer_cache
export SINGULARITY_CACHEDIR=$APPTAINER_CACHEDIR

nextflow run WEHI-SODA-Hub/sp_segment \
    -r 0.5.0 \
    -profile apptainer,medium \
    --input samplesheet.csv \
    --outdir results \
    --deepcell_cache_dir /vast/scratch/users/$USER/deepcell \
    -resume
```

```bash
sbatch run_sp_segment.sh
```

Set `--time` longer than the whole run is expected to take, and within the
partition's time limit. If the Nextflow job ends early, submit it again: `-resume`
reuses every task that already finished.

The `medium` profile submits tasks to the `regular` partition and GPU tasks to
`gpuq` with `--gres=gpu:1 --constraint=Ampere`. Choose the profile from the
image size and cell density. See [Resource profiles](reference/profiles.md) and
[Large and dense images](how-to/large-images.md).

### Secrets

Set your DeepCell token once. Nextflow stores it in `~/.nextflow/secrets`:

```bash
nextflow secrets set DEEPCELL_ACCESS_TOKEN <your-token>
```

### KRONOS2 weights

Download the KRONOS2 snapshot to a project or scratch directory, not under
`/opt`. See [KRONOS2 embeddings](how-to/kronos2.md#get-the-model).

## Seqera Platform

The pipeline is available on the
[Seqera Platform](https://seqera.services.biocommons.org.au/). The launch form
is built from the pipeline's parameter schema, so each field has the same name
and description as in [Parameters](reference/parameters.md). Parameters marked
hidden (mask smoothing, Cellpose preprocessing) appear only when the form is set
to show hidden parameters.
