# Resource profiles

The resource profiles run every task as a Slurm job. Combine one with a
container profile, for example `-profile apptainer,medium`.

All resource profiles:

- submit CPU tasks to the `regular` partition and GPU tasks (`process_gpu`) to
  `gpuq` with `--gres=gpu:1 --constraint=Ampere`
- give GPU tasks the `--nv` (Apptainer, Singularity) or `--gpus all` (Docker)
  container option
- set a time limit of 48 hours per task
- submit at most 10 jobs at a time (`queueSize = 10`), one every 10 seconds
- multiply memory by the attempt number, so a retried task gets twice the memory

The partition names and GPU constraint are those of WEHI's Milton cluster. On
another cluster, write your own profile or override `queue` and
`clusterOptions`.

## Which profile

| Profile      | Use for                                | Target node                            |
| ------------ | -------------------------------------- | -------------------------------------- |
| `small`      | images up to 150 GB                    | Milton sml (28 cores, 113 GB)          |
| `medium`     | images up to 300 GB                    | Milton med/il (28 or 48 cores, 491 GB) |
| `large`      | images up to 600 GB                    | Milton lrg (64 cores, 1498 GB)         |
| `wehi_small` | dense samples on sml nodes             | Milton sml (28 cores, 113 GB)          |
| `wehi_med`   | dense samples on med/il nodes          | Milton med/il (28 or 48 cores, 491 GB) |
| `wehi_large` | very dense samples (millions of cells) | Milton lrg (64 cores, 1498 GB)         |

Memory use depends more on the number of cells than on the file size. The
`wehi_*` profiles request most of a node and suit images with many cells. See
[Large and dense images](../how-to/large-images.md).

## Requests per label

CPUs / memory on the first attempt. Every label has a 48-hour time limit.

| Label            | `small`    | `medium`    | `large`     | `wehi_small` | `wehi_med`  | `wehi_large` |
| ---------------- | ---------- | ----------- | ----------- | ------------ | ----------- | ------------ |
| `process_single` | 1 / 36 GB  | 1 / 120 GB  | 1 / 450 GB  | 1 / 56 GB    | 1 / 244 GB  | 1 / 744 GB   |
| `process_medium` | 1 / 24 GB  | 1 / 48 GB   | 1 / 192 GB  | 1 / 56 GB    | 1 / 244 GB  | 1 / 744 GB   |
| `process_high`   | 1 / 36 GB  | 1 / 120 GB  | 1 / 450 GB  | 1 / 56 GB    | 1 / 244 GB  | 1 / 744 GB   |
| `process_multi`  | 14 / 36 GB | 24 / 120 GB | 48 / 450 GB | 14 / 56 GB   | 24 / 244 GB | 48 / 744 GB  |
| `process_gpu`    | 12 / 48 GB | 16 / 120 GB | 24 / 450 GB | 12 / 244 GB  | 16 / 244 GB | 24 / 244 GB  |

`process_low` is not set by these profiles and uses the default: 1 CPU, 6 GB,
6 hours.

## Default requests

Without a resource profile, tasks run on the local machine with these requests
(from `conf/base.config`):

| Label            | CPUs | Memory | Time |
| ---------------- | ---- | ------ | ---- |
| `process_single` | 1    | 6 GB   | 6 h  |
| `process_low`    | 1    | 6 GB   | 6 h  |
| `process_medium` | 1    | 12 GB  | 48 h |
| `process_high`   | 1    | 24 GB  | 48 h |
| `process_multi`  | 8    | 48 GB  | 48 h |
| `process_gpu`    | 8    | 48 GB  | 48 h |

Memory and the 6-hour limits are multiplied by the attempt number.

## Labels used by each process

| Process                     | Label            |
| --------------------------- | ---------------- |
| `EXTRACTMARKERS`            | `process_low`    |
| `BACKSUB`                   | `process_single` |
| `COMBINECHANNELS`           | `process_high`   |
| `SOPA_CONVERT`              | `process_high`   |
| `SOPA_PATCHIFYIMAGE`        | `process_low`    |
| `SOPA_SEGMENTATIONCELLPOSE` | `process_gpu`    |
| `SOPA_RESOLVECELLPOSE`      | `process_low`    |
| `PARQUETTOTIFF`             | `process_medium` |
| `CELLPOSEMODEL`             | `process_single` |
| `MESMERSEGMENT`             | `process_multi`  |
| `CELLSAMSEGMENT`            | `process_gpu`    |
| `SMOOTHMASKS`               | `process_medium` |
| `CELLMEASUREMENT`           | `process_multi`  |
| `KRONOS2EMBEDDINGS`         | `process_gpu`    |
| `SEGMENTATIONREPORT`        | `process_high`   |

To change one process, see
[Changing resources](../usage.md#changing-resources).
