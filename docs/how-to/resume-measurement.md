# Resume a killed measurement task

On whole-slide images with millions of cells, polygon extraction is the
slowest part of `CELLMEASUREMENT`. If the task is killed by its time limit or
by running out of memory, Nextflow's retry and `-resume` start it again from
the beginning.

Set `geometry_checkpoint_dir` to save progress. The task then writes each
completed batch of cells to that directory, and a retried or resumed task
continues from the last completed batch.

## Set it up

```bash
nextflow run WEHI-SODA-Hub/sp_segment \
    -profile apptainer,medium \
    --input samplesheet.csv \
    --outdir results \
    --geometry_checkpoint_dir /vast/scratch/users/$USER/sp_segment_checkpoints \
    -resume
```

- The directory must be outside the Nextflow work directory. Each retry runs in
  a new work directory, so a checkpoint inside the work directory is not found.
- Each sample and each mask gets its own subdirectory.
- `geometry_batch_size` (default 2000) sets the number of cells per batch.
  Smaller batches save progress more often but add overhead. A killed task loses
  at most one batch per worker.

## When checkpoints are not used

A checkpoint is ignored, and extraction starts again, if any of these differ
from the run that wrote it:

- polygon simplification settings
- `geometry_batch_size`
- mask dimensions
- number of cells

## Clean up

Checkpoints are not deleted automatically. Delete the directory after the run
has finished successfully.
