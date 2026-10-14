# Testing

## What CI runs

| Workflow                         | Job                       | What it runs                                                                  |
| -------------------------------- | ------------------------- | ----------------------------------------------------------------------------- |
| `.github/workflows/ci.yml`       | `test`                    | `nextflow run . -profile test,singularity` on Nextflow 24.04.2 and the latest release |
| `.github/workflows/linting.yml`  | `pre-commit`              | prettier and editorconfig-checker                                             |
| `.github/workflows/linting.yml`  | `pytest`                  | `pytest tests/python -q`                                                      |
| `.github/workflows/linting.yml`  | `nf-core`                 | `nf-core pipelines lint`                                                      |
| `.github/workflows/docs.yml`     | `build`                   | `generate_reference.py --check` and `mkdocs build --strict`                   |

The `test` profile runs Cellpose with background subtraction on one synthetic
image. CI has no GPU, so Mesmer, CellSAM and KRONOS2 are not run in CI. Run
their test profiles locally on a GPU node before merging a change that affects
them. nf-test is not run in CI (see below).

## Pipeline test profiles

Run every profile affected by a change. See
[Test profiles](../getting-started.md#test-profiles) for the list. From a clone:

```bash
nextflow run . -profile test,singularity --outdir results/test
nextflow run . -profile test_mesmer,singularity --outdir results/test_mesmer
nextflow run . -profile test_cellsam,singularity --outdir results/test_cellsam
KRONOS2_MODEL_DIR=/path/to/KRONOS2 nextflow run . -profile test_mesmer_kronos2,singularity --outdir results/kronos2
```

At WEHI, submit these with `sbatch` and add a resource profile, for example
`-profile test_cellsam,apptainer,small`, so GPU tasks go to `gpuq`.

## Python unit tests

```bash
pip install pytest numpy pandas tifffile rasterio shapely pyarrow typer
pytest tests/python -q
```

`tests/python/conftest.py` adds `bin/` to the import path.

## nf-test

Modules and subworkflows have nf-test tests (`tests/main.nf.test`) with
committed snapshots (`main.nf.test.snap`). Run them with:

```bash
nf-test test --profile test,singularity
nf-test test subworkflows/local/sopa_segment_compartment --profile test,singularity
```

To avoid downloading the 1.2 GB Cellpose weights on every run, export
`CELLPOSE_MODELS_DIR` with a pre-downloaded model directory.
`tests/nextflow.config` passes it to `cellpose_models_dir`.

### Updating snapshots

A change to a process graph changes the `versions` output and invalidates the
snapshot. Update it in two steps:

1. `nf-test test <path> --profile test,singularity --update-snapshot`
2. Run the same test again without `--update-snapshot`. It must pass.

The second run is the check: `--update-snapshot` always agrees with what it
just wrote. Update GPU-dependent snapshots (CellSAM) on a GPU node.

### Known state

- `cellposemodel`, `kronos2embeddings` and `smoothmasks` have no snapshot.
- `sopa_cellpose` has no test.
- `cellsamsegment - nuclear` and `- whole-cell` (module level): mask md5s
  differ from the snapshot. Not diagnosed. Do not regenerate them before the
  CPU determinism question below is answered.
- `segmentationreport - basic functionality`: the HTML md5 differs. Quarto may
  write the render date into the HTML. Not confirmed.
- CellSAM determinism on CPU is not known. `process_gpu` sets
  `accelerator = 1`, so Docker test runs that looked like CPU runs used the GPU.
  To test CPU, set `env.CUDA_VISIBLE_DEVICES = ''` and
  `process { withLabel: process_gpu { accelerator = 0 } }` in
  `tests/nextflow.config` (do not commit this), run `cellsam_segment - tiff`
  twice and compare the outputs.
- The backsub container writes a new UUID into each OME-TIFF, so the
  `backgroundsubtract` test checks the file name, not its md5.

### Before adding nf-test to CI

- GitHub's `ubuntu-latest` runners have no GPU. The `process_gpu` processes are
  `SOPA_SEGMENTATIONCELLPOSE`, `CELLSAMSEGMENT` and `KRONOS2EMBEDDINGS`.
  Decide how they run before enabling tests that reach them.
- `.nf-core.yml` lists `.github/workflows/nf-test.yml` under
  `lint.files_exist`, so lint expects it to be absent. Adding it also turns on
  the `actions_nf_test` lint check. Add it once the suite passes.
- Cache the Cellpose model directory keyed on the sopa container tag, and set
  `CELLPOSE_MODELS_DIR`:

```yaml
- name: Restore Cellpose model cache
  id: cellpose-cache
  uses: actions/cache@v4
  with:
    path: ${{ env.CELLPOSE_MODELS_DIR }}
    key: cellpose-models-${{ env.CELLPOSE_CONTAINER }}

- name: Fetch Cellpose model weights
  if: steps.cellpose-cache.outputs.cache-hit != 'true'
  run: |
    mkdir -p "${CELLPOSE_MODELS_DIR}"
    docker run --rm -u "$(id -u):$(id -g)" \
      -e CELLPOSE_LOCAL_MODELS_PATH=/models \
      -v "${CELLPOSE_MODELS_DIR}:/models" \
      "${CELLPOSE_CONTAINER}" \
      python -c "from cellpose import models; models.CellposeModel(gpu=False)"
```

- CellSAM downloads its 1.7 GB weights in every task unless
  `deepcell_cache_dir` is set. One measured download took 41 minutes.
