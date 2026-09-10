# AVITI: Cellpose 3.x membrane-model support (per-sample-row model choice)

## Problem statement

The AVITI whole-cell/membrane path currently always uses Cellpose v4 SAM
(`cpsam_v2`, `AVITIWHOLECELLSEGMENT`). Elembio also ships several pretrained
Cellpose 3.x membrane models, and users don't always know up front which
model suits their cell type. We need a way to select a Cellpose 3.x membrane
model per sample, reusing the existing Cellpose-3.x script/environment
(`bin/aviti_nuclear_segment.py`, same container as nuclear segmentation)
rather than adding a third environment.

## Decisions confirmed with user

- **One model per samplesheet row, not a list.** Add a single optional
  `membrane_model` column (a model _name_, not a list). Blank/absent means
  Cellpose v4 SAM (`cpsam_v2`, unchanged default behaviour).
- **Multi-model comparison is achieved by adding multiple samplesheet rows**,
  not by fan-out inside the pipeline. Since `sample` must already be unique
  per row, a user wanting to compare models on the same run directory adds
  one row per model, e.g.:

  ```
  sample,run_dir,wells,membrane_model
  run1_jurkat,/path/to/run,A1,jurkat_v1
  run1_pbmc,/path/to/run,A1,pbmc_v1
  run1_default,/path/to/run,A1,
  ```

  All three rows point at the same `run_dir`/well, each gets its own
  independent, unmodified, full pipeline run (tile discovery re-runs per
  row, which is cheap relative to segmentation — no attempt to dedupe tile
  discovery/merge-channels across rows sharing a `run_dir`).

- This means **no in-pipeline fan-out, no channel-plumbing risk, no output
  path/naming changes** for the existing case — a huge simplification over
  the earlier list-based design. Each row behaves exactly like an AVITI
  sample does today; only the whole-cell segmentation module selected for
  that row's tiles differs (branch on `meta.membrane_model` blank vs. set).
- **Composite channel order for the v3 membrane models** (confirmed by
  user): `composite[:, :, 0] = cell`, `composite[:, :, 1] = nuclear`,
  `composite[:, :, 2] = actin` (3-channel mode only). This is a **different
  order** from the existing v4 SAM stack in `aviti_wholecell_segment.py`
  (`nucleus, membrane[, actin]`) — the two scripts build their stacks
  independently and must not share `build_stack()` blindly.
  **Revised during implementation** (matches Elembio's own notebook, and
  confirmed against the Cellpose 3.x `eval()` docstring): rather than
  `channels=[1, 2]`, the actual implementation passes `channels=None,
channel_axis=-1` and instead matches the _number_ of stacked channels
  (2 or 3) to what the selected model's checkpoint was trained with (see
  `model_input_channels()`/`reconcile_membrane_stack()` in
  `bin/aviti_cellpose3_segment.py`). `channel_axis` is **not** v4-only —
  Cellpose 3.x's `eval()` accepts it too (default `None`, meaning
  auto-detected); passing `-1` explicitly is safe and matches how the
  composite is actually stacked (`np.stack(channels, axis=-1)`). With
  `channels=None`, Cellpose does no channel reordering/selection at all, so
  the stack's build order **is** the model's expected channel order —
  correctness now depends entirely on `build_membrane_stack()` matching
  training order, not on a `channels=` mapping.
- **cellpose_min_area / diameter / flow / cellprob thresholds**: reuse the
  same global `cellpose_*` param surface already used for the v4 SAM path —
  no new per-model parameter set unless a specific model needs an override
  later.
- **Model file resolution**: samplesheet carries only the model _name_; a
  new global param `aviti_models_dir` supplies the directory the
  name is resolved against (`${aviti_models_dir}/${name}`),
  mirroring how `aviti_nuclear_model_path` is staged once and reused.
  This is now a shared convention for both nuclear and membrane models:
  `aviti_nuclear_model` and `membrane_model` are both resolved from the
  same `aviti_models_dir`, and the direct-path staging of the nuclear model
  is intentionally retired in favour of consistent name-based lookup.

## Current-state findings (this session's exploration)

- `bin/aviti_nuclear_segment.py`: single-channel input, `channels=[0,0]`,
  Cellpose 3.x, always binarizes to uint8 for the published `Nuclear.tif`,
  optionally also emits an instance-labeled `--output-label` uint16 mask.
  Runs in the `python_numpy_pandas_scikit-image_pruned` container (Cellpose
  3.0.7 + torch 2.10.0+cu128 pinned in
  `modules/local/aviti/nuclearsegment/environment.yml`) — **this is the
  environment the new membrane-model path will reuse.**
- `bin/aviti_wholecell_segment.py`: multi-channel stack builder
  (`nucleus, membrane[, actin]`), Cellpose v4 SAM
  (`models.CellposeModel(pretrained_model=...)`, `channel_axis=-1`, no
  `channels=` argument — v4 SAM doesn't use it), always instance-labeled
  uint16 output, no binarize step. This module/script is untouched by this
  feature; it remains the default path when `membrane_model` is blank.
- `subworkflows/local/aviti_segment/main.nf`: `ch_tiles` is built once per
  (sample, well, tile) via `AVITIDISCOVERTILES.out.manifest.splitCsv(...)
.combine(ch_aviti_samplesheet, by: 0)`. `sample_meta` (from
  `ch_aviti_samplesheet`) is already available at this point, so adding
  `membrane_model` to that meta and copying it onto `meta_tile` is a small,
  local change — no restructuring of the manifest/discovery step needed.
- `assets/schema_input_aviti.json` / samplesheet destructuring in
  `subworkflows/local/utils_nfcore_sp_segment_pipeline/main.nf`: AVITI rows
  are `[sample, run_dir, wells]`, manually destructured into a `meta` map —
  schema property order matters, and adding `membrane_model` here needs the
  same treatment as `wells` (fold into `meta`, default `''`/`null` when
  absent).
- `workflows/sp_segment.nf`: stages `ch_aviti_nuclear_model` once
  (`channel.value(file(...))`) before calling `AVITI_SEGMENT`, with an
  explicit comment about "stage a shared model once" being deliberately kept
  at this level. The membrane-models directory should be staged the same
  way (once, as a value channel), not re-resolved per tile.
- `docs/output.md` / `docs/usage.md`: document the current output layout
  (`avitisegmentation/<sample>/CellSegmentation/Well<well>/`,
  `avitistitched/<sample>/Well<well>/`) — **this layout is unchanged by this
  feature**, since `sample` is already the uniqueness key the user controls.

## Proposed design

### 1) Samplesheet + schema

Add an optional `membrane_model` column to `assets/schema_input_aviti.json`
(single model name, not a list):

```
sample,run_dir,wells,membrane_model
run1_jurkat,/path/to/run,A1,jurkat_v1
run1_default,/path/to/run,A1,
```

- Update the samplesheet destructuring in
  `subworkflows/local/utils_nfcore_sp_segment_pipeline/main.nf` to read the
  4th column and fold it into `meta.membrane_model` (string; `null`/`''`
  when absent, treated as "use v4 SAM"). Optionally, allow the user to
  specify 'cpsam' as the model, which uses the v4 SAM model, just to make
  this extra-clear.
- New global param `aviti_models_dir`: directory containing the
  named Cellpose 3.x model files. This is now required (see parameters section).
- Add an optional `cell_diameter` parameter to the samplesheet params that is
  passed to the segment script. The default value here is 0, which means use
  the model's pretrained cell diameter size.

### 2) Generalise the Cellpose-3.x segmentation script

Rename `bin/aviti_nuclear_segment.py` → `bin/aviti_cellpose3_segment.py`
(update the docstring accordingly), and extend it to cover both current
nuclear behaviour and the new membrane behaviour, selected by a `--mode
{nuclear,membrane}` flag (default `nuclear` to keep existing nuclear
invocation/tests unchanged):

- `nuclear` (current behaviour, unchanged): single positional image input,
  `channels=[0, 0]`, binarize + optional instance-labeled output, filenames
  stay `--output`/`--output-label`.
- `membrane`: takes `--cell-tif`, `--nucleus-tif`, optional `--actin-tif`
  (mirrors `aviti_wholecell_segment.py`'s CLI shape for consistency), builds
  the composite in the **confirmed order** `[cell, nucleus, actin]`,
  calls `model.eval(composite, channels=[1, 2], diameter=..., ...)` (no
  `channel_axis=-1` — that's a v4-only argument), always instance-labeled
  uint16 output (no binarize, no `--output-label` needed — mirrors
  `aviti_wholecell_segment.py`'s single `--output`).
- Both modes keep `remove_small_cells()`, `MASK_COMPRESSION`, `log()`, model
  loading (`models.CellposeModel(gpu=use_gpu, pretrained_model=str(model_path))`)
  shared — only channel-stack construction, `channels=` value, and
  binarize/label-output behaviour differ.
- Update `modules/local/aviti/nuclearsegment/main.nf`'s script line
  (`aviti_nuclear_segment.py` → `aviti_cellpose3_segment.py --mode nuclear`).
  Keep the `AVITINUCLEARSEGMENT` process name (it still only does nuclear
  segmentation).

### 3) New module: `AVITIMEMBRANESEGMENT`

New module (`modules/local/aviti/membranesegment/`), modelled closely on
`AVITINUCLEARSEGMENT`'s container/environment (reuses the same
Cellpose-3.x/torch-cu128 environment.yml and Wave container — no new image
to build) but `AVITIWHOLECELLSEGMENT`'s input shape (`cell/nucleus/actin`
tuple) and output contract (`${meta.tile}_Cell.tif`, matching Cytocanvas
naming exactly like the v4 SAM path does):

```
input:
tuple val(meta), path(nucleus_tif), path(membrane_tif), path(actin_tif)
path model_path   // the one model file this row's meta.membrane_model resolves to

output:
tuple val(meta), path("${meta.tile}_Cell.tif"), emit: cell_mask
path "versions.yml", emit: versions

script:
aviti_cellpose3_segment.py --mode membrane \
    --cell-tif ${membrane_tif} --nucleus-tif ${nucleus_tif} ${actin_arg} \
    --output ${meta.tile}_Cell.tif --model-path ${model_path} ${args}
```

(Naming note: the AVITI manifest's existing `membrane_tif` column is the
"cell paint"/membrane channel file — it maps to the v3 script's `--cell-tif`
argument per the confirmed channel order; no manifest/discovery changes
needed here, just argument-name translation in the module.)

### 4) Subworkflow routing: branch per row on `meta.membrane_model`, no fan-out

Because model selection is now a per-sample-row attribute (not a list),
routing is a simple `branch{}`/`mix()` on the existing `ch_tiles`, with no
fan-out, no reduced-key joins, and no change to how nuclear segmentation,
merged-image, stitching, or downstream (CELLMEASUREMENT/report/KRONOS) work:

- In the `ch_tiles` construction (`AVITIDISCOVERTILES.out.manifest.splitCsv
(...).combine(ch_aviti_samplesheet, by: 0).map { ... }`), copy
  `sample_meta.membrane_model` onto `meta_tile.membrane_model` (same pattern
  already used for `channel_mode`).
- `ch_tiles.branch { v3: it[0].membrane_model; v4: true }` (or equivalent) to
  split into two lanes based on whether `membrane_model` is set.
- v4 lane → `AVITIWHOLECELLSEGMENT(ch_tiles_v4)` (unchanged, as today).
- v3 lane → `AVITIMEMBRANESEGMENT(ch_tiles_v3, resolved_model_path)`, where
  `resolved_model_path` maps `meta.membrane_model` to
  `file("${params.aviti_models_dir}/${meta.membrane_model}")` — if
  multiple distinct model names appear across different rows in one run,
  each `AVITIMEMBRANESEGMENT` task still only stages the one file its row
  needs (Nextflow handles this per-task regardless of other rows' models).
- `cell_mask = AVITIWHOLECELLSEGMENT.out.cell_mask.mix(AVITIMEMBRANESEGMENT.out.cell_mask)`
  — a plain `mix()`, not a `join()`: both branches originated from the same
  `ch_tiles`/`meta_tile`, so no key needs reconciling.
- Everything below this point (`join` with nuclear + image by `meta_tile`,
  `groupTuple` by well, stitching, CELLMEASUREMENT, KRONOS, report) is
  **completely unchanged** — the join still works because every `meta_tile`
  that reaches it, from either branch, is the identical map object shape
  already used by the nuclear/image channels today (`membrane_model` is just
  one more static field mirrored on both sides via the same `ch_tiles`
  source, not something that needs to match on the nuclear/image side).

### 5) Parameters

- `aviti_models_dir` now becomes the shared model directory for both nuclear
  and membrane models (user will most likely keep them in the same
  directory). `aviti_nuclear_model_path` now becomes `aviti_nuclear_model`
  and specifies the nuclear model name only (not the full path). Likewise,
  `membrane_model` is a name resolved under `aviti_models_dir`.
  Make the default `aviti_nuclear_model = 20250212_cellpose_nuc_8diam`.
- Validate in `PIPELINE_INITIALISATION`: for any AVITI samplesheet row using
  either a non-empty `membrane_model` or a nuclear model name, `aviti_models_dir`
  must be set and exist, and the named model file
  (`${aviti_models_dir}/${model_name}`) must exist for every referenced model
  — fail fast, matching the existing nuclear-model-path error style.

### 6) Output layout / docs impact

**No output path changes.** `avitisegmentation/<sample>/CellSegmentation/
Well<well>/` and `avitistitched/<sample>/Well<well>/` are unchanged —
`<sample>` is already whatever unique row name the user chose (e.g.
`run1_jurkat` vs `run1_pbmc`), so comparing models across rows just means
looking in two different `<sample>` directories, no new naming convention
needed.

- Update `docs/usage.md`: document `membrane_model` (per-row, singular),
  `aviti_models_dir`, and the "add one row per model, give each a
  unique sample name" pattern for comparing models on the same run
  directory. Also update the aviti args that we have chaned in this feature
  implementation.
- Update `docs/output.md`: brief note that whole-cell segmentation may come
  from either the v4 SAM path or a named Cellpose 3.x model, selected by the
  samplesheet's `membrane_model` column, with no change to the output
  contract either way.

### 7) Tests

- `tests/python/test_aviti_cellpose3_segment.py` (renamed from
  `test_aviti_nuclear_segment.py`): keep all existing nuclear-mode tests
  (unchanged behaviour with `--mode nuclear` default), add membrane-mode
  tests covering composite channel order `[cell, nucleus, actin]`,
  `channels=[1, 2]`, no-binarize instance-labeled output, and 2-channel
  (no actin) composite.
- `modules/local/aviti/membranesegment/tests/main.nf.test`: stub-mode
  nf-test mirroring `AVITINUCLEARSEGMENT`'s and `AVITIWHOLECELLSEGMENT`'s
  existing test structure.
- `subworkflows/local/aviti_segment` nf-test: extend with a case where
  `membrane_model` is set (routes through `AVITIMEMBRANESEGMENT`) alongside
  the existing blank case (routes through `AVITIWHOLECELLSEGMENT`),
  confirming the `branch`/`mix` doesn't disturb the nuclear/image joins or
  well-level grouping.
- `assets/schema_input_aviti.json` validation tests (if any exist) extended
  for the new optional column.

## Key risks

- **`branch{}` predicate correctness**: must treat blank string, `null`, and
  missing key consistently as "use v4 SAM" — get this wrong and a blank
  `membrane_model` could silently route to the v3 branch with no model
  file. Cover with an explicit nf-test case using an empty string.
- **Channel-order mismatch risk between the two whole-cell scripts**:
  `aviti_wholecell_segment.py` (v4 SAM) stacks `[nucleus, membrane, actin]`;
  the new v3 membrane script stacks `[cell, nucleus, actin]`. These must
  never be unified into one shared `build_stack()` without also unifying
  channel semantics — keep them as two clearly-documented, independent
  functions to avoid an accidental silent swap.
- **Repeated tile discovery per row sharing a `run_dir`**: comparing N
  models on the same run directory now means AVITIDISCOVERTILES/
  AVITIMERGETILECHANNELS/AVITINUCLEARSEGMENT all re-run N times (once per
  sample row), not just whole-cell segmentation. This is accepted as a
  reasonable, simple tradeoff per the user's own proposed samplesheet
  pattern, but worth calling out in docs so users understand the cost of
  comparing many models (N full tile-discovery + nuclear-segmentation
  passes, not just N whole-cell segmentation passes). For this reason,
  specify in the docs that the recommendation is that the user should
  only run this kind of many-model test run on a _subset_ of wells (not
  in the whole sample with all wells).

## Explicitly out of scope for this iteration

- Deduplicating tile discovery / nuclear segmentation / merged-image
  generation across multiple sample rows that share the same `run_dir`.
- Any automatic "best model" selection/scoring — comparing rows/models is a
  manual/downstream activity.
- Changing the existing Cellpose v4 SAM path's behaviour, script, or
  environment.
