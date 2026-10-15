# KRONOS2 embeddings

KRONOS2 is a foundation model for multiplex spatial proteomics. The pipeline
uses it to add a 768-value embedding to every cell in the cell measurement
GeoJSON. Embeddings can be used for clustering and other downstream analysis.

KRONOS2 runs after cell measurement, on the same image the cells were measured
on (the background-subtracted image when `run_backsub = true`). It needs a GPU.

## Get the model

The weights are gated on Hugging Face. Request access to
[MahmoodLab/KRONOS2](https://huggingface.co/MahmoodLab/KRONOS2), then download
the full snapshot:

```bash
hf auth login
hf download MahmoodLab/KRONOS2 --local-dir /path/to/KRONOS2 --exclude "demo_image/*"
```

- Download the whole snapshot. The loader imports the `dinov2` package that is
  part of the snapshot, so a partial download does not work.
- Do not put the snapshot under `/opt`. The host `/opt` would be mounted over
  the container's `/opt`, and the task would fail with
  `/usr/bin/env: 'python3': No such file or directory`.

## Run

```bash
nextflow run WEHI-SODA-Hub/sp_segment \
    -profile apptainer,medium \
    --input samplesheet.csv \
    --outdir results \
    --enable_kronos true \
    --kronos_model_path /path/to/KRONOS2
```

Outputs:

- `cellmeasurement/<sample>.geojson.gz` with `kronos_emb_0` to
  `kronos_emb_767` added to every cell.
- `kronos2embeddings/<sample>_marker_report.txt`, which lists how each channel
  was matched.

See [Outputs](../output.md#kronos2-embeddings).

## Fix unmatched marker names

KRONOS2 has a vocabulary of 288 marker names. It matches a channel name to the
vocabulary exactly, ignoring separators: `CD-8`, `CD_8` and `CD8` match each
other. Spelling variants do not match: `Cytokeritin` does not match
`CYTOKERATIN`.

If any channel name does not match, the run fails and lists the names. Read the
marker report, then do one of the following for each name.

### The channel is a marker under a different name

Map it to the vocabulary name with a JSON file:

```json title="kronos_marker_mapping.json"
{
  "Collagen 4": "COLLAGENiv",
  "Cytokeritin": "CYTOKERATIN",
  "VISA": "VISTA"
}
```

```bash
--kronos_marker_mapping /path/to/kronos_marker_mapping.json
```

The value can also be inline JSON: `--kronos_marker_mapping '{"VISA": "VISTA"}'`.
Keys that start with `_` are ignored, so they can hold comments. A value that
is neither an existing file nor valid JSON fails the run. Mapping keys that
match no channel are listed in the marker report.

The mapping changes only the names given to KRONOS2. Channel names in the image
and in the GeoJSON do not change. An example file is in
[`docs/examples/kronos_marker_mapping.json`](../examples/kronos_marker_mapping.json).

### The channel is not a marker

For autofluorescence, blank or empty channels, stop KRONOS2 from reading the
channel:

```bash
--kronos_exclude_markers 'Autofluorescence'
```

- Separate several names with commas.
- Names are case sensitive and must match the image's channel names, before any
  mapping.
- A name that matches no channel fails the run.
- The nuclear channel cannot be excluded.
- Only KRONOS2 skips the channel. The image, segmentation and the channel's cell
  measurements are not changed.

With background subtraction, KRONOS2 reads the background-subtracted image, so
channels in `remove_markers` are already gone. Do not list them in
`kronos_exclude_markers`, or the run fails with "no such channel".

### Accept default statistics

`--kronos_allow_novel_defaults true` lets the run continue. KRONOS2 then
normalises each unmatched channel with default statistics instead of
marker-specific ones, which makes that channel's contribution less reliable.
Use mapping or exclusion where possible.

## Nuclear channel

KRONOS2 treats one channel as its DAPI channel. The pipeline uses the
samplesheet `nuclear_channel` for this, so slides stained with DRAQ5 or Hoechst
get DAPI statistics. Set `kronos_nuclear_marker` to use a different channel.

## Keep embeddings comparable across a cohort

Use the same settings for every image you want to compare:

- the same `kronos_exclude_markers` list. KRONOS2 uses all channels together,
  so excluding a channel changes every cell's embedding.
- the same `kronos_patch_size` (default 64) and `kronos_isolate_cell` (default
  `true`).
- the same `kronos_batch_size` (default 16). Embeddings change very slightly
  with batch size; 16 reproduces the published values.

Embeddings from KRONOS1 (pipeline versions before 0.5.0, 384 values per cell)
cannot be compared with KRONOS2 embeddings.

For details of normalisation and patch content, see
[KRONOS2 marker handling](../background/kronos2.md).
