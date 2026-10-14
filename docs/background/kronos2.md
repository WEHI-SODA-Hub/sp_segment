# KRONOS2 marker handling

This page explains how the pipeline prepares cells for KRONOS2. For how to run
it, see [KRONOS2 embeddings](../how-to/kronos2.md).

## Which cells

The cells come from the cell measurement GeoJSON. Embedding row _i_ is computed
for GeoJSON feature _i_, so embeddings are joined to cells by position, not by
spatial matching.

## Patch content

For each cell, the pipeline cuts a square patch of `kronos_patch_size` pixels
(default 64) centred on the cell's centroid, from every channel KRONOS2 reads.

With `kronos_isolate_cell = true` (the default), pixels outside the cell outline
are set to 0, so the embedding describes that cell and not its neighbours. A
64-pixel patch at 0.28 µm/pixel covers about 18 µm. Set it to `false` to keep
the surrounding pixels.

## Intensity scaling

Pixel values are divided by a number derived from the image data type, as
KRONOS2 does itself:

| Image data type | Divisor |
| --------------- | ------- |
| `uint8`         | 255     |
| `uint16`        | 65535   |
| float           | 400     |

Set `kronos_max_value` to use a different divisor.

## Marker names and normalisation

KRONOS2 normalises each channel with statistics stored for each marker in its
288-name vocabulary. It finds a channel's statistics by matching the channel
name to the vocabulary exactly, ignoring separators. It has no alias or fuzzy
matching.

A channel that does not match gets default statistics. KRONOS2 does not report
this, and the result for that channel is less reliable. The pipeline therefore
checks every name against the vocabulary first, using KRONOS2's own matching
code, and fails the run when a name does not match. You then:

- map the name to a vocabulary name with `kronos_marker_mapping`
- leave the channel out with `kronos_exclude_markers`
- or accept default statistics with `kronos_allow_novel_defaults`

The marker report lists the outcome for every channel.

## DAPI

KRONOS2's only alias is for the nuclear stain: the channel passed as its
`preferred_dapi` gets DAPI statistics whatever its name. The pipeline passes the
samplesheet `nuclear_channel`, or `kronos_nuclear_marker` if set. A DRAQ5 or
Hoechst channel therefore needs no mapping.

## Excluded channels

`kronos_exclude_markers` stops KRONOS2 reading a channel. The pipeline opens only
the remaining channels. KRONOS2 uses all channels it reads together, so the set
of channels affects every cell's embedding. Embeddings made with different
channel sets cannot be compared.

`remove_markers` is different: it removes channels from the background-subtracted
image, so no later step sees them.

## Reproducibility

KRONOS2 runs in 32-bit floating point. The GPU library chooses calculation
methods by batch size, so embeddings change very slightly (about 1e-4) with
`kronos_batch_size`. The default, 16, reproduces the published values.
