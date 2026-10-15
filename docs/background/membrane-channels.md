# Combining membrane channels

When `membrane_channels` lists more than one channel, the pipeline combines
them into one membrane channel before segmentation. This page describes how,
and what was found when the two methods were compared.

## What the pipeline does

`combine_method` sets how the channels are combined:

- `max` (default): the per-pixel maximum of the raw intensities.
- `norm_mean`: each channel is rescaled so its 50th percentile is 0 and its
  99.5th percentile is 1 (values below 0 are set to 0), then the per-pixel mean
  is taken and mapped linearly so its maximum is 65535.

`COMBINECHANNELS` does the combination, and Cellpose and CellSAM receive the
combined channel. Mesmer combines the channels inside `mesmer-segment`, which
has only `max`. A sample that runs Mesmer with `norm_mean` and more than one
membrane channel stops the run at startup.

Each segmenter then rescales the channel it receives. Cellpose rescales to the
1st and 99th percentiles; Mesmer and CellSAM also apply contrast enhancement
(CLAHE). This rescaling happens after the combination.

## Why there are two methods

Channels have different intensity ranges, set by exposure time, antibody and
dilution. In a raw maximum, the channel with the highest values wins most
pixels. A dim marker can lose everywhere, including in the cells it marks, so
the segmenter may not see the membrane of cells marked only by that marker.

`norm_mean` puts every channel on the same scale before combining, so a dim
marker is not hidden by a bright one. It also lowers the contrast of the
brightest markers, and where several dim markers overlap it can turn their
signal into an even haze with no edges between cells.

## What was found

`max`, `norm_mean`, the maximum of the rescaled channels and a single marker
were compared with Cellpose, Mesmer and CellSAM on multiplexed tissue images,
scored against manually labelled cell centroids and against outlines built
from each labelled cell's own marker and nucleus.

- No method was better than `max` for all three segmenters. With Cellpose,
  `max` found the most labelled cells.
- With Cellpose, `norm_mean` gave about half as many objects without a nucleus
  as `max`. Most of those objects were small, bright membrane specks and
  fragments of fibres.
- With Cellpose, `norm_mean` missed nucleated cells in dense regions where
  several dim markers overlap.
- The method did not change the cell outlines. With every method, Cellpose
  outlines lay inside the marker-based outlines and were close to nucleus-sized.

So `max` stays the default. Use `norm_mean` when `max` outlines non-cell
objects created by a bright marker, and check dense regions for missed cells.

Limitations: the outlines used for scoring were built from marker signal, not
drawn by hand, and the comparison covered a small number of images.

## Choosing membrane channels

- A membrane marker present on most cells, such as Na⁺/K⁺-ATPase,
  β-catenin or HLA class I, outlines more cells than a mix of cell-type
  markers.
- Each extra channel in a maximum can only raise the background. List the
  channels that outline cells, not every marker in the panel.
- Check the combined image in `combinechannels/<sample>_combined_channels.tiff`
  before running a cohort. With `max`, if one marker dominates it, the other
  listed markers contribute little.
