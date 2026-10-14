# Combining membrane channels

When `membrane_channels` lists more than one channel, the pipeline combines
them into one membrane channel before segmentation. This page describes how,
and records three tests of whether rescaling each channel before combining
gives better segmentation.

## What the pipeline does

`combine_method` sets how the channels are combined:

- `max` (default): the per-pixel maximum of the raw intensities.
- `norm_mean`: each channel is rescaled so its 50th percentile is 0 and its
  99.5th percentile is 1 (values below 0 are set to 0), then the per-pixel mean
  is taken and mapped linearly so its maximum is 65535. This is the rescaled
  mean tested below.

`COMBINECHANNELS` does the combination, and Cellpose and CellSAM receive the
combined channel. Mesmer combines the channels inside `mesmer-segment`, which
has only `max`. A sample that runs Mesmer with `norm_mean` and more than one
membrane channel stops the run at startup.

Each segmenter then rescales the channel it receives. Cellpose rescales to the
1st and 99th percentiles; Mesmer and CellSAM also apply contrast enhancement
(CLAHE). This rescaling happens after the combination.

## The concern

Channels have different intensity ranges, set by exposure time, antibody and
dilution. In a raw maximum, the channel with the highest values wins most
pixels. A dim marker can lose everywhere, including in the cells it marks.

Contrast of marker-only pixels in each membrane input, measured on the brain
crop (test 1). For each marker, the pixels where only that marker is bright are
compared with the marker's dim pixels; contrast is the difference in means
divided by the standard deviation of the dim pixels. The value does not change
under a linear rescale, so the segmenters' own percentile rescaling does not
change it.

| Marker      | Pixels | Raw max | Rescaled max | Rescaled mean | Marker alone |
| ----------- | ------ | ------- | ------------ | ------------- | ------------ |
| CD45        | 1378   | 0.50    | 11.59        | 9.06          | 7.31         |
| MHCII       | 439    | 1.42    | 9.75         | 9.22          | 10.67        |
| Iba1        | 490    | 3.60    | 9.70         | 8.87          | 11.40        |
| GFAP        | 4576   | 25.56   | 10.79        | 8.06          | 21.28        |
| Integrin β1 | 6501   | 100.70  | 10.12        | 6.85          | 215.26       |

In the raw maximum, a cell marked only by CD45 is 0.5 background standard
deviations above background, which is within the noise. Rescaling each channel
first gives every marker a contrast of about 7 to 12, but lowers the contrast
of the brightest markers.

## Membrane inputs tested

| Name          | Membrane channel                                                                                 |
| ------------- | ------------------------------------------------------------------------------------------------ |
| Raw max       | per-pixel maximum of the raw channels (the pipeline's method)                                    |
| Rescaled max  | each channel rescaled (median → 0, 99.5th percentile → 1, values below 0 set to 0), then maximum |
| Rescaled mean | the same rescaling, then the per-pixel mean (`combine_method norm_mean`)                         |
| Single        | one marker on its own                                                                            |

Each input was segmented by Cellpose (`cpsam_v2`), Mesmer and CellSAM with the
pipeline defaults (`use_whole_cell_only = true`, `cellpose_min_area = 100`,
`cellsam_min_area = 200`). The final masks from `cellmeasurement/` were scored:

- Cells: number of segmented cells.
- Median area: median cell area in µm².
- Nuclei inside a cell: share of nuclei, found by Laplacian-of-Gaussian blob
  detection on DAPI (radius 2–6 µm), whose centre falls inside a segmented
  cell.
- Ground truth inside a cell (test 2 only): share of manually labelled cells
  whose centroid falls inside a segmented cell.
- Ground truth alone (test 2 only): share of manually labelled cells whose
  centroid falls inside a segmented cell that contains no other labelled
  centroid. Two labelled cells merged into one segmented cell both count as
  failures.

## Test 1: brain

### Data

- One 1024 × 1024 pixel crop of a background-subtracted COMET image of brain
  tissue, 0.28 µm/pixel, 42 channels.
- Nuclear channel: `DAPI`.
- Membrane channels: `VCW_CD45`, `VCW_MHCII`, `VCW_Iba1_AIF1`,
  `VCW_GFAP_diluted`, `VCW_Integrin_b1`.
- Single marker: `VCW_Integrin_b1`.
- 153 nuclei. No ground truth.

### Results

| Segmenter | Input         | Cells | Median area (µm²) | Nuclei inside a cell |
| --------- | ------------- | ----- | ----------------- | -------------------- |
| Cellpose  | Raw max       | 92    | 61.0              | 60.8%                |
| Cellpose  | Rescaled max  | 81    | 56.1              | 55.6%                |
| Cellpose  | Rescaled mean | 77    | 56.8              | 49.7%                |
| Cellpose  | Integrin β1   | 89    | 64.7              | 58.8%                |
| Mesmer    | Raw max       | 104   | 141.4             | 66.0%                |
| Mesmer    | Rescaled max  | 125   | 130.3             | 76.5%                |
| Mesmer    | Rescaled mean | 114   | 142.2             | 75.2%                |
| Mesmer    | Integrin β1   | 94    | 140.9             | 66.0%                |
| CellSAM   | Raw max       | 92    | 107.5             | 64.1%                |
| CellSAM   | Rescaled max  | 75    | 143.3             | 63.4%                |
| CellSAM   | Rescaled mean | 72    | 155.2             | 64.1%                |
| CellSAM   | Integrin β1   | 95    | 84.8              | 66.7%                |

![Cell outlines (white) for each membrane input and segmenter on the brain crop, on CD45 (red), GFAP (green) and DAPI (blue). The first column shows the membrane input.](../images/membrane-combine/overlays.jpg)

### Observations

- The raw maximum is close to Integrin β1 alone. Integrin β1 has the highest
  values, so it wins most pixels, and the segmentations from the two inputs are
  similar.
- Rescaling raises the background noise and the autofluorescent granules,
  which appear in every channel in this tissue. Mesmer and CellSAM draw larger
  outlines around these granules (median area up to 155 µm²).
- Rescaling helped Mesmer (125 cells and 76.5% of nuclei inside a cell, against
  104 and 66.0%), and made Cellpose and CellSAM find fewer cells.
- Every input left a third or more of the nuclei outside a cell. None of the
  five markers outlines most cells, and that limits all four inputs more than
  the combination method does.

## Test 2: atherosclerosis

### Data

- Five 2048 × 2048 pixel crops (573 × 573 µm) from two background-subtracted
  COMET images of human atherosclerotic plaque, 0.28 µm/pixel:

  | Crop      | Image                     | Top-left (y, x) | Single marker   | Nuclei | Labelled cells |
  | --------- | ------------------------- | --------------- | --------------- | ------ | -------------- |
  | b6_immune | batch 6, 25_1592_1A, AS_1 | 18912, 12668    | `R1_CSL_S100A9` | 1796   | 146            |
  | b6_macro  | batch 6, 25_1592_1A, AS_1 | 14816, 12156    | `R1_CSL_S100A9` | 594    | 140            |
  | b2_endo   | batch 2, 24_2721_4A, AS_1 | 3072, 28272     | `R1_CSL_aSMA`   | 907    | 34             |
  | b2_vsmc   | batch 2, 24_2721_4A, AS_1 | 2560, 12656     | `R1_CSL_FAP`    | 497    | 44             |
  | b2_tcell  | batch 2, 24_2721_4A, AS_1 | 9728, 35440     | `R1_CSL_S100A9` | 1243   | 120            |

- Nuclear channel: `R1_DAPI`.
- Membrane channels: the 12 channels used for the existing segmentation of
  this cohort: `R2_CSL_CD11c`, `R2_CSL_VECadherin`, `R1_CSL_aSMA`,
  `R2_CSL_Adipophilin`, `R1_CSL_S100A9`, `R1_CSL_FAP`, `R1_CSL_MMP9`,
  `R1_CSL_CD45`, `R1_CSL_MPO`, `R1_CSL_HistoneH3`, `R1_CSL_CD3`, `R2_CSL_CD68`.
- Single marker: the marker that has the highest value in the most pixels of
  the raw maximum in that crop.
- Ground truth: 484 manually labelled cells (`ground_truth_8_10_25.csv`) inside
  the crops: 143 Other, 135 neutrophils, 45 fibroblasts, 44 VSMCs, 42
  monocytes, 34 endothelial cells, 21 macrophages, 20 T cells. The labels are
  cell centroids, not outlines.
- The existing segmentation of this cohort used Cellpose on the raw maximum of
  the same 12 channels. Its median cell area is 37.5 µm².

### Contrast

Median over the five crops of the contrast defined in
[The concern](#the-concern). `R1_CSL_MMP9` and `R1_CSL_CD45` had fewer than 20
marker-only pixels in b2_vsmc, so their median is over four crops.

| Marker               | Raw max | Rescaled max | Rescaled mean | Marker alone |
| -------------------- | ------- | ------------ | ------------- | ------------ |
| `R2_CSL_CD68`        | 1.4     | 15.4         | 6.0           | 29.2         |
| `R2_CSL_CD11c`       | 1.8     | 10.4         | 5.5           | 25.1         |
| `R2_CSL_Adipophilin` | 2.7     | 10.4         | 5.9           | 30.7         |
| `R1_CSL_CD3`         | 3.8     | 15.9         | 7.3           | 32.0         |
| `R1_CSL_CD45`        | 5.1     | 16.7         | 9.6           | 34.0         |
| `R2_CSL_VECadherin`  | 9.8     | 14.5         | 6.6           | 127.5        |
| `R1_CSL_MMP9`        | 10.1    | 17.4         | 12.6          | 34.8         |
| `R1_CSL_MPO`         | 10.6    | 16.3         | 13.6          | 21.1         |
| `R1_CSL_aSMA`        | 17.4    | 14.5         | 11.0          | 55.4         |
| `R1_CSL_S100A9`      | 19.1    | 14.2         | 11.2          | 28.6         |
| `R1_CSL_HistoneH3`   | 19.2    | 18.6         | 15.8          | 75.0         |
| `R1_CSL_FAP`         | 20.3    | 13.2         | 8.4           | 36.3         |

### Results

Pooled over the five crops: 5037 nuclei and 484 labelled cells. The median
area is the range of the five per-crop medians.

| Segmenter | Input         | Cells | Median area (µm²) | Nuclei inside a cell | Ground truth inside a cell | Ground truth alone |
| --------- | ------------- | ----- | ----------------- | -------------------- | -------------------------- | ------------------ |
| Cellpose  | Raw max       | 3834  | 25.0–76.4         | 77.2%                | 95.7%                      | 95.7%              |
| Cellpose  | Rescaled max  | 3597  | 24.5–72.4         | 74.4%                | 91.9%                      | 91.5%              |
| Cellpose  | Rescaled mean | 3480  | 25.6–75.3         | 72.7%                | 87.6%                      | 87.2%              |
| Cellpose  | Single        | 3689  | 27.0–65.9         | 77.0%                | 94.8%                      | 93.2%              |
| Mesmer    | Raw max       | 4089  | 65.9–139.9        | 90.2%                | 96.3%                      | 90.5%              |
| Mesmer    | Rescaled max  | 4141  | 70.9–128.6        | 91.4%                | 95.5%                      | 91.3%              |
| Mesmer    | Rescaled mean | 3940  | 69.0–129.1        | 88.7%                | 94.8%                      | 89.9%              |
| Mesmer    | Single        | 4069  | 70.1–124.9        | 90.1%                | 96.1%                      | 91.1%              |
| CellSAM   | Raw max       | 3454  | 53.2–115.5        | 80.1%                | 91.3%                      | 86.0%              |
| CellSAM   | Rescaled max  | 3597  | 53.5–108.6        | 82.7%                | 93.4%                      | 88.8%              |
| CellSAM   | Rescaled mean | 3319  | 57.2–121.9        | 78.0%                | 88.0%                      | 83.5%              |
| CellSAM   | Single        | 3547  | 57.3–110.4        | 81.8%                | 91.3%                      | 85.5%              |

Ground truth alone, by cell type (number of labelled cells in brackets):

| Segmenter | Input         | Endothelial (34) | Fibroblasts (45) | Macrophages (21) | Monocytes (42) | Neutrophils (135) | Other (143) | T cells (20) | VSMCs (44) |
| --------- | ------------- | ---------------- | ---------------- | ---------------- | -------------- | ----------------- | ----------- | ------------ | ---------- |
| Cellpose  | Raw max       | 100.0%           | 100.0%           | 100.0%           | 92.9%          | 88.9%             | 98.6%       | 95.0%        | 100.0%     |
| Cellpose  | Rescaled max  | 100.0%           | 97.8%            | 100.0%           | 81.0%          | 80.7%             | 95.8%       | 100.0%       | 100.0%     |
| Cellpose  | Rescaled mean | 100.0%           | 97.8%            | 95.2%            | 61.9%          | 74.8%             | 93.7%       | 95.0%        | 100.0%     |
| Cellpose  | Single        | 100.0%           | 97.8%            | 100.0%           | 88.1%          | 80.7%             | 99.3%       | 100.0%       | 100.0%     |
| Mesmer    | Raw max       | 97.1%            | 91.1%            | 81.0%            | 78.6%          | 88.1%             | 93.0%       | 100.0%       | 95.5%      |
| Mesmer    | Rescaled max  | 97.1%            | 91.1%            | 81.0%            | 81.0%          | 88.9%             | 95.1%       | 95.0%        | 95.5%      |
| Mesmer    | Rescaled mean | 97.1%            | 91.1%            | 81.0%            | 78.6%          | 87.4%             | 92.3%       | 95.0%        | 95.5%      |
| Mesmer    | Single        | 97.1%            | 91.1%            | 81.0%            | 76.2%          | 88.9%             | 95.8%       | 95.0%        | 95.5%      |
| CellSAM   | Raw max       | 97.1%            | 95.6%            | 85.7%            | 66.7%          | 80.7%             | 90.9%       | 95.0%        | 81.8%      |
| CellSAM   | Rescaled max  | 97.1%            | 95.6%            | 85.7%            | 69.0%          | 84.4%             | 93.7%       | 95.0%        | 90.9%      |
| CellSAM   | Rescaled mean | 94.1%            | 93.3%            | 85.7%            | 61.9%          | 72.6%             | 89.5%       | 95.0%        | 93.2%      |
| CellSAM   | Single        | 97.1%            | 95.6%            | 85.7%            | 64.3%          | 75.6%             | 92.3%       | 90.0%        | 93.2%      |

![Cell outlines (white) for each membrane input and segmenter on a 400 × 400 pixel window of the b6_immune crop, on the rescaled maximum of the 12 membrane channels (red) and DAPI (blue). Yellow crosses are labelled cells. The first column shows the membrane input.](../images/membrane-combine/overlays_csl.jpg)

### Observations

- In the raw maximum, `R2_CSL_CD68`, `R2_CSL_CD11c`, `R2_CSL_Adipophilin`,
  `R1_CSL_CD3` and `R1_CSL_CD45` have a contrast of 5.1 or less. Rescaling
  first raises them to 10.4–16.7 and lowers the four brightest markers
  (`R1_CSL_aSMA`, `R1_CSL_S100A9`, `R1_CSL_HistoneH3`, `R1_CSL_FAP`) from
  17.4–20.3 to 13.2–18.6.
- The higher contrast changed the segmentation little. Pooled ground truth
  alone changed from raw max to rescaled max by −4.2 points for Cellpose, +0.8
  for Mesmer and +2.8 for CellSAM.
- Cellpose lost most with rescaling on monocytes (92.9% to 81.0%) and
  neutrophils (88.9% to 80.7%). Most of these are in b6_immune (39 of 42
  monocytes, 92 of 135 neutrophils).
- CellSAM gained most with rescaling in b2_vsmc: nuclei inside a cell rose from
  60.0% to 70.6%, and VSMCs alone from 81.8% to 90.9%.
- The rescaled mean scored lower than the raw max with all three segmenters.
- With all three segmenters, every pooled score of the single marker was
  within 2.5 points of the raw max.
- Cellpose cells in the three batch 2 crops have a median area of 24.5–32.3 µm²
  with every input, which is nucleus-sized, as in the existing segmentation.
  The membrane input does not change this, so the cause is not the combination
  method.
- Mesmer gave the largest cells (median 65.9–139.9 µm²) and the most nuclei
  inside a cell (88.7–91.4%) with every input.

## Test 3: objects without a nucleus, missed cells and outlines

Test 3 used the Cellpose masks of test 2 and asked what the ground-truth
centroids do not show: objects that are not cells, cells that are missed, and
whether outlines follow the cell. Scores are over the five crops.

Objects without a nucleus. A cell counts as low-DAPI when the 90th percentile
of DAPI inside it is below 0.2 times the median of the production cells in the
same crop:

| Input         | Cells | Low-DAPI cells |
| ------------- | ----- | -------------- |
| Raw max       | 3834  | 147 (3.8%)     |
| Rescaled max  | 3597  | 90 (2.5%)      |
| Rescaled mean | 3480  | 71 (2.0%)      |
| Single        | 3689  | 42 (1.1%)      |

Most low-DAPI raw-max cells are small, bright membrane specks and fragments of
fibres in b6_immune. The single marker had the fewest, but in b6_macro it
outlined a few large, smooth objects in diffuse background, which have faint
DAPI and are not counted here.

Missed cells. In b6_immune, 227 raw-max cells with clear DAPI (90th percentile
at least 0.5 times the median production cell) have no matching cell (IoU 0.3
or higher) in the rescaled mean; at 200 of them the rescaled mean has no cell
at all. They are in the dense immune region, where rescaling lifts several dim
markers into an even haze with no edges between cells. In the other four crops
the rescaled mean lacks 1 to 33 such cells. Ground truth inside a cell in
b6_immune: raw max 87.0%, rescaled max 78.1%, rescaled mean 66.4%.

Outlines. For 118 labelled cells (macrophages, T cells, endothelial cells,
fibroblasts, monocytes and 37 of 135 neutrophils), a reference outline was
built from the cell's own marker (CD68, CD3, VE-cadherin, FAP, CD45, S100A9)
and its nucleus, cut at the midline to neighbouring nuclei and at 12 µm from
the centroid. VSMCs were left out: aSMA is a continuous sheet and gives no cell
edge. Mean IoU with the reference, paired over the same cells:

| Comparison              | Difference in IoU | 95% interval   |
| ----------------------- | ----------------- | -------------- |
| Rescaled mean − raw max | −0.006            | −0.029, +0.019 |
| Rescaled max − raw max  | −0.007            | −0.028, +0.016 |

The membrane input did not change the outlines. With every input, Cellpose
outlines lie inside the reference (98% or more of each outline is on it) and
cover 52% to 59% of it, so they are close to nucleus-sized.

Cellpose settings tested on the same crops, with Cellpose run directly on
each crop (default settings: IoU 0.51, 3849 cells), did not give a better
default: `cellpose_cellprob_threshold` −1 raised the IoU to 0.57; a flow
threshold of 2 with cell probability threshold −2 raised it to 0.62 but added
many objects without a nucleus (6117 cells); sopa's own
preprocessing (Gaussian sigma 1 and CLAHE clip limit 0.2) with its defaults
(flow threshold 2, cell probability threshold −6) gave 14594 cells, 69% of
them low-DAPI.

## Conclusion

Rescaling each channel before combining helped some segmenters and made others
worse:

| Segmenter | Brain: nuclei inside a cell, raw max → rescaled max | Atherosclerosis: ground truth alone, raw max → rescaled max |
| --------- | --------------------------------------------------- | ----------------------------------------------------------- |
| Cellpose  | 60.8% → 55.6%                                       | 95.7% → 91.5%                                               |
| Mesmer    | 66.0% → 76.5%                                       | 90.5% → 91.3%                                               |
| CellSAM   | 64.1% → 63.4%                                       | 86.0% → 88.8%                                               |

No input was better than the raw maximum for all three segmenters on both
datasets, so `max` stays the default. For Cellpose, the raw maximum gave the
best ground-truth scores of the four inputs on both datasets.

`norm_mean` is available for Cellpose and CellSAM. In test 3 it gave about half
as many objects without a nucleus as `max` and the same outlines, but missed
nucleated cells where several dim markers overlap in dense tissue. Use it when
`max` outlines non-cell objects that a bright marker creates, and check dense
regions for missed cells.

Limitations:

- Two tissues: one brain crop, and five crops from two atherosclerosis images.
- The test 3 outlines are a marker-based reference, not hand-drawn outlines,
  and cover 118 cells of six types.
- No ground-truth outlines. The atherosclerosis ground truth is cell centroids,
  so it measures whether cells are found and kept apart, not whether their
  outlines are correct.
- Per-marker measures (share of each marker's signal inside cells) were also
  computed for the brain crop, but CD45 was above the positivity threshold on
  about 12% of the image, which is mostly diffuse background in this tissue,
  so they are not reported here. They are in
  `experiments/membrane_combine/brain_per_marker.csv`.
- A panel with one membrane marker per cell type and clear membrane staining
  (for example CD45 with pan-cytokeratin in a tumour) has not been tested.

## Choosing membrane channels

- A membrane marker present on most cells, such as Na⁺/K⁺-ATPase,
  β-catenin or HLA class I, outlines more cells than a mix of cell-type
  markers.
- Each extra channel in a maximum can only raise the background. List the
  channels that outline cells, not every marker in the panel.
- Check the combined image in `combinechannels/<sample>_combined_channels.tiff`
  before running a cohort. If one marker dominates it, the other listed markers
  contribute little.

## Repeating the test

The scripts are in `experiments/membrane_combine/`. Each test is described by
a JSON config: `config_brain_run5.json` (test 1) and `config_csl.json`
(test 2).

1. `make_inputs.py --config <config> --outdir <inputs>` writes the four
   membrane inputs for each region as 2-channel OME-TIFFs and a samplesheet with
   one row per input and segmenter.
2. Run the pipeline on `<inputs>/samplesheet.csv` with `--deepcell_cache_dir`
   set.
3. `analyse.py --config <config> --inputs <inputs> --results <outdir> --out <dir>`
   writes `scores.csv`, `contrast.csv`, `gt_by_type.csv` (if the config has
   ground truth) and one `overlays_<region>.png` per region.

For another image or panel, write a new config. The keys are described at the
top of `make_inputs.py`.

Test 3 scripts and result tables are in `experiments/membrane_combine/test3/`;
see its `README.md`. They read the test 2 inputs and pipeline results from
WEHI paths.
