# Choose a segmentation method

Choose the method per sample with the `run_cellpose`, `run_mesmer` or
`run_cellsam` column of the samplesheet. Each row can use one method.

## Comparison

|                 | Cellpose                                                   | Mesmer                                                 | CellSAM                                        |
| --------------- | ---------------------------------------------------------- | ------------------------------------------------------ | ---------------------------------------------- |
| Input           | OME-TIFF                                                   | OME-TIFF or MIBI TIFF                                  | OME-TIFF or MIBI TIFF                          |
| Hardware        | GPU                                                        | CPU                                                    | GPU                                            |
| How it runs     | sopa splits the image into patches; one GPU task per patch | one task per mask                                      | one GPU task per mask, tiled inside the task   |
| Default model   | `cpsam_v2` (Cellpose-SAM)                                  | Mesmer (DeepCell)                                      | CellSAM generalist                             |
| Token           | none                                                       | `DEEPCELL_ACCESS_TOKEN` required                       | `DEEPCELL_ACCESS_TOKEN` optional               |
| Size filter     | `cellpose_min_area` (200 pixels²)                          | `mesmer_min_nuclei_area` (40 pixels², nuclei)          | `cellsam_min_area` (200 pixels²)               |
| Edge cells      | kept                                                       | removed (`mesmer_remove_border_cells`)                 | removed (`cellsam_remove_border_cells`)        |
| Main parameters | `cellpose_pretrained_model`, `patch_width_pixel`           | `mesmer_maxima_threshold`, `mesmer_interior_threshold` | `cellsam_bbox_threshold`, `cellsam_block_size` |

All three segment whole cells from the nuclear channel plus the combined
membrane channels. With `use_whole_cell_only = false`, they also segment nuclei
from the nuclear channel alone, and cell measurement matches each nucleus to a
cell.

## Which to use

- Start with Cellpose for COMET and other OME-TIFF images. The defaults need no
  tuning. See [Cellpose settings](../background/cellpose.md).
- Use Mesmer or CellSAM for MIBI images, because Cellpose needs an OME-TIFF.
- Use Mesmer if no GPU is available.
- CellSAM draws larger cell boundaries around membrane signal than Cellpose. On
  one COMET test region, the median CellSAM cell area was about 30% larger than
  the Cellpose `cpsam_v2` median.

Check results in QuPath before choosing for a whole cohort. To compare methods
on the same image, see [Run several methods on one image](multiple-methods.md).

## Run time

- Cellpose runs one short GPU task per patch (`patch_width_pixel`, default
  3000). Many patches can run at once if GPUs are free.
- CellSAM runs one GPU task per mask for the whole image. On whole-slide COMET
  images this has taken from about 1 to 22 hours per image, depending on the
  number of cells.
- Mesmer runs one CPU task per mask.

## Cellpose models

Set `cellpose_pretrained_model`:

| Value         | Model                                                                              |
| ------------- | ---------------------------------------------------------------------------------- |
| `cpsam_v2`    | Default. The current Cellpose-SAM model.                                           |
| `cpsam`       | The previous Cellpose-SAM model. Use it to match results from before cellpose 4.2. |
| `cpdino`      | DINOv3 backbone (ViT-L, about 303 million parameters).                             |
| `cpdino-vitb` | Smaller DINOv3 backbone (ViT-B).                                                   |
| a file path   | Your own trained Cellpose model. The file must exist.                              |

On the 119-cell synthetic test image, compared with the earlier pipeline
version (which used `cpsam`):

| Model      | Cells | Foreground IoU | Cells matched at IoU ≥ 0.90 |
| ---------- | ----- | -------------- | --------------------------- |
| `cpsam`    | 119   | 0.9987         | 119 of 119 (100%)           |
| `cpsam_v2` | 119   | 0.9552         | 112 of 119 (94%)            |
| `cpdino`   | 123   | 0.9251         | 106 of 119 (89%)            |

This shows the models differ. It does not show which is more accurate for your
tissue.
