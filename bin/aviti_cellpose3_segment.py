#!/usr/bin/env python
'''
Module      : aviti_cellpose3_segment
Description : Cellpose 3.x segmentation for AVITI tiles in either nucleus or
              membrane mode. This keeps AVITI's custom nuclear model and the
              Cellpose 3.x membrane models in the same script/environment, while
              preserving the existing nuclear path semantics and adding the
              per-sample membrane model selection path used by AVITI whole-cell
              segmentation.
Copyright   : (c) WEHI SODA Hub, 2026
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import sys
from pathlib import Path
from typing import Annotated, Optional

import numpy as np
import tifffile
import typer

app = typer.Typer(add_completion=False)

MASK_COMPRESSION = "zlib"


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def read_channel(path: Path) -> np.ndarray:
    img = tifffile.imread(path)
    if img.ndim != 2:
        raise ValueError(f"Expected a single-page 2D image, got shape {img.shape} from {path}")
    return img


def model_input_channels(model_path: Path) -> Optional[int]:
    '''
    Read how many input channels the Cellpose checkpoint at ``model_path`` was
    trained with, taken from the first conv weight found in the state dict
    (shape ``[out_ch, in_ch, k, k]``, so ``in_ch`` is the channel count).
    Deliberately does not hardcode a specific parameter name -- Cellpose's
    exact internal layer naming is version-dependent and not part of its
    public API, whereas "the first conv weight registered on the network is
    the input-facing one" holds for any standard CNN state dict, since
    ``state_dict()`` preserves module-registration (i.e. definition) order.
    Used to set ``nchan`` when constructing the model so it matches the
    checkpoint rather than the tile. Returns None if the file cannot be
    inspected, in which case the caller falls back to the input's own
    channel count.
    '''
    import torch

    state_dict = None
    for weights_only in (True, False):
        try:
            state_dict = torch.load(model_path, map_location="cpu", weights_only=weights_only)
            break
        except Exception:
            continue
    if state_dict is None:
        return None
    if isinstance(state_dict, dict) and "model_state_dict" in state_dict:
        state_dict = state_dict["model_state_dict"]
    if not isinstance(state_dict, dict):
        return None

    for value in state_dict.values():
        if getattr(value, "ndim", None) == 4:
            return int(value.shape[1])
    return None


def reconcile_membrane_stack(stack: np.ndarray, model_nchan: int) -> np.ndarray:
    '''
    Make the composite's channel count match what the membrane model was
    trained with. Element ships a '_2ch' and a '_3ch' model per cell type:
    the '_3ch' model expects [membrane, nucleus, actin], the '_2ch' model
    expects [membrane, nucleus]. If the run's cell-paint mode does not match
    the chosen model we adapt rather than crash mid-run, but the right fix is
    to point membrane_model at the model that matches the run.
    '''
    have = stack.shape[-1]
    if model_nchan == have:
        return stack

    log("")
    log("=" * 78)
    log(f"WARNING: membrane model expects {model_nchan} input channel(s) but this "
        f"tile has {have}.")
    if model_nchan == 3 and have == 2:
        log("  This run has no actin channel. Padding a zero actin plane so the")
        log("  '_3ch' model can run -- prefer the matching '_2ch' model for a")
        log("  nucleus + membrane (Cell Paint only) run.")
        stack = np.concatenate([stack, np.zeros_like(stack[..., :1])], axis=-1)
    elif model_nchan == 2 and have == 3:
        log("  Dropping the actin plane so the '_2ch' model can run -- prefer the")
        log("  matching '_3ch' model to make use of the actin channel.")
        stack = stack[..., :2]
    else:
        log("  Cannot reconcile automatically; Cellpose will likely fail to load "
            "this model.")
    log("=" * 78)
    log("")
    return stack


def remove_small_cells(mask: np.ndarray, min_area: int) -> np.ndarray:
    if min_area <= 0:
        return mask
    labels, counts = np.unique(mask, return_counts=True)
    small = labels[(labels != 0) & (counts < min_area)]
    if small.size == 0:
        return mask
    out = mask.copy()
    out[np.isin(out, small)] = 0
    return out


def binarize(mask: np.ndarray) -> np.ndarray:
    return (mask > 0).astype(np.uint8)


def build_membrane_stack(cell_tif: Path, nucleus_tif: Path, actin_tif: Optional[Path]) -> np.ndarray:
    '''
    Stack channels in the confirmed training order [cell, nucleus[, actin]].
    This order is load-bearing: the membrane models are called with
    ``channels=None``, so Cellpose does no channel reordering/selection --
    whatever order this function stacks in is exactly what the model sees.
    This is a *different* order from the v4 SAM stack built independently by
    aviti_wholecell_segment.py ([nucleus, membrane[, actin]]) -- do not
    unify the two without also unifying channel semantics.
    '''
    cell = read_channel(cell_tif)
    nucleus = read_channel(nucleus_tif)
    if cell.shape != nucleus.shape:
        raise ValueError(f"Channel shape mismatch: cell {cell.shape} vs nucleus {nucleus.shape}")

    channels = [cell, nucleus]
    if actin_tif is not None:
        actin = read_channel(actin_tif)
        if actin.shape != cell.shape:
            raise ValueError(f"Channel shape mismatch: cell {cell.shape} vs actin {actin.shape}")
        channels.append(actin)

    return np.stack(channels, axis=-1)


@app.command()
def main(
    mode: Annotated[str, typer.Option(help="Cellpose 3.x segmentation mode: 'nuclear' or 'membrane'.")] = 'nuclear',
    nucleus_tif: Annotated[Optional[Path], typer.Option('--nucleus-tif', exists=True, help="Nucleus channel TIFF for the nuclear mode.")] = None,
    cell_tif: Annotated[Optional[Path], typer.Option('--cell-tif', exists=True, help="Cell-paint / membrane channel TIFF for membrane mode.")] = None,
    actin_tif: Annotated[Optional[Path], typer.Option('--actin-tif', exists=True, help="Optional actin TIFF for 3-channel membrane mode.")] = None,
    output: Annotated[Optional[Path], typer.Option(help="Output TIFF path.")] = None,
    output_label: Annotated[Optional[Path], typer.Option(help="Optional instance-labeled output TIFF for nuclear mode.")] = None,
    model_path: Annotated[Optional[Path], typer.Option(exists=True, help="Path to the Cellpose 3.x model checkpoint.")] = None,
    diameter: Annotated[float, typer.Option(help="Expected cell or nucleus diameter in pixels. 0 uses the model's own trained diameter.")] = 0.0,
    flow_threshold: Annotated[float, typer.Option()] = 0.4,
    cellprob_threshold: Annotated[float, typer.Option()] = 0.0,
    min_area: Annotated[int, typer.Option(help="Discard objects smaller than this many px^2. 0 disables.")] = 0,
    gpu: Annotated[bool, typer.Option(help="Run on GPU.")] = True,
):
    if mode not in {'nuclear', 'membrane'}:
        raise typer.BadParameter("--mode must be one of: nuclear, membrane")
    if output is None:
        raise typer.BadParameter("--output is required")
    if model_path is None:
        raise typer.BadParameter("--model-path is required")

    import torch
    from cellpose import models

    log(f"PyTorch version: {torch.__version__}")
    log(f"CUDA available: {torch.cuda.is_available()}")
    use_gpu = gpu and torch.cuda.is_available()
    if gpu and not use_gpu:
        log("WARNING: --gpu requested but no CUDA device is available; running on CPU.")

    if mode == 'nuclear':
        if nucleus_tif is None:
            raise typer.BadParameter("--nucleus-tif is required in nuclear mode")
        img = read_channel(nucleus_tif)
        # nchan must match the checkpoint, not the input: the AVITI nuclear
        # models (e.g. 20250212_cellpose_nuc_8diam) are 2-channel networks fed
        # a grayscale image via channels=[0, 0], as in Elembio's own notebook.
        model_nchan = model_input_channels(model_path) or 2
        model = models.CellposeModel(gpu=use_gpu, pretrained_model=str(model_path), nchan=model_nchan)
        masks, _flows, _styles = model.eval(
            img,
            channels=[0, 0],
            diameter=diameter if diameter > 0 else None,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
        )
        masks = remove_small_cells(masks, min_area)
        n_objects = len(np.unique(masks)) - 1
        if output_label is not None:
            tifffile.imwrite(output_label, masks.astype(np.uint16), imagej=True, compression=MASK_COMPRESSION)
        binary = binarize(masks)
        tifffile.imwrite(output, binary, imagej=True, compression=MASK_COMPRESSION)
        log(f"Nuclear segmentation: {n_objects} nuclei written to {output}")
        return

    if cell_tif is None or nucleus_tif is None:
        raise typer.BadParameter("--cell-tif and --nucleus-tif are required in membrane mode")

    stack = build_membrane_stack(cell_tif, nucleus_tif, actin_tif)
    log(f"Built {stack.shape[-1]}-channel stack with shape {stack.shape} for membrane segmentation")

    # nchan follows the checkpoint, not the tile. Element ships a matched
    # '_2ch'/'_3ch' pair per cell type; reconcile the composite to whichever
    # was selected and feed it as-is (channels=None), like Elembio's notebook.
    model_nchan = model_input_channels(model_path)
    if model_nchan is not None:
        stack = reconcile_membrane_stack(stack, model_nchan)
    nchan = model_nchan or stack.shape[-1]

    model = models.CellposeModel(gpu=use_gpu, pretrained_model=str(model_path), nchan=nchan)
    masks, _flows, _styles = model.eval(
        stack,
        channels=None,
        channel_axis=-1,
        diameter=diameter if diameter > 0 else None,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )
    masks = remove_small_cells(masks.astype(np.uint16), min_area)
    tifffile.imwrite(output, masks, imagej=True, compression=MASK_COMPRESSION)
    n_cells = len(np.unique(masks)) - 1
    log(f"Membrane segmentation: {n_cells} cell(s) written to {output}")


if __name__ == "__main__":
    app()
