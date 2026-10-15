#!/usr/bin/env python3
import argparse
import fcntl
import sys
import os
import shutil

import numpy as np
import tifffile
from cellSAM.cellsam_pipeline import cellsam_pipeline
import cellSAM

from tiff_channels import get_channel_names, read_channels


# This is an intermediate mask, smoothed downstream, so compress it: label masks
# are long runs of a repeated ID and deflate shrinks them ~90x. Only the final
# smoothed mask is left uncompressed for viewers.
MASK_COMPRESSION = "zlib"
MASK_COMPRESSION_ARGS = {"level": 1}


def download_model_weights():
    '''Download latest CellSAM model weights (v1.2) from users.deepcell.org.'''
    from pathlib import Path

    if 'DEEPCELL_ACCESS_TOKEN' in os.environ:
        print("Downloading/verifying latest CellSAM model weights (v1.2)...")

        # The whole-cell and nuclear tasks share this cache and run at the same
        # time. Without the lock, one task could see the other's download in
        # progress as a missing file and delete it. MESMERSEGMENT uses the
        # same pattern.
        deepcell_dir = Path.home() / ".deepcell"
        deepcell_dir.mkdir(parents=True, exist_ok=True)
        with open(deepcell_dir / "cellsam_download.lock", "a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)

            # Clear corrupted model cache if the .pt file exists but is invalid
            model_dir = deepcell_dir / "models" / "cellsam_v1.2"
            model_file = model_dir / "cellsam_general.pt"

            model_file_corrupted = model_file.exists() and model_file.stat().st_size < 1_000_000
            model_file_missing = model_dir.exists() and not model_file.exists()
            if model_file_corrupted or model_file_missing:
                print(f"Warning: {model_file} appears corrupted or missing. Removing and re-downloading...")
                shutil.rmtree(model_dir)

            cellSAM.get_model()
    else:
        print("Warning: DEEPCELL_ACCESS_TOKEN not set. Using default model weights.")


def find_channel_index(channel_name, channel_names):
    '''
    Find the index of a channel by name.

    First tries exact match, then case-insensitive match.
    Returns None if not found.
    '''
    # Exact match
    if channel_name in channel_names:
        return channel_names.index(channel_name)

    # Case-insensitive match
    for i, name in enumerate(channel_names):
        if name.lower() == channel_name.lower():
            return i

    return None


def create_cellsam_input(img, nuc_idx, mem_idx=None):
    '''
    Create a 3-channel array for CellSAM input.

    Format: [blank, nuclear, membrane]
    If membrane channel is not provided, it remains blank (zeros).
    '''
    if img.ndim == 2:
        h, w = img.shape
        cellsam_img = np.zeros((h, w, 3), dtype=img.dtype)
        cellsam_img[..., 1] = img
    else:
        h, w = img.shape[1], img.shape[2]
        cellsam_img = np.zeros((h, w, 3), dtype=img.dtype)
        cellsam_img[..., 1] = img[nuc_idx, ...]

        if mem_idx is not None:
            cellsam_img[..., 2] = img[mem_idx, ...]

    return cellsam_img


def remove_border_cells(mask):
    '''
    Remove cells that touch the image borders.

    Args:
        mask: 2D numpy array with cell labels

    Returns:
        Filtered mask with border-touching cells removed (set to 0)
    '''
    h, w = mask.shape
    border_cells = set()

    # Find cells touching top and bottom edges
    border_cells.update(mask[0, :])
    border_cells.update(mask[h-1, :])

    # Find cells touching left and right edges
    border_cells.update(mask[:, 0])
    border_cells.update(mask[:, w-1])

    # Remove background (0) from border cells set
    border_cells.discard(0)

    # Create filtered mask
    filtered_mask = mask.copy()
    if border_cells:
        filtered_mask[np.isin(mask, list(border_cells))] = 0

    n_removed = len(border_cells)
    n_total = len(np.unique(mask)) - 1  # -1 for background
    print(f"Removed {n_removed} border-touching cells (kept {n_total - n_removed}/{n_total})")

    return filtered_mask


def remove_small_cells(mask, min_area):
    '''
    Remove cells smaller than a minimum pixel area.

    Args:
        mask: 2D numpy array with cell labels
        min_area: Minimum number of pixels for a cell to be retained

    Returns:
        Filtered mask with small cells removed (set to 0)
    '''
    if min_area <= 0:
        return mask

    labels, counts = np.unique(mask, return_counts=True)
    small = labels[(labels != 0) & (counts < min_area)]

    if small.size == 0:
        print(f"No cells below min-area {min_area}px")
        return mask

    filtered = mask.copy()
    filtered[np.isin(filtered, small)] = 0

    n_total = int((labels != 0).sum())
    print(f"Removed {len(small)} cells below {min_area}px (kept {n_total - len(small)}/{n_total})")
    return filtered


def extract_channels(tiff_path, nuclear_channel, membrane_channels, compartment):
    '''
    Extract and format channels from a multi-channel TIFF for CellSAM.

    Args:
        tiff_path: Path to the input TIFF file
        nuclear_channel: Name of the nuclear channel
        membrane_channels: List of membrane channel names (for whole-cell segmentation)
        compartment: 'nuclear' or 'whole-cell'

    Returns:
        3-channel numpy array formatted for CellSAM [blank, nuclear, membrane]
    '''
    channel_names = get_channel_names(tiff_path)
    print(f"Available channels: {channel_names}")

    nuc_idx = find_channel_index(nuclear_channel, channel_names)
    if nuc_idx is None:
        raise ValueError(f"Nuclear channel '{nuclear_channel}' not found in {channel_names}")

    # Find membrane channel (for whole-cell segmentation)
    mem_idx = None
    if compartment == "whole-cell" and membrane_channels:
        for mem_chan in membrane_channels:
            mem_idx = find_channel_index(mem_chan, channel_names)
            if mem_idx is not None:
                break

        if mem_idx is None:
            others = [i for i in range(len(channel_names)) if i != nuc_idx]
            if others:
                print("Warning: Membrane channel not found. Using mean of non-nuclear channels.")
                nuclear = read_channels(tiff_path, [nuc_idx])[0]
                # Cast back to the image dtype so the stack is not promoted to float64.
                membrane = read_channels(tiff_path, others).mean(axis=0).astype(nuclear.dtype)
                return create_cellsam_input(np.stack([nuclear, membrane]), 0, 1)

    # Only the channels CellSAM uses are decoded.
    indices = [nuc_idx] + ([mem_idx] if mem_idx is not None else [])
    planes = read_channels(tiff_path, indices)
    return create_cellsam_input(planes, 0, 1 if mem_idx is not None else None)


def parse_arguments():
    '''Parse command-line arguments for CellSAM segmentation.'''
    parser = argparse.ArgumentParser(description='CellSAM segmentation')
    parser.add_argument('tiff', help='Input TIFF file path')
    parser.add_argument('--output', required=True, help='Output mask file path')
    parser.add_argument('--compartment', required=True, choices=['nuclear', 'whole-cell'],
                      help='Segmentation compartment')
    parser.add_argument('--nuclear-channel', required=True, help='Nuclear channel name')
    parser.add_argument('--membrane-channel', action='append', default=[],
                      help='Membrane channel name(s)')
    parser.add_argument('--bbox-threshold', type=float, default=0.4,
                      help='Bounding box threshold')
    parser.add_argument('--block-size', type=int, default=712,
                      help='Block size for tiling')
    parser.add_argument('--overlap', type=int, default=120,
                      help='Overlap between tiles')
    parser.add_argument('--iou-depth', type=int, default=120,
                      help='IOU depth parameter')
    parser.add_argument('--iou-threshold', type=float, default=0.5,
                      help='IOU threshold for merging')
    parser.add_argument('--use-wsi', action='store_true',
                      help='Use WSI mode')
    parser.add_argument('--gauge-cell-size', action='store_true',
                      help='Gauge cell size')
    parser.add_argument('--low-contrast-enhancement', action='store_true',
                      help='Apply low contrast enhancement')
    parser.add_argument('--model-path', type=str, default=None,
                      help='Custom model path (optional)')
    parser.add_argument('--remove-border-cells', action=argparse.BooleanOptionalAction, default=True,
                      help='Remove cells touching image borders (default: on)')
    parser.add_argument('--min-area', type=int, default=0,
                      help='Minimum cell area in pixels. Cells smaller than this are removed (default: 0 = disabled)')

    return parser.parse_args()


def run_segmentation(args):
    '''Run CellSAM segmentation pipeline.'''
    import torch
    print(f"PyTorch version: {torch.__version__}", file=sys.stderr)
    print(f"CUDA available: {torch.cuda.is_available()}", file=sys.stderr)
    print(f"CUDA version: {torch.version.cuda}", file=sys.stderr)
    if torch.cuda.is_available():
        print(f"GPU available: {torch.cuda.get_device_name(0)}", file=sys.stderr)
    else:
        print("WARNING: No GPU detected - CellSAM will run on CPU and may be very slow.", file=sys.stderr)

    # Extract and format channels
    img = extract_channels(
        args.tiff,
        args.nuclear_channel,
        args.membrane_channel,
        args.compartment
    )

    # Run CellSAM segmentation
    mask = cellsam_pipeline(
        img,
        bbox_threshold=args.bbox_threshold,
        block_size=args.block_size,
        overlap=args.overlap,
        iou_depth=args.iou_depth,
        iou_threshold=args.iou_threshold,
        use_wsi=args.use_wsi,
        gauge_cell_size=args.gauge_cell_size,
        low_contrast_enhancement=args.low_contrast_enhancement,
        model_path=args.model_path
    )

    return mask


def main():
    '''Main execution function.'''
    # Download/verify model weights
    download_model_weights()

    # Parse arguments
    args = parse_arguments()

    # Run segmentation
    mask = run_segmentation(args)

    # cellSAM turns a block that raises into an all-zero block and carries on,
    # so a run where every block failed produces an empty mask and exit 0.
    if mask.max() == 0:
        sys.exit(
            "ERROR: CellSAM found no cells in the image. If the log shows "
            "'Error segmenting chunk' lines, segmentation failed in every block; "
            "otherwise check the nuclear/membrane channel names and cellsam_bbox_threshold."
        )

    if args.remove_border_cells:
        mask = remove_border_cells(mask)

    # Remove small cells if threshold set
    if args.min_area > 0:
        mask = remove_small_cells(mask, args.min_area)

    if mask.max() == 0:
        sys.exit(
            "ERROR: no cells left after removing border cells and cells below "
            f"cellsam_min_area ({args.min_area} px)."
        )

    # Save segmentation mask
    tifffile.imwrite(
        args.output, mask.astype(np.uint32),
        compression=MASK_COMPRESSION, compressionargs=MASK_COMPRESSION_ARGS,
    )
    print(f"Segmentation mask saved to {args.output}")


if __name__ == "__main__":
    main()
