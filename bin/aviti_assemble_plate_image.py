#!/usr/bin/env python
'''
Module      : aviti_assemble_plate_image
Description : Assembles every stitched AVITI well image into one plate-level,
              tiled, pyramidal BigTIFF OME-TIFF that QuPath can open and zoom
              smoothly, plus the plate layout CSV that
              aviti_merge_plate_geojson.py uses to place cell polygons in the
              same coordinate frame.

              The plate canvas is far too large to hold in memory -- a 12-well
              plate of 12 tiles each is roughly 14 GB at full resolution for
              two channels -- so nothing here ever allocates it. Each source
              well image is memory-mapped (bin/aviti_stitch.py writes them
              uncompressed, one strip per plane, so they are memmappable) and
              the output is streamed tile by tile straight into
              tifffile.TiffWriter. Reduced pyramid levels are produced by
              re-reading the same memmaps in bounded row sub-chunks, so peak
              memory is set by --read-budget-mb and not by the plate size or
              the number of levels.

              Well placement comes from aviti_plate_common (letter = column,
              number = row -- see that module for why that inverts the usual
              microplate convention).
Copyright   : (c) WEHI SODA Hub, 2026
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import csv
import sys
from pathlib import Path
from typing import Annotated, Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import tifffile
import typer

from aviti_plate_common import (
    PlacedWell,
    compute_plate_layout,
    microns_to_px,
    parse_well_id,
    write_layout_csv,
)
from aviti_stitch import get_channel_names

app = typer.Typer(add_completion=False)

# tifffile requires tile dimensions to be multiples of 16.
TILE_MULTIPLE = 16


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def read_well_rows(manifest: Path) -> List[dict]:
    '''
    Read the well manifest (columns: well, image), ordered by plate position.

    Sorting here rather than trusting the manifest order means Nextflow's
    groupTuple emission order -- which follows task completion -- can never
    change the assembled result.
    '''
    with open(manifest, newline="") as fh:
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in csv.DictReader(fh)
        ]
    if not rows:
        raise ValueError(f"No rows found in plate manifest {manifest}")
    return sorted(rows, key=lambda r: parse_well_id(r["well"]))


def open_well_readers(rows: Sequence[dict]) -> Dict[str, np.ndarray]:
    '''
    Memory-map each well's stitched image, normalising a 2-D (single-channel)
    image to (1, H, W) so callers always see CYX.
    '''
    readers: Dict[str, np.ndarray] = {}
    for row in rows:
        path = Path(row["image"])
        try:
            array = tifffile.memmap(path)
        except Exception as exc:
            raise ValueError(
                f"Could not memory-map well image {path}: {exc}. The plate assembler "
                "relies on aviti_stitch.py writing the stitched well image "
                "uncompressed; if compression was added there, this step needs to "
                "read the images eagerly instead."
            ) from exc
        if array.ndim == 2:
            array = array[np.newaxis, :, :]
        if array.ndim != 3:
            raise ValueError(
                f"Expected a 2-D or CYX 3-D well image, got shape {array.shape} from {path}"
            )
        readers[row["well"]] = array
    return readers


def check_channels_consistent(
    rows: Sequence[dict], readers: Dict[str, np.ndarray]
) -> Tuple[int, np.dtype, List[str]]:
    '''
    Every well must agree on channel count, dtype and channel names -- a plate
    whose channel 2 means different things in different wells would be
    silently wrong rather than obviously broken.
    '''
    first_well = rows[0]["well"]
    n_channels = readers[first_well].shape[0]
    dtype = readers[first_well].dtype

    for row in rows:
        array = readers[row["well"]]
        if array.shape[0] != n_channels:
            raise ValueError(
                f"Well {row['well']} has {array.shape[0]} channel(s) but well "
                f"{first_well} has {n_channels}; cannot assemble a coherent plate."
            )
        if array.dtype != dtype:
            raise ValueError(
                f"Well {row['well']} has dtype {array.dtype} but well {first_well} "
                f"has {dtype}; cannot assemble a coherent plate."
            )

    channel_names = get_channel_names(Path(rows[0]["image"]), n_channels)
    for row in rows[1:]:
        names = get_channel_names(Path(row["image"]), n_channels)
        if names != channel_names:
            raise ValueError(
                f"Well {row['well']} has channel names {names} but well {first_well} "
                f"has {channel_names}; cannot assemble a coherent plate."
            )
    return n_channels, dtype, channel_names


def read_plate_region(
    readers: Dict[str, np.ndarray],
    placements: Sequence[PlacedWell],
    y: int,
    x: int,
    height: int,
    width: int,
    n_channels: int,
    dtype: np.dtype,
) -> np.ndarray:
    '''
    Composite the plate rectangle (y, x, height, width) from every well that
    intersects it, zero-filling gaps and anything past the canvas edge.

    This deliberately handles a region spanning several wells. Aligning well
    origins to the output tile size only guarantees one-well-per-tile at full
    resolution: at pyramid level L the origin becomes x0 / 2**L, so a single
    reduced-resolution tile can straddle two or four wells.
    '''
    out = np.zeros((n_channels, height, width), dtype=dtype)

    for placed in placements:
        src_y0 = max(y, placed.y0)
        src_y1 = min(y + height, placed.y0 + placed.height)
        src_x0 = max(x, placed.x0)
        src_x1 = min(x + width, placed.x0 + placed.width)
        if src_y0 >= src_y1 or src_x0 >= src_x1:
            continue

        array = readers[placed.well]
        out[
            :,
            src_y0 - y: src_y1 - y,
            src_x0 - x: src_x1 - x,
        ] = array[
            :,
            src_y0 - placed.y0: src_y1 - placed.y0,
            src_x0 - placed.x0: src_x1 - placed.x0,
        ]

    return out


def downsample2x(array: np.ndarray) -> np.ndarray:
    '''
    Rounded 2x2 block mean, computed in uint32 so a uint16 input cannot
    overflow and no float conversion is needed. Trailing odd rows/columns are
    dropped; callers always pass even-sized blocks.
    '''
    even = array[..., : array.shape[-2] // 2 * 2, : array.shape[-1] // 2 * 2]
    wide = even.astype(np.uint32)
    total = (
        wide[..., 0::2, 0::2] + wide[..., 0::2, 1::2]
        + wide[..., 1::2, 0::2] + wide[..., 1::2, 1::2]
    )
    return ((total + 2) >> 2).astype(array.dtype)


def downsample_region(
    readers: Dict[str, np.ndarray],
    placements: Sequence[PlacedWell],
    level: int,
    out_y: int,
    out_x: int,
    out_h: int,
    out_w: int,
    n_channels: int,
    dtype: np.dtype,
    read_budget_bytes: int,
) -> np.ndarray:
    '''
    Produce an (n_channels, out_h, out_w) block at ``level`` by reading the
    corresponding full-resolution region and halving it ``level`` times.

    The source region is read in row sub-chunks whose height is a multiple of
    2**level, so peak memory stays near ``read_budget_bytes`` regardless of
    level. Without this an output tile at level 6 would pull a 2 GB block per
    channel.
    '''
    if level == 0:
        return read_plate_region(
            readers, placements, out_y, out_x, out_h, out_w, n_channels, dtype
        )

    factor = 1 << level
    src_y = out_y * factor
    src_x = out_x * factor
    src_w = out_w * factor

    itemsize = np.dtype(dtype).itemsize
    row_bytes = max(1, n_channels * src_w * itemsize)
    chunk_rows = max(factor, (read_budget_bytes // row_bytes) // factor * factor)

    out = np.zeros((n_channels, out_h, out_w), dtype=dtype)
    for out_row in range(0, out_h, chunk_rows // factor):
        out_rows = min(chunk_rows // factor, out_h - out_row)
        block = read_plate_region(
            readers, placements,
            src_y + out_row * factor, src_x,
            out_rows * factor, src_w,
            n_channels, dtype,
        )
        for _ in range(level):
            block = downsample2x(block)
        out[:, out_row: out_row + out_rows, :] = block

    return out


def pyramid_shapes(
    height: int, width: int, max_levels: int, min_level_size: int
) -> List[Tuple[int, int]]:
    '''
    Full-resolution shape followed by successive halvings, stopping once both
    dimensions fit within ``min_level_size`` or ``max_levels`` is reached.
    '''
    shapes = [(height, width)]
    while len(shapes) <= max_levels:
        h, w = shapes[-1]
        if h <= min_level_size and w <= min_level_size:
            break
        if h <= 1 and w <= 1:
            break
        shapes.append((max(1, (h + 1) // 2), max(1, (w + 1) // 2)))
    return shapes


def iter_level_tiles(
    readers: Dict[str, np.ndarray],
    placements: Sequence[PlacedWell],
    level: int,
    shape: Tuple[int, int, int],
    tile: int,
    dtype: np.dtype,
    read_budget_bytes: int,
) -> Iterator[np.ndarray]:
    '''
    Yield every tile of one pyramid level in the order tifffile's iterator
    writes expect: channel-major, then tile-row, then tile-column. Each item
    is exactly (tile, tile) -- tifffile zero-pads nothing for us when an
    iterator is used, so edge tiles are padded here.
    '''
    n_channels, height, width = shape
    empty = np.zeros((tile, tile), dtype=dtype)

    for channel in range(n_channels):
        for y in range(0, height, tile):
            for x in range(0, width, tile):
                out_h = min(tile, height - y)
                out_w = min(tile, width - x)
                block = downsample_region(
                    readers, placements, level, y, x, out_h, out_w,
                    n_channels, dtype, read_budget_bytes,
                )[channel]
                if out_h == tile and out_w == tile:
                    yield block
                else:
                    padded = empty.copy()
                    padded[:out_h, :out_w] = block
                    yield padded


@app.command()
def main(
    manifest: Annotated[Path, typer.Argument(
        exists=True,
        help="CSV manifest of the stitched per-well images (columns: well, image)."
    )],
    output_image: Annotated[Path, typer.Option(
        help="Path to write the pyramidal plate OME-TIFF."
    )],
    output_layout: Annotated[Path, typer.Option(
        help="Path to write the plate layout CSV consumed by the GeoJSON merge."
    )],
    output_overview: Annotated[Optional[Path], typer.Option(
        help="Optional path for a small single-resolution overview image (the "
             "smallest pyramid level), for quick-look QC."
    )] = None,
    well_gap_microns: Annotated[float, typer.Option(
        help="Visual gap inserted between adjacent wells."
    )] = 500.0,
    pixel_size_microns: Annotated[float, typer.Option(
        help="Pixel size written into the OME metadata as PhysicalSizeX/Y, which is "
             "what QuPath shows as its scale bar, and used to convert the well gap."
    )] = 0.48,
    tile_size: Annotated[int, typer.Option(
        help="Output TIFF tile edge in pixels; must be a multiple of 16."
    )] = 512,
    max_levels: Annotated[int, typer.Option(
        help="Maximum number of reduced-resolution pyramid levels."
    )] = 8,
    min_level_size: Annotated[int, typer.Option(
        help="Stop adding pyramid levels once both dimensions are at most this."
    )] = 1024,
    compression: Annotated[str, typer.Option(
        help="TIFF compression for every level ('zlib', 'none', 'lzw')."
    )] = "zlib",
    compression_level: Annotated[int, typer.Option(
        help="Compression level, where the codec supports one."
    )] = 1,
    well_align: Annotated[int, typer.Option(
        help="Round well origins up to a multiple of this many pixels. Purely a "
             "read-efficiency convenience at full resolution; correctness never "
             "depends on it."
    )] = 512,
    read_budget_mb: Annotated[int, typer.Option(
        help="Approximate cap on the source block read for one output tile."
    )] = 512,
):
    '''
    Assemble the stitched per-well images into one pyramidal plate OME-TIFF.
    '''
    if tile_size % TILE_MULTIPLE != 0:
        raise typer.BadParameter(
            f"--tile-size must be a multiple of {TILE_MULTIPLE}, got {tile_size}"
        )

    rows = read_well_rows(manifest)
    readers = open_well_readers(rows)
    n_channels, dtype, channel_names = check_channels_consistent(rows, readers)

    gap_px = microns_to_px(well_gap_microns, pixel_size_microns)
    placements, canvas_h, canvas_w = compute_plate_layout(
        [(row["well"], readers[row["well"]].shape[1], readers[row["well"]].shape[2])
         for row in rows],
        gap_px=gap_px,
        align=well_align,
    )
    write_layout_csv(placements, output_layout)

    modal_size = {(p.height, p.width) for p in placements}
    if len(modal_size) > 1:
        log(
            "WARNING: wells differ in stitched size "
            f"({sorted(modal_size)}); this usually means a well is missing tiles. "
            "Columns/rows are sized to the largest well, so smaller wells leave a "
            "ragged edge rather than overlapping."
        )

    shapes = pyramid_shapes(canvas_h, canvas_w, max_levels, min_level_size)
    n_sub = len(shapes) - 1
    compression_args = {"level": compression_level} if compression not in (None, "none") else None
    read_budget_bytes = max(1, read_budget_mb) * 1024 * 1024

    log(
        f"Assembling {len(placements)} well(s) into a {canvas_w} x {canvas_h} plate "
        f"({n_channels} channel(s), {dtype}), {len(shapes)} pyramid level(s), "
        f"{tile_size} px tiles, {gap_px} px well gap"
    )

    with tifffile.TiffWriter(output_image, bigtiff=True, ome=True) as writer:
        writer.write(
            data=iter_level_tiles(
                readers, placements, 0, (n_channels, canvas_h, canvas_w),
                tile_size, dtype, read_budget_bytes,
            ),
            shape=(n_channels, canvas_h, canvas_w),
            dtype=dtype,
            tile=(tile_size, tile_size),
            photometric="minisblack",
            subifds=n_sub,
            compression=compression,
            compressionargs=compression_args,
            metadata={
                "axes": "CYX",
                "Channel": {"Name": channel_names},
                "PhysicalSizeX": pixel_size_microns,
                "PhysicalSizeXUnit": "µm",
                "PhysicalSizeY": pixel_size_microns,
                "PhysicalSizeYUnit": "µm",
            },
        )
        # Sub-resolutions carry no metadata of their own: OME-XML lives only in
        # the first page, and writing it again here would corrupt the series.
        for level, (level_h, level_w) in enumerate(shapes[1:], start=1):
            writer.write(
                data=iter_level_tiles(
                    readers, placements, level, (n_channels, level_h, level_w),
                    tile_size, dtype, read_budget_bytes,
                ),
                shape=(n_channels, level_h, level_w),
                dtype=dtype,
                tile=(tile_size, tile_size),
                photometric="minisblack",
                subfiletype=1,
                compression=compression,
                compressionargs=compression_args,
            )

    log(f"Plate image written to {output_image}")

    if output_overview is not None:
        overview_level = len(shapes) - 1
        overview_h, overview_w = shapes[overview_level]
        overview = downsample_region(
            readers, placements, overview_level, 0, 0, overview_h, overview_w,
            n_channels, dtype, read_budget_bytes,
        )
        overview_pixel_size = pixel_size_microns * (1 << overview_level)
        tifffile.imwrite(
            output_overview,
            overview,
            ome=True,
            metadata={
                "axes": "CYX",
                "Channel": {"Name": channel_names},
                "PhysicalSizeX": overview_pixel_size,
                "PhysicalSizeXUnit": "µm",
                "PhysicalSizeY": overview_pixel_size,
                "PhysicalSizeYUnit": "µm",
            },
        )
        log(f"Plate overview ({overview_w} x {overview_h}) written to {output_overview}")


if __name__ == "__main__":
    app()
