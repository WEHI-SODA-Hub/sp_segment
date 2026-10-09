#!/usr/bin/env python
'''
Module      : aviti_plate_common
Description : Shared plate-layout geometry for the AVITI plate-level
              artefacts: the all-wells pyramidal OME-TIFF
              (aviti_assemble_plate_image.py) and the merged plate GeoJSON
              (aviti_merge_plate_geojson.py). Both must place every well at
              byte-identical coordinates or the cell polygons land off the
              cells they describe, so the layout maths lives here once and
              is passed between the two as a layout CSV.

              **Well naming convention.** The AVITI plate layout used here is
              NOT the standard microplate convention. The letter is the
              COLUMN and the number is the ROW -- illustrated below with a
              minimal, arbitrary-sized example (this code makes no
              assumption anywhere about how many wells, columns, or rows an
              actual plate has: compute_plate_layout() derives the grid
              purely from whichever well ids are actually present):

                  A1  B1
                  A2  B2

              so A1 is top-left, A2 sits directly *below* A1 (same column,
              next row), and B1 is in the next column across (same row,
              first row). A standard microplate would have A2 to the right
              of A1 and B1 below it instead. This mirrors how the wells are
              presented for these runs; getting it backwards transposes the
              whole plate, which is why parse_well_id() is covered by an
              explicit regression test.

              Wells are placed on a grid rather than by stage coordinates:
              the inter-well stage jumps are an order of magnitude larger
              than the intra-well tile steps, and only the within-well
              placement (bin/aviti_stitch.py) is coordinate-derived.
Copyright   : (c) WEHI SODA Hub, 2026
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import csv
import re
from collections import namedtuple
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

# One well's placement on the plate canvas. x0/y0 are the pixel coordinates of
# the well's top-left corner; width/height are that well's own stitched size
# (wells can differ when a tile is missing).
PlacedWell = namedtuple("PlacedWell", "well col row x0 y0 width height")

WELL_ID_PATTERN = re.compile(r"^([A-Za-z]+)([0-9]+)$")


def parse_well_id(well: str) -> Tuple[int, int]:
    '''
    Map a WellLocation string to (column, row) zero-based indices.

    The letter is the column and the number is the row -- see this module's
    docstring for why that is the opposite of the standard microplate
    convention. Multi-letter prefixes are read as bijective base-26 ("Z" ->
    25, "AA" -> 26), so plates wider than 26 columns still work.
    '''
    match = WELL_ID_PATTERN.match((well or "").strip())
    if match is None:
        raise ValueError(
            f"Could not parse well id {well!r}: expected letters followed by digits, e.g. 'A1'."
        )
    letters, digits = match.groups()

    col = 0
    for char in letters.upper():
        col = col * 26 + (ord(char) - ord("A") + 1)
    col -= 1

    row = int(digits) - 1
    if row < 0:
        raise ValueError(f"Could not parse well id {well!r}: well numbers start at 1.")
    return col, row


def _round_up(value: int, multiple: int) -> int:
    if multiple <= 1:
        return value
    return ((value + multiple - 1) // multiple) * multiple


def compute_plate_layout(
    wells: Sequence[Tuple[str, int, int]],
    gap_px: int,
    align: int = 512,
) -> Tuple[List[PlacedWell], int, int]:
    '''
    Place wells on the plate grid and return (placements, canvas_h, canvas_w).

    ``wells`` is a sequence of (well_id, height, width) for each well actually
    processed. Column widths are the maximum width over the wells in that
    column and row heights the maximum height over that row, so a well that
    stitched smaller (a missing tile) leaves a ragged edge rather than
    overlapping its neighbour.

    Origins advance cumulatively with ``gap_px`` between neighbours and are
    rounded up to a multiple of ``align``; ``align=0`` or ``1`` packs tightly.
    Alignment is a read-efficiency convenience at full resolution only -- it
    does not hold at reduced pyramid levels, so consumers must never assume
    an output tile intersects at most one well.

    The grid may be sparse: wells that were not processed simply leave holes,
    and the remaining wells keep their true plate positions.
    '''
    if not wells:
        raise ValueError("No wells to place on the plate.")

    _Parsed = namedtuple("_Parsed", "well col row height width")
    parsed = []
    for well, height, width in wells:
        col, row = parse_well_id(well)
        parsed.append(_Parsed(well, col, row, int(height), int(width)))

    col_width: dict = {}
    row_height: dict = {}
    for entry in parsed:
        col_width[entry.col] = max(col_width.get(entry.col, 0), entry.width)
        row_height[entry.row] = max(row_height.get(entry.row, 0), entry.height)

    # Cumulative origins over the occupied columns/rows only, so a plate whose
    # first column was not processed does not carry a leading empty column.
    col_x0: dict = {}
    cursor = 0
    for col in sorted(col_width):
        col_x0[col] = _round_up(cursor, align)
        cursor = col_x0[col] + col_width[col] + gap_px

    row_y0: dict = {}
    cursor = 0
    for row in sorted(row_height):
        row_y0[row] = _round_up(cursor, align)
        cursor = row_y0[row] + row_height[row] + gap_px

    placements = [
        PlacedWell(
            well=entry.well, col=entry.col, row=entry.row,
            x0=col_x0[entry.col], y0=row_y0[entry.row],
            width=entry.width, height=entry.height,
        )
        for entry in parsed
    ]
    placements.sort(key=lambda p: (p.col, p.row))

    canvas_w = max(p.x0 + p.width for p in placements)
    canvas_h = max(p.y0 + p.height for p in placements)
    return placements, canvas_h, canvas_w


LAYOUT_FIELDNAMES = ["well", "col", "row", "x0", "y0", "width", "height"]


def write_layout_csv(placements: Sequence[PlacedWell], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=LAYOUT_FIELDNAMES)
        writer.writeheader()
        for placed in placements:
            writer.writerow(placed._asdict())


def read_layout_csv(path: Path) -> List[PlacedWell]:
    placements = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            placements.append(PlacedWell(
                well=row["well"].strip(),
                col=int(row["col"]), row=int(row["row"]),
                x0=int(row["x0"]), y0=int(row["y0"]),
                width=int(row["width"]), height=int(row["height"]),
            ))
    if not placements:
        raise ValueError(f"Plate layout {path} contains no wells.")
    return placements


def microns_to_px(value_microns: float, pixel_size_microns: float) -> int:
    '''
    Convert a distance in microns to whole pixels. Mirrors the helper of the
    same name in bin/aviti_stitch.py so the plate and well gaps are computed
    identically.
    '''
    return round(float(value_microns) / float(pixel_size_microns))


def find_placement(placements: Sequence[PlacedWell], well: str) -> Optional[PlacedWell]:
    for placed in placements:
        if placed.well == well:
            return placed
    return None
