"""Unit tests for ``bin/aviti_plate_common.py`` -- the shared plate-layout
geometry used by both the plate image assembler and the plate GeoJSON merge.

Pure-python/pure-arithmetic, no image or Nextflow dependencies.
"""

import pytest

from aviti_plate_common import (
    compute_plate_layout,
    find_placement,
    microns_to_px,
    parse_well_id,
    read_layout_csv,
    write_layout_csv,
)


# --- well id parsing -------------------------------------------------------


def test_parse_well_id_uses_letter_as_column_and_number_as_row():
    # This deliberately INVERTS the standard microplate convention: on this
    # plate A1 is top-left, A2 sits directly below A1, and B1 is in the next
    # column across. Getting it backwards transposes the whole plate, so pin
    # both axes explicitly.
    assert parse_well_id("A1") == (0, 0)
    assert parse_well_id("A2") == (0, 1)  # A2 is BELOW A1, not to its right
    assert parse_well_id("B1") == (1, 0)  # B1 is the NEXT COLUMN, not below
    assert parse_well_id("B2") == (1, 1)
    assert parse_well_id("D3") == (3, 2)


def test_parse_well_id_handles_multi_digit_and_multi_letter():
    assert parse_well_id("A10") == (0, 9)
    assert parse_well_id("Z1") == (25, 0)
    assert parse_well_id("AA1") == (26, 0)


def test_parse_well_id_is_case_insensitive_and_strips():
    assert parse_well_id("  b2 ") == (1, 1)
    assert parse_well_id("b2") == (1, 1)


@pytest.mark.parametrize("bad", ["1A", "", "A", "A1B", "-1", "A0"])
def test_parse_well_id_rejects_malformed_ids(bad):
    with pytest.raises(ValueError):
        parse_well_id(bad)


# --- plate layout ----------------------------------------------------------


def _synthetic_wells(height=100, width=200):
    """
    An arbitrary-sized (4 columns x 3 rows) set of well ids, purely to
    exercise compute_plate_layout()'s general multi-row/multi-column
    placement. This is NOT modelling any real AVITI plate's dimensions --
    the function itself makes no assumption about well count, so these
    numbers are picked only to be "more than one" in each axis.
    """
    return [
        (f"{letter}{number}", height, width)
        for letter in "ABCD"
        for number in (1, 2, 3)
    ]


def test_compute_plate_layout_places_wells_by_column_and_row():
    placements, canvas_h, canvas_w = compute_plate_layout(
        _synthetic_wells(), gap_px=10, align=0
    )

    assert len(placements) == 12
    assert {p.col for p in placements} == {0, 1, 2, 3}
    assert {p.row for p in placements} == {0, 1, 2}

    # 4 columns of 200 px + 3 gaps; 3 rows of 100 px + 2 gaps.
    assert canvas_w == 4 * 200 + 3 * 10
    assert canvas_h == 3 * 100 + 2 * 10

    a1 = find_placement(placements, "A1")
    assert (a1.x0, a1.y0) == (0, 0)
    # A2 directly below A1 (same column), B1 directly right of A1 (same row).
    assert find_placement(placements, "A2")._replace(well="") == a1._replace(
        well="", row=1, y0=110
    )
    assert find_placement(placements, "B1").x0 == 210
    assert find_placement(placements, "B1").y0 == 0


def test_compute_plate_layout_aligns_origins_and_packs_tightly_when_disabled():
    aligned, _h, _w = compute_plate_layout(_synthetic_wells(), gap_px=10, align=512)
    assert all(p.x0 % 512 == 0 for p in aligned)
    assert all(p.y0 % 512 == 0 for p in aligned)

    tight, _h, _w = compute_plate_layout(_synthetic_wells(), gap_px=10, align=0)
    assert find_placement(tight, "B1").x0 == 210


def test_compute_plate_layout_never_overlaps_when_wells_differ_in_size():
    # A well that lost a tile stitches smaller; columns/rows are sized to the
    # largest member so the ragged well leaves a gap rather than overlapping.
    wells = _synthetic_wells()
    wells[0] = ("A1", 60, 150)
    wells[5] = ("B3", 130, 260)

    placements, _h, _w = compute_plate_layout(wells, gap_px=10, align=0)

    for i, first in enumerate(placements):
        for second in placements[i + 1:]:
            overlap_x = (
                first.x0 < second.x0 + second.width
                and second.x0 < first.x0 + first.width
            )
            overlap_y = (
                first.y0 < second.y0 + second.height
                and second.y0 < first.y0 + first.height
            )
            assert not (overlap_x and overlap_y), f"{first.well} overlaps {second.well}"


def test_compute_plate_layout_leaves_holes_for_unprocessed_wells():
    # Only two wells processed; each keeps its own plate position and the
    # missing ones simply leave holes.
    placements, _h, _w = compute_plate_layout(
        [("A1", 100, 200), ("C2", 100, 200)], gap_px=10, align=0
    )
    assert [p.well for p in placements] == ["A1", "C2"]
    assert find_placement(placements, "A1").col == 0
    assert find_placement(placements, "C2").col == 2
    assert find_placement(placements, "C2").row == 1
    # Origins are cumulative over *occupied* columns only, so C2 does not carry
    # a leading empty column B.
    assert find_placement(placements, "C2").x0 == 210


def test_compute_plate_layout_rejects_an_empty_plate():
    with pytest.raises(ValueError):
        compute_plate_layout([], gap_px=10)


# --- layout CSV round trip -------------------------------------------------


def test_layout_csv_round_trip(tmp_path):
    placements, _h, _w = compute_plate_layout(_synthetic_wells(), gap_px=10, align=512)
    path = tmp_path / "layout.csv"

    write_layout_csv(placements, path)
    assert read_layout_csv(path) == placements


def test_read_layout_csv_rejects_an_empty_file(tmp_path):
    path = tmp_path / "layout.csv"
    path.write_text("well,col,row,x0,y0,width,height\n")
    with pytest.raises(ValueError):
        read_layout_csv(path)


# --- unit conversion -------------------------------------------------------


def test_microns_to_px_matches_the_stitch_helper():
    assert microns_to_px(500.0, 0.48) == 1042
    assert microns_to_px(32.0, 0.48) == 67
    assert microns_to_px(0.0, 0.48) == 0
