"""Unit tests for ``bin/aviti_assemble_plate_image.py``.

These cover the streaming/compositing logic that keeps the plate assembler's
memory bounded, and pin the tile-iterator contract that tifffile's
iterator-based writes depend on (channel-major, then tile-row, then
tile-column, every item exactly tile-sized).

The end-to-end test also acts as the on-container verification that
``ome=True`` + ``subifds`` + an iterator really does produce a readable
pyramid with the physical pixel size QuPath needs.
"""

import numpy as np
import pytest
import tifffile

from aviti_assemble_plate_image import (
    check_channels_consistent,
    downsample2x,
    downsample_region,
    iter_level_tiles,
    open_well_readers,
    pyramid_shapes,
    read_plate_region,
    read_well_rows,
)
from aviti_plate_common import PlacedWell


# --- downsampling ----------------------------------------------------------


def test_downsample2x_takes_the_rounded_block_mean():
    block = np.array([[0, 2, 10, 20], [4, 6, 30, 40]], dtype=np.uint16)
    # means: (0+2+4+6)/4 = 3, (10+20+30+40)/4 = 25
    np.testing.assert_array_equal(downsample2x(block), np.array([[3, 25]], dtype=np.uint16))


def test_downsample2x_preserves_dtype_and_does_not_overflow_at_uint16_max():
    block = np.full((2, 2), 65535, dtype=np.uint16)
    out = downsample2x(block)
    assert out.dtype == np.uint16
    np.testing.assert_array_equal(out, np.array([[65535]], dtype=np.uint16))


def test_downsample2x_drops_a_trailing_odd_row_and_column():
    block = np.ones((3, 3), dtype=np.uint16)
    assert downsample2x(block).shape == (1, 1)


# --- compositing wells into plate regions ---------------------------------


def _two_well_plate():
    """Two 4x4 single-channel wells side by side with a 2 px gap."""
    readers = {
        "A1": np.full((1, 4, 4), 11, dtype=np.uint16),
        "B1": np.full((1, 4, 4), 22, dtype=np.uint16),
    }
    placements = [
        PlacedWell("A1", 0, 0, x0=0, y0=0, width=4, height=4),
        PlacedWell("B1", 1, 0, x0=6, y0=0, width=4, height=4),
    ]
    return readers, placements


def test_read_plate_region_inside_one_well():
    readers, placements = _two_well_plate()
    out = read_plate_region(readers, placements, 1, 1, 2, 2, 1, np.uint16)
    np.testing.assert_array_equal(out, np.full((1, 2, 2), 11, dtype=np.uint16))


def test_read_plate_region_zero_fills_the_gap_between_wells():
    readers, placements = _two_well_plate()
    out = read_plate_region(readers, placements, 0, 4, 1, 2, 1, np.uint16)
    np.testing.assert_array_equal(out, np.zeros((1, 1, 2), dtype=np.uint16))


def test_read_plate_region_composites_across_two_wells():
    # The case that tile-alignment does NOT protect against: at reduced
    # pyramid levels a single output tile can straddle wells, so the region
    # reader must pull from every intersecting well.
    readers, placements = _two_well_plate()
    out = read_plate_region(readers, placements, 0, 2, 1, 8, 1, np.uint16)
    expected = np.array([[[11, 11, 0, 0, 22, 22, 22, 22]]], dtype=np.uint16)
    np.testing.assert_array_equal(out, expected)


def test_read_plate_region_zero_fills_past_the_canvas_edge():
    readers, placements = _two_well_plate()
    out = read_plate_region(readers, placements, 2, 8, 4, 4, 1, np.uint16)
    expected = np.zeros((1, 4, 4), dtype=np.uint16)
    expected[0, :2, :2] = 22
    np.testing.assert_array_equal(out, expected)


def test_downsample_region_matches_a_direct_full_read_and_halve():
    rng = np.random.default_rng(0)
    readers = {"A1": rng.integers(0, 4000, size=(2, 16, 16), dtype=np.uint16)}
    placements = [PlacedWell("A1", 0, 0, x0=0, y0=0, width=16, height=16)]

    direct = downsample2x(downsample2x(readers["A1"]))
    chunked = downsample_region(
        readers, placements, level=2, out_y=0, out_x=0, out_h=4, out_w=4,
        n_channels=2, dtype=np.uint16,
        # Tiny budget forces multiple sub-chunks -- the mechanism that keeps
        # peak memory bounded at high pyramid levels.
        read_budget_bytes=64,
    )
    np.testing.assert_array_equal(chunked, direct)


# --- pyramid shapes --------------------------------------------------------


def test_pyramid_shapes_halves_until_within_min_level_size():
    assert pyramid_shapes(4096, 8192, max_levels=8, min_level_size=1024) == [
        (4096, 8192), (2048, 4096), (1024, 2048), (512, 1024),
    ]


def test_pyramid_shapes_returns_a_single_level_for_a_small_canvas():
    assert pyramid_shapes(100, 100, max_levels=8, min_level_size=1024) == [(100, 100)]


def test_pyramid_shapes_respects_max_levels():
    assert len(pyramid_shapes(100000, 100000, max_levels=3, min_level_size=1)) == 4


# --- the tile iterator contract -------------------------------------------


def test_iter_level_tiles_yields_channel_major_row_major_full_size_tiles():
    rng = np.random.default_rng(1)
    well = rng.integers(0, 4000, size=(2, 6, 6), dtype=np.uint16)
    readers = {"A1": well}
    placements = [PlacedWell("A1", 0, 0, x0=0, y0=0, width=6, height=6)]
    tile = 4

    tiles = list(iter_level_tiles(
        readers, placements, 0, (2, 6, 6), tile, np.uint16, 1 << 20
    ))

    # 2 channels x 2 tile-rows x 2 tile-cols
    assert len(tiles) == 2 * 2 * 2
    # tifffile requires every item to match the tile shape exactly; edge tiles
    # are zero-padded by us, not by tifffile, when writing from an iterator.
    assert all(t.shape == (tile, tile) for t in tiles)

    # Reassembling in the documented order must reproduce the source, which is
    # what pins the ordering contract.
    rebuilt = np.zeros((2, 8, 8), dtype=np.uint16)
    index = 0
    for channel in range(2):
        for y in range(0, 8, tile):
            for x in range(0, 8, tile):
                rebuilt[channel, y:y + tile, x:x + tile] = tiles[index]
                index += 1
    np.testing.assert_array_equal(rebuilt[:, :6, :6], well)
    # Padding is zero, not wrapped data.
    assert rebuilt[:, 6:, :].sum() == 0
    assert rebuilt[:, :, 6:].sum() == 0


# --- manifest / reader plumbing -------------------------------------------


def _write_well_image(path, array, channel_names):
    """Write a well image the way bin/aviti_stitch.py writes one.

    aviti_stitch.py's own call omits ``photometric=``; the container it runs
    in is pinned to tifffile 2025.6.11, which tolerates that. A newer
    tifffile (as this dev/test environment happens to resolve) requires
    ``photometric="minisblack"`` to accept a >2-channel, non-RGB CYX stack --
    passed explicitly here so these fixtures aren't at the mercy of which
    tifffile `uv` picks, independent of what's under test.
    """
    tifffile.imwrite(
        path, array, ome=True, photometric="minisblack",
        metadata={"axes": "CYX", "Channel": {"Name": channel_names}},
    )


def test_open_well_readers_memmaps_stitch_style_well_images(tmp_path):
    # Regression guard: the plate assembler streams from memory-mapped well
    # images, which only works while aviti_stitch.py writes them uncompressed.
    path = tmp_path / "w.tif"
    _write_well_image(path, np.zeros((2, 8, 8), dtype=np.uint16), ["Nucleus", "Cell-Membrane"])

    readers = open_well_readers([{"well": "A1", "image": str(path)}])
    assert readers["A1"].shape == (2, 8, 8)
    assert isinstance(readers["A1"], np.memmap)


def test_read_well_rows_orders_by_plate_position_not_manifest_order(tmp_path):
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("well,image\nB1,b1.tif\nA2,a2.tif\nA1,a1.tif\n")
    assert [r["well"] for r in read_well_rows(manifest)] == ["A1", "A2", "B1"]


def test_check_channels_consistent_rejects_a_channel_count_mismatch(tmp_path):
    one = tmp_path / "one.tif"
    two = tmp_path / "two.tif"
    _write_well_image(one, np.zeros((2, 4, 4), dtype=np.uint16), ["Nucleus", "Cell-Membrane"])
    _write_well_image(two, np.zeros((3, 4, 4), dtype=np.uint16), ["Nucleus", "Cell-Membrane", "Actin"])

    rows = [{"well": "A1", "image": str(one)}, {"well": "B1", "image": str(two)}]
    with pytest.raises(ValueError, match="channel"):
        check_channels_consistent(rows, open_well_readers(rows))


# --- end to end ------------------------------------------------------------


def test_end_to_end_writes_a_readable_pyramidal_ome_tiff(tmp_path):
    from typer.testing import CliRunner

    from aviti_assemble_plate_image import app

    names = ["Nucleus", "Cell-Membrane"]
    a1 = np.full((2, 64, 64), 100, dtype=np.uint16)
    b1 = np.full((2, 64, 64), 200, dtype=np.uint16)
    _write_well_image(tmp_path / "a1.tif", a1, names)
    _write_well_image(tmp_path / "b1.tif", b1, names)

    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        f"well,image\nA1,{tmp_path / 'a1.tif'}\nB1,{tmp_path / 'b1.tif'}\n"
    )
    out_image = tmp_path / "plate.ome.tif"
    out_layout = tmp_path / "plate_layout.csv"
    out_overview = tmp_path / "plate_overview.tif"

    result = CliRunner().invoke(app, [
        str(manifest),
        "--output-image", str(out_image),
        "--output-layout", str(out_layout),
        "--output-overview", str(out_overview),
        "--well-gap-microns", "16",
        "--pixel-size-microns", "0.5",   # gap becomes 32 px
        "--tile-size", "16",
        "--min-level-size", "16",
        "--well-align", "16",
    ])
    assert result.exit_code == 0, result.output

    with tifffile.TiffFile(out_image) as tif:
        assert tif.is_ome
        assert tif.is_bigtiff
        series = tif.series[0]
        assert series.axes == "CYX"
        # A 64x160 canvas halves until BOTH dims are <= min_level_size (16):
        # (64,160) -> (32,80) -> (16,40) -> (8,20) -> (4,10); width lags height
        # since the canvas isn't square, so this needs one more halving than
        # height alone would.
        assert len(series.levels) == 5
        assert series.levels[0].shape == (2, 64, 160)
        assert series.levels[-1].shape == (2, 4, 10)

        level0 = series.levels[0].asarray()
        # A1 at x 0-63, 32 px gap, B1 at x 96-159.
        np.testing.assert_array_equal(level0[:, :, :64], a1)
        np.testing.assert_array_equal(level0[:, :, 64:96], np.zeros((2, 64, 32), dtype=np.uint16))
        np.testing.assert_array_equal(level0[:, :, 96:], b1)

        ome = tif.ome_metadata
        assert "Nucleus" in ome and "Cell-Membrane" in ome
        assert 'PhysicalSizeX="0.5"' in ome

    assert out_layout.exists()
    with tifffile.TiffFile(out_overview) as tif:
        assert tif.series[0].shape == (2, 4, 10)
