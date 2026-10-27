"""Unit tests for ``bin/expand_cells.py``.

The strip-by-strip expansion must give the same mask as expanding the whole
image at once (``skimage.segmentation.expand_labels``), and with a nuclear mask
it must leave nuclei that have no whole-cell label untouched.

Pure python -- no pipeline run required.
"""

import numpy as np
import pytest
import tifffile
from skimage.segmentation import expand_labels

from expand_cells import expand_block, expand_mask, main


def _random_mask(shape, n_cells, radius, seed):
    """Label mask of ``n_cells`` discs at random positions; later discs overwrite earlier ones."""
    rng = np.random.default_rng(seed)
    mask = np.zeros(shape, dtype=np.uint32)
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    for label in range(1, n_cells + 1):
        cy, cx = rng.integers(0, shape[0]), rng.integers(0, shape[1])
        mask[(yy - cy) ** 2 + (xx - cx) ** 2 <= radius ** 2] = label
    return mask


@pytest.mark.parametrize("distance_px", [1, 3, 6])
@pytest.mark.parametrize("strip_rows", [7, 50, 1000])
def test_strips_match_whole_image_expansion(distance_px, strip_rows):
    mask = _random_mask((211, 173), n_cells=60, radius=4, seed=distance_px)
    out = np.zeros_like(mask)
    expand_mask(mask, distance_px, out, strip_rows=strip_rows)
    np.testing.assert_array_equal(out, expand_labels(mask, distance_px))


def test_growth_is_limited_to_distance():
    mask = np.zeros((41, 41), dtype=np.uint32)
    mask[20, 20] = 5
    out = expand_block(mask, 3)
    yy, xx = np.nonzero(out)
    assert np.all(np.hypot(yy - 20, xx - 20) <= 3)
    assert set(np.unique(out)) == {0, 5}
    assert np.count_nonzero(out) == np.count_nonzero(expand_labels(mask, 3))


def test_cells_stop_where_they_meet():
    mask = np.zeros((1, 21), dtype=np.uint32)
    mask[0, 5] = 1
    mask[0, 15] = 2
    out = expand_block(mask, 10)
    assert list(out[0, :10]) == [1] * 10
    assert list(out[0, 11:]) == [2] * 10


def test_nuclei_without_a_cell_are_not_taken():
    whole_cell = np.zeros((30, 30), dtype=np.uint32)
    whole_cell[10:20, 2:8] = 1
    nuclear = np.zeros_like(whole_cell)
    nuclear[12:18, 4:7] = 1      # nucleus inside cell 1
    nuclear[12:18, 9:12] = 2     # nucleus with no whole-cell label, next to cell 1
    out = np.zeros_like(whole_cell)
    expand_mask(whole_cell, 5, out, nuclear=nuclear)
    assert np.all(out[nuclear == 2] == 0)
    assert np.all(out[whole_cell == 1] == 1)
    # Cell 1 still grows into free pixels around the unmatched nucleus.
    assert out[11, 10] == 1


def test_empty_mask_stays_empty():
    mask = np.zeros((20, 20), dtype=np.uint16)
    out = np.ones_like(mask)
    expand_mask(mask, 4, out, strip_rows=8)
    assert not out.any()


def test_cli_round_trip(tmp_path):
    mask = _random_mask((120, 90), n_cells=25, radius=3, seed=1).astype(np.uint16)
    wc_path = tmp_path / "wc.tiff"
    tifffile.imwrite(wc_path, mask, compression="zlib")
    out_path = tmp_path / "out.tiff"
    main(["--whole-cell-mask", str(wc_path), "--distance-px", "3",
          "--output", str(out_path), "--strip-rows", "16"])
    out = tifffile.imread(out_path)
    assert out.dtype == np.uint16
    assert out.shape == mask.shape
    # Uncompressed and contiguous, so cellmeasurement can memory-map it.
    assert tifffile.memmap(out_path, mode="r").shape == mask.shape
    assert np.count_nonzero(out) > np.count_nonzero(mask)
    assert set(np.unique(out)) == set(np.unique(mask))


def test_cli_rejects_zero_distance(tmp_path):
    wc_path = tmp_path / "wc.tiff"
    tifffile.imwrite(wc_path, np.zeros((5, 5), dtype=np.uint32))
    with pytest.raises(SystemExit):
        main(["--whole-cell-mask", str(wc_path), "--distance-px", "0",
              "--output", str(tmp_path / "o.tiff")])
