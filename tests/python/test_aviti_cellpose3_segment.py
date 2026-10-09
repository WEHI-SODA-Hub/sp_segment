"""Unit tests for ``bin/aviti_cellpose3_segment.py`` (nuclear and membrane
modes).

These do not import cellpose/torch (those imports live inside ``main()``
and ``model_input_channels()``), so they run without a GPU or a Cellpose
environment. ``model_input_channels()`` is exercised by injecting a fake
``torch`` module into ``sys.modules`` before it runs its lazy ``import
torch`` -- see ``_install_fake_torch`` below.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

from aviti_cellpose3_segment import (
    binarize,
    build_membrane_stack,
    model_input_channels,
    reconcile_membrane_stack,
    remove_small_cells,
)


# --- nuclear mode: mask post-processing ------------------------------------


def test_remove_small_cells_zeroes_labels_below_min_area():
    mask = np.array([[1, 1, 2], [1, 1, 2], [3, 3, 3]], dtype=np.int32)
    # label 1: area 4, label 2: area 2, label 3: area 3
    out = remove_small_cells(mask, min_area=3)
    expected = np.array([[1, 1, 0], [1, 1, 0], [3, 3, 3]], dtype=np.int32)
    np.testing.assert_array_equal(out, expected)


def test_remove_small_cells_noop_when_min_area_is_zero():
    mask = np.array([[1, 2], [0, 2]], dtype=np.int32)
    out = remove_small_cells(mask, min_area=0)
    np.testing.assert_array_equal(out, mask)


def test_binarize_collapses_instance_labels_to_zero_one():
    # Per Elembio's notebook, Nuclear.tif is a uint8 0/1 presence mask, not
    # instance-labeled -- distinct nucleus IDs (1, 2, 3, ...) must all
    # collapse to 1.
    mask = np.array([[0, 1, 2], [3, 0, 7]], dtype=np.int32)
    out = binarize(mask)
    expected = np.array([[0, 1, 1], [1, 0, 1]], dtype=np.uint8)
    np.testing.assert_array_equal(out, expected)
    assert out.dtype == np.uint8


def test_binarize_handles_a_high_instance_count_without_wraparound():
    # A real AVITI tile can contain thousands of nuclei -- confirming that
    # binarizing (rather than narrowing instance IDs to uint8 directly)
    # avoids the label-collision-via-modulo-wraparound risk that motivated
    # this fix: label values here deliberately exceed uint8's range.
    mask = np.array([[0, 300, 1000, 65000]], dtype=np.int32)
    out = binarize(mask)
    np.testing.assert_array_equal(out, np.array([[0, 1, 1, 1]], dtype=np.uint8))


def test_remove_small_cells_then_binarize_preserves_presence_after_filtering():
    mask = np.array([[1, 1, 2, 3]], dtype=np.int32)  # label areas: 1->2, 2->1, 3->1
    filtered = remove_small_cells(mask, min_area=2)  # drops labels 2 and 3
    out = binarize(filtered)
    np.testing.assert_array_equal(out, np.array([[1, 1, 0, 0]], dtype=np.uint8))


def test_label_mask_retains_instance_ids_that_binary_mask_collapses():
    # This is the pair of artifacts AVITINUCLEARSEGMENT writes from one
    # filtering pass: --output-label keeps remove_small_cells' output as-is
    # (instance IDs intact, for stitching/CELLMEASUREMENT), while --output
    # binarizes it (Elembio's Nuclear.tif viewer/cells2stats contract).
    mask = np.array([[0, 1, 1, 2, 2]], dtype=np.int32)
    filtered = remove_small_cells(mask, min_area=1)

    label_mask = filtered.astype(np.uint16)
    binary_mask = binarize(filtered)

    # The label mask distinguishes nucleus 1 from nucleus 2...
    assert set(np.unique(label_mask)) == {0, 1, 2}
    # ...while the binary mask does not.
    assert set(np.unique(binary_mask)) == {0, 1}


# --- membrane mode: composite construction ----------------------------------


def test_build_membrane_stack_orders_channels_cell_nucleus_actin(tmp_path):
    # Confirmed order: [cell, nucleus, actin] -- a different order from the
    # v4 SAM path's [nucleus, membrane, actin] stack. This is load-bearing:
    # the membrane models are called with channels=None, so Cellpose does no
    # reordering -- whatever order this function stacks in is what the model
    # sees.
    cell = np.full((4, 4), 1, dtype=np.uint16)
    nucleus = np.full((4, 4), 2, dtype=np.uint16)
    actin = np.full((4, 4), 3, dtype=np.uint16)

    cell_tif, nucleus_tif, actin_tif = _write_tiles(tmp_path, cell=cell, nucleus=nucleus, actin=actin)

    stack = build_membrane_stack(cell_tif, nucleus_tif, actin_tif)

    assert stack.shape == (4, 4, 3)
    np.testing.assert_array_equal(stack[..., 0], cell)
    np.testing.assert_array_equal(stack[..., 1], nucleus)
    np.testing.assert_array_equal(stack[..., 2], actin)


def test_build_membrane_stack_omits_actin_plane_when_actin_tif_is_none(tmp_path):
    cell = np.full((4, 4), 1, dtype=np.uint16)
    nucleus = np.full((4, 4), 2, dtype=np.uint16)

    cell_tif, nucleus_tif, _ = _write_tiles(tmp_path, cell=cell, nucleus=nucleus)

    stack = build_membrane_stack(cell_tif, nucleus_tif, None)

    assert stack.shape == (4, 4, 2)
    np.testing.assert_array_equal(stack[..., 0], cell)
    np.testing.assert_array_equal(stack[..., 1], nucleus)


def test_build_membrane_stack_rejects_mismatched_channel_shapes(tmp_path):
    cell = np.zeros((4, 4), dtype=np.uint16)
    nucleus = np.zeros((5, 5), dtype=np.uint16)

    cell_tif, nucleus_tif, _ = _write_tiles(tmp_path, cell=cell, nucleus=nucleus)

    with pytest.raises(ValueError, match="Channel shape mismatch"):
        build_membrane_stack(cell_tif, nucleus_tif, None)


def _write_tiles(tmp_path, *, cell=None, nucleus=None, actin=None):
    import tifffile

    def _write(name, arr):
        if arr is None:
            return None
        path = tmp_path / name
        tifffile.imwrite(path, arr)
        return path

    return _write("cell.tif", cell), _write("nucleus.tif", nucleus), _write("actin.tif", actin)


# --- membrane mode: reconciling composite channel count with the model -----


def test_reconcile_membrane_stack_noop_when_channel_counts_already_match():
    stack = np.zeros((4, 4, 2), dtype=np.uint16)
    out = reconcile_membrane_stack(stack, model_nchan=2)
    assert out is stack


def test_reconcile_membrane_stack_pads_zero_actin_for_a_3ch_model_on_a_2ch_run():
    stack = np.stack([np.full((2, 2), 1), np.full((2, 2), 2)], axis=-1).astype(np.uint16)
    out = reconcile_membrane_stack(stack, model_nchan=3)

    assert out.shape == (2, 2, 3)
    np.testing.assert_array_equal(out[..., 0], stack[..., 0])
    np.testing.assert_array_equal(out[..., 1], stack[..., 1])
    np.testing.assert_array_equal(out[..., 2], np.zeros((2, 2), dtype=np.uint16))


def test_reconcile_membrane_stack_drops_actin_for_a_2ch_model_on_a_3ch_run():
    stack = np.stack(
        [np.full((2, 2), 1), np.full((2, 2), 2), np.full((2, 2), 3)], axis=-1
    ).astype(np.uint16)
    out = reconcile_membrane_stack(stack, model_nchan=2)

    assert out.shape == (2, 2, 2)
    np.testing.assert_array_equal(out[..., 0], stack[..., 0])
    np.testing.assert_array_equal(out[..., 1], stack[..., 1])


def test_reconcile_membrane_stack_leaves_unreconcilable_mismatch_unchanged():
    # No pad/drop rule covers e.g. a 1-channel model on a 2-channel stack --
    # the function logs a warning and hands the stack back as-is rather than
    # guessing, so Cellpose itself raises on the mismatch.
    stack = np.zeros((2, 2, 2), dtype=np.uint16)
    out = reconcile_membrane_stack(stack, model_nchan=1)
    assert out.shape == (2, 2, 2)


# --- model_input_channels: reading nchan from a checkpoint's state dict ----


class _FakeTensor:
    """Minimal stand-in for a torch.Tensor: just .ndim and .shape."""

    def __init__(self, shape):
        self.shape = shape
        self.ndim = len(shape)


def _install_fake_torch(monkeypatch, load_result=None, raises=False):
    """
    Inject a fake ``torch`` module into sys.modules so
    ``model_input_channels()``'s internal ``import torch`` binds it instead
    of importing the real (heavyweight, GPU-oriented) package -- the same
    "imports live inside the function" isolation this test file already
    relies on for cellpose/torch generally.
    """

    def fake_load(path, map_location=None, weights_only=None):
        if raises:
            raise RuntimeError("cannot load")
        return load_result

    fake_torch = type(sys)("torch")
    fake_torch.load = fake_load
    monkeypatch.setitem(sys.modules, "torch", fake_torch)


def test_model_input_channels_finds_in_channels_from_first_4d_conv_weight(monkeypatch):
    # Deliberately does not hardcode a specific parameter name (see the
    # function's docstring) -- the 1D batchnorm weight ("...conv_0.0...")
    # is skipped because it isn't 4D, and the first 4D tensor's second
    # dimension (in_channels) is returned.
    state_dict = {
        "downsample.down.res_down_0.conv.conv_0.0.weight": _FakeTensor((32,)),
        "downsample.down.res_down_0.conv.conv_0.1.weight": _FakeTensor((32, 2, 3, 3)),
        "downsample.down.res_down_0.conv.conv_0.1.bias": _FakeTensor((32,)),
    }
    _install_fake_torch(monkeypatch, load_result=state_dict)

    assert model_input_channels(Path("fake_model.pt")) == 2


def test_model_input_channels_unwraps_model_state_dict_key(monkeypatch):
    checkpoint = {
        "model_state_dict": {
            "conv1.weight": _FakeTensor((16, 3, 3, 3)),
        },
        "optimizer_state_dict": {},
    }
    _install_fake_torch(monkeypatch, load_result=checkpoint)

    assert model_input_channels(Path("fake_model.pt")) == 3


def test_model_input_channels_returns_none_when_no_4d_tensor_present(monkeypatch):
    _install_fake_torch(monkeypatch, load_result={"some_scalar": _FakeTensor((1,))})
    assert model_input_channels(Path("fake_model.pt")) is None


def test_model_input_channels_returns_none_when_torch_load_fails(monkeypatch):
    _install_fake_torch(monkeypatch, raises=True)
    assert model_input_channels(Path("fake_model.pt")) is None


def test_model_input_channels_returns_none_for_a_non_dict_checkpoint(monkeypatch):
    _install_fake_torch(monkeypatch, load_result=_FakeTensor((4, 3, 3, 3)))
    assert model_input_channels(Path("fake_model.pt")) is None
