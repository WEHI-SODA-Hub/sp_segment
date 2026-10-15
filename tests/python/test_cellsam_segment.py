"""
Tests for bin/cellsam_segment.py that do not need cellSAM, torch or a GPU.

The cellSAM package is replaced by a stub whose pipeline returns a fixed mask,
so these cover the script's own logic: channel selection, the empty-mask
failure, and the border-cell flag.
"""

import sys
import types

import numpy as np
import pytest
import tifffile

# Stub cellSAM before importing the script.
_cellsam = types.ModuleType("cellSAM")
_cellsam.get_model = lambda *a, **k: None
_pipeline = types.ModuleType("cellSAM.cellsam_pipeline")
_pipeline.cellsam_pipeline = lambda img, **kwargs: _pipeline.next_mask
_cellsam.cellsam_pipeline = _pipeline
sys.modules.setdefault("cellSAM", _cellsam)
sys.modules.setdefault("cellSAM.cellsam_pipeline", _pipeline)

_torch = types.ModuleType("torch")
_torch.__version__ = "stub"
_torch.version = types.SimpleNamespace(cuda=None)
_torch.cuda = types.SimpleNamespace(is_available=lambda: False)
sys.modules.setdefault("torch", _torch)

import cellsam_segment  # noqa: E402


@pytest.fixture
def combined_tiff(tmp_path):
    path = tmp_path / "combined.ome.tif"
    planes = np.stack([np.full((20, 30), 100, np.uint16), np.full((20, 30), 200, np.uint16)])
    tifffile.imwrite(path, planes, ome=True,
                     metadata={"axes": "CYX", "Channel": {"Name": ["DAPI", "combined_membrane"]}})
    return path


def _run(monkeypatch, tiff, out, mask, *extra):
    _pipeline.next_mask = mask
    monkeypatch.delenv("DEEPCELL_ACCESS_TOKEN", raising=False)
    monkeypatch.setattr(sys, "argv", [
        "cellsam_segment.py", str(tiff), "--output", str(out),
        "--compartment", "whole-cell", "--nuclear-channel", "DAPI",
        "--membrane-channel", "combined_membrane", *extra,
    ])
    cellsam_segment.main()
    return tifffile.imread(out)


def _mask_with_edge_cell():
    mask = np.zeros((20, 30), np.uint32)
    mask[0:3, 0:3] = 1      # touches the edge
    mask[8:12, 10:14] = 2   # interior
    return mask


def test_extract_channels_orders_blank_nuclear_membrane(combined_tiff):
    img = cellsam_segment.extract_channels(combined_tiff, "DAPI", ["combined_membrane"], "whole-cell")
    assert img.shape == (20, 30, 3)
    assert img[0, 0].tolist() == [0, 100, 200]


def test_nuclear_compartment_leaves_membrane_blank(combined_tiff):
    img = cellsam_segment.extract_channels(combined_tiff, "DAPI", ["combined_membrane"], "nuclear")
    assert img[0, 0].tolist() == [0, 100, 0]


def test_border_cells_removed_by_default(monkeypatch, combined_tiff, tmp_path):
    out = _run(monkeypatch, combined_tiff, tmp_path / "m.tiff", _mask_with_edge_cell())
    assert set(np.unique(out)) == {0, 2}


def test_no_remove_border_cells_keeps_them(monkeypatch, combined_tiff, tmp_path):
    out = _run(monkeypatch, combined_tiff, tmp_path / "m.tiff", _mask_with_edge_cell(),
               "--no-remove-border-cells")
    assert set(np.unique(out)) == {0, 1, 2}


def test_empty_mask_fails(monkeypatch, combined_tiff, tmp_path):
    with pytest.raises(SystemExit) as exc:
        _run(monkeypatch, combined_tiff, tmp_path / "m.tiff", np.zeros((20, 30), np.uint32))
    assert "no cells" in str(exc.value)


def test_mask_empty_after_filtering_fails(monkeypatch, combined_tiff, tmp_path):
    mask = np.zeros((20, 30), np.uint32)
    mask[8:10, 10:12] = 1  # 4 px, below --min-area
    with pytest.raises(SystemExit) as exc:
        _run(monkeypatch, combined_tiff, tmp_path / "m.tiff", mask, "--min-area", "200")
    assert "no cells left" in str(exc.value)
