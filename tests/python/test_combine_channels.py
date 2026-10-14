"""
Tests for bin/combine_channels.py and bin/tiff_channels.py.

Covers the product scaling, both channel layouts (one CYX series, and one 2D
series per channel as in the COMET _fixed export), MIBI detection, and the
OME-XML written for the 2-channel output.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import tifffile

from combine_channels import CombineMethod, combine, to_uint16
from tiff_channels import channel_layout, get_channel_names, read_channels

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "combine_channels.py"
NAMES = ["DAPI", "CD45", "HLA I", "Pan CK"]


def _planes(n=len(NAMES), shape=(32, 48)):
    # Channel c holds the constant value (c + 1) * 1000, so reads are checkable.
    return np.stack([np.full(shape, (c + 1) * 1000, dtype=np.uint16) for c in range(n)])


@pytest.fixture
def stacked_tiff(tmp_path):
    path = tmp_path / "stacked.ome.tif"
    tifffile.imwrite(
        path, _planes(), ome=True,
        metadata={"axes": "CYX", "Channel": {"Name": NAMES},
                  "PhysicalSizeX": 0.28, "PhysicalSizeY": 0.28},
    )
    return path


@pytest.fixture
def split_tiff(tmp_path):
    # One OME Image (series) per channel, each with a single named channel.
    path = tmp_path / "split.ome.tif"
    planes = _planes()
    with tifffile.TiffWriter(path, ome=True) as tif:
        for name, plane in zip(NAMES, planes):
            tif.write(plane, metadata={"axes": "YX", "Channel": {"Name": name}})
    return path


@pytest.fixture
def mibi_tiff(tmp_path):
    path = tmp_path / "mibi.tif"
    with tifffile.TiffWriter(path) as tif:
        for name, plane in zip(NAMES, _planes()):
            tif.write(plane, description=json.dumps({"channel.target": name}), metadata=None)
    return path


# ---------------------------------------------------------------------------
# combine()
# ---------------------------------------------------------------------------


def test_prod_rescales_instead_of_saturating():
    planes = np.stack([
        np.array([[1000, 10]], dtype=np.uint16),
        np.array([[2000, 10]], dtype=np.uint16),
    ])
    out = combine(planes, CombineMethod.PROD)
    # 2,000,000 is the peak and maps to 65535; 100 scales by the same factor.
    assert out.dtype == np.uint16
    assert out[0, 0] == 65535
    assert out[0, 1] == round(100 * 65535 / 2_000_000)


def test_prod_does_not_wrap_with_many_channels():
    # 60000**5 overflows uint64; the float64 product must stay ordered.
    planes = np.stack([np.array([[60000, 30000]], dtype=np.uint16)] * 5)
    out = combine(planes, CombineMethod.PROD)
    assert out[0, 0] == 65535
    assert 0 < out[0, 1] < out[0, 0]


def test_prod_within_range_is_unchanged():
    planes = np.stack([np.array([[10, 20]], dtype=np.uint16)] * 2)
    assert combine(planes, CombineMethod.PROD).tolist() == [[100, 400]]


def test_max():
    planes = np.stack([np.array([[1, 9]], dtype=np.uint16), np.array([[5, 2]], dtype=np.uint16)])
    assert combine(planes, CombineMethod.MAX).tolist() == [[5, 9]]


def test_to_uint16_clips_negative_floats():
    out = to_uint16(np.array([-3.0, 1.6, 70000.0], dtype=np.float32))
    assert out.tolist() == [0, 2, 65535]


# ---------------------------------------------------------------------------
# tiff_channels
# ---------------------------------------------------------------------------


def test_stacked_layout_names_and_reads(stacked_tiff):
    with tifffile.TiffFile(stacked_tiff) as tif:
        assert channel_layout(tif) == ("stacked", 4)
    assert get_channel_names(stacked_tiff) == NAMES
    planes = read_channels(stacked_tiff, [2, 0])
    assert planes.shape == (2, 32, 48)
    assert planes[0, 0, 0] == 3000 and planes[1, 0, 0] == 1000


def test_split_layout_names_and_reads(split_tiff):
    with tifffile.TiffFile(split_tiff) as tif:
        assert channel_layout(tif) == ("split", 4)
    assert get_channel_names(split_tiff) == NAMES
    planes = read_channels(split_tiff, [3, 1])
    assert planes[0, 0, 0] == 4000 and planes[1, 0, 0] == 2000


def test_mibi_names(mibi_tiff):
    assert get_channel_names(mibi_tiff) == NAMES


def test_mibi_series_with_singleton_axis():
    # The repository's MIBI test image reads as one (1, Y, X) series per page.
    path = Path(__file__).resolve().parents[1] / "data" / "mesmer" / "test_data.tiff"
    with tifffile.TiffFile(path) as tif:
        assert [s.shape for s in tif.series] == [(1, 512, 512), (1, 512, 512)]
        assert channel_layout(tif) == ("split", 2)
    assert get_channel_names(path) == ["nuclear", "membrane"]
    assert read_channels(path, [1, 0]).shape == (2, 512, 512)


def test_plain_json_description_is_not_mibi(tmp_path):
    # tifffile's default description is {"shape": [...]}, which parses as JSON.
    path = tmp_path / "plain.tif"
    tifffile.imwrite(path, _planes(), photometric="minisblack")
    assert get_channel_names(path) == [f"Channel_{i}" for i in range(4)]


def test_interleaved_layout(tmp_path):
    path = tmp_path / "yxc.tif"
    tifffile.imwrite(path, np.moveaxis(_planes(n=3), 0, -1), photometric="rgb")
    with tifffile.TiffFile(path) as tif:
        assert channel_layout(tif) == ("stacked", 3)
    assert read_channels(path, [1])[0, 0, 0] == 2000


# ---------------------------------------------------------------------------
# combine_channels.py end to end
# ---------------------------------------------------------------------------


def _run(tiff, out, *args):
    with open(out, "wb") as fh:
        return subprocess.run(
            [sys.executable, str(SCRIPT), str(tiff), *args],
            stdout=fh, stderr=subprocess.PIPE, text=True,
        )


@pytest.mark.parametrize("fixture", ["stacked_tiff", "split_tiff", "mibi_tiff"])
def test_cli_combines_and_writes_valid_ome(fixture, request, tmp_path):
    src = request.getfixturevalue(fixture)
    out = tmp_path / "combined.tiff"
    result = _run(src, out, "--nuclear-channel", "DAPI",
                  "--membrane-channel", "CD45", "--membrane-channel", "HLA I",
                  "--pixel-size-microns", "0.39")
    assert result.returncode == 0, result.stderr

    assert get_channel_names(out) == ["DAPI", "combined_membrane"]
    data = tifffile.imread(out)
    assert data.shape == (2, 32, 48)
    assert data[0, 0, 0] == 1000
    assert data[1, 0, 0] == 3000  # max(2000, 3000)

    with tifffile.TiffFile(out) as tif:
        xml = tif.ome_metadata
    assert 'SizeC="2"' in xml
    # The stacked source carries 0.28 um; the others fall back to the option.
    expected = "0.28" if fixture == "stacked_tiff" else "0.39"
    assert f'PhysicalSizeX="{expected}"' in xml


def test_cli_single_membrane_keeps_its_name(stacked_tiff, tmp_path):
    out = tmp_path / "combined.tiff"
    result = _run(stacked_tiff, out, "--nuclear-channel", "DAPI", "--membrane-channel", "Pan CK")
    assert result.returncode == 0, result.stderr
    assert get_channel_names(out) == ["DAPI", "Pan CK"]


def test_cli_unknown_channel_lists_available(stacked_tiff, tmp_path):
    result = _run(stacked_tiff, tmp_path / "x.tiff",
                  "--nuclear-channel", "DAPI", "--membrane-channel", "CD99")
    assert result.returncode != 0
    assert "CD99" in result.stderr and "Pan CK" in result.stderr
