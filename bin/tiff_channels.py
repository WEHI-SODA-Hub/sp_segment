#!/usr/bin/env python3
"""
Channel layout, channel names and channel reads for multichannel TIFFs.

Shared by the bin/ scripts so they agree on what a channel is. Needs only the
standard library, numpy and tifffile, so it can be imported from any of the
pipeline's containers.

Two layouts occur in practice. The usual one keeps every channel in a single
``CYX`` series. The other gives each channel its own 2D series, which some OME
writers produce -- notably the inForm unmixing export behind the COMET
``_fixed`` images. A reader that only consults ``series[0]`` or the first OME
``Image`` sees the second as a one-channel image.
"""

import json
import xml.etree.ElementTree as ET

import numpy as np
import tifffile


def _local_tag_name(tag):
    """XML tag name without its namespace."""
    return tag.split("}", 1)[-1]


def _plane_shape(series):
    """Series shape without singleton axes; MIBI pages read as ``(1, Y, X)``."""
    return tuple(d for d in series.shape if d != 1)


def channel_layout(tif):
    """How a TIFF stores its channels: ``("stacked" | "split", n_channels)``."""
    series = list(getattr(tif, "series", None) or [])
    if (
        len(series) > 1
        and all(len(_plane_shape(s)) == 2 for s in series)
        and len({_plane_shape(s) for s in series}) == 1
    ):
        return "split", len(series)
    if series:
        shape = tuple(series[0].shape)
        if len(shape) >= 3:
            # Interleaved single-page images store channels last (YXC/YXS).
            if is_channel_last(tif):
                return "stacked", shape[-1]
            return "stacked", shape[0]
        return "stacked", 1
    return "stacked", (len(tif.pages) if len(tif.pages) > 1 else 1)


def is_channel_last(tif):
    """True for a single interleaved series stored as ``YXC``/``YXS``."""
    series = list(getattr(tif, "series", None) or [])
    if not series or len(series[0].shape) < 3:
        return False
    axes = series[0].axes
    return axes[-1] in "CS" and axes[0] not in "CIQ"


def _mibi_channel_names(tif):
    """``channel.target`` from each page's JSON description, or ``[]``."""
    names = []
    for page in tif.pages:
        try:
            desc = json.loads(page.description)
        except (json.JSONDecodeError, TypeError):
            return []
        # tifffile's own {"shape": [...]} description also parses as JSON, so
        # require the MIBI key rather than treating any JSON as MIBI.
        if not isinstance(desc, dict) or "channel.target" not in desc:
            return []
        names.append(desc["channel.target"])
    return names


def _ome_xml(tif):
    """The OME-XML string, also when it is only in the first page description."""
    xml = tif.ome_metadata
    if not xml and tif.pages:
        desc = tif.pages[0].description
        if isinstance(desc, str) and desc.lstrip().startswith("<"):
            xml = desc
    return xml


def ome_channel_names(xml):
    """Every OME ``Channel`` name, in document order, in any OME namespace.

    Taking channels from every ``Image`` covers the split layout, where each
    ``Image`` holds one channel.
    """
    root = ET.fromstring(xml)
    return [
        elem.get("Name") or elem.get("ID")
        for elem in root.iter()
        if _local_tag_name(elem.tag) == "Channel"
    ]


def ome_physical_size(xml):
    """``(PhysicalSizeX, PhysicalSizeY, unit)`` of the first ``Pixels``, or ``None``."""
    root = ET.fromstring(xml)
    for elem in root.iter():
        if _local_tag_name(elem.tag) == "Pixels":
            x, y = elem.get("PhysicalSizeX"), elem.get("PhysicalSizeY")
            if x is None or y is None:
                return None
            return float(x), float(y), elem.get("PhysicalSizeXUnit", "µm")
    return None


def get_channel_names(tiff_path):
    """Channel names from MIBI JSON, OME-XML or ImageJ labels, else numbered.

    Reads metadata only. The result has exactly one entry per image channel: a
    channel the metadata leaves unnamed gets a synthesised name rather than
    being dropped, because dropping would slide every later name onto the
    wrong channel index.
    """
    with tifffile.TiffFile(tiff_path) as tif:
        _, n_channels = channel_layout(tif)

        names = _mibi_channel_names(tif)
        if not names:
            xml = _ome_xml(tif)
            if xml:
                try:
                    names = ome_channel_names(xml)
                except ET.ParseError:
                    names = []
        if not names and getattr(tif, "imagej_metadata", None):
            names = [str(label) for label in tif.imagej_metadata.get("Labels", [])]

    names = [name or f"Channel_{i}" for i, name in enumerate(names[:n_channels])]
    names += [f"Channel_{i}" for i in range(len(names), n_channels)]
    return names


def read_channels(tiff_path, indices):
    """Full-resolution planes for the given channel indices: ``(len(indices), Y, X)``.

    Only the requested channels are decoded where the layout allows it (one
    page or one series per channel); an interleaved single-page image has to
    be decoded whole.
    """
    indices = list(indices)
    with tifffile.TiffFile(tiff_path) as tif:
        layout, n_channels = channel_layout(tif)
        bad = [i for i in indices if not 0 <= i < n_channels]
        if bad:
            raise IndexError(f"channel indices {bad} out of range for {n_channels} channels")

        if layout == "split":
            planes = [tif.series[i].asarray() for i in indices]
            return np.stack([p.reshape(p.shape[-2:]) for p in planes])

        series = tif.series[0]
        pages = series.pages
        if len(pages) == n_channels and all(p is not None and len(p.shape) == 2 for p in pages):
            return np.stack([pages[i].asarray() for i in indices])

        data = series.asarray()
        if data.ndim == 2:
            data = data[np.newaxis]
        elif data.ndim == 3 and data.shape[0] != n_channels and data.shape[-1] == n_channels:
            data = np.moveaxis(data, -1, 0)
        data = data.reshape(-1, data.shape[-2], data.shape[-1])
        return data[indices]
