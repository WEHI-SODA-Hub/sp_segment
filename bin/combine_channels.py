#!/usr/bin/env python
'''
Module      : combine_channels
Description : Takes a multichannel TIFF (OME-TIFF including OPAL OME/QPTIFF,
              MIBI TIFF, or ImageJ TIFF) and writes a 2-channel OME-TIFF to
              stdout: the nuclear channel, and a membrane channel that is
              either one input channel or a combination of several (the
              per-pixel maximum, or the mean of percentile-rescaled channels).
Copyright   : (c) WEHI SODA Hub, 2025
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import sys
import xml.etree.ElementTree as ET
from enum import Enum
from pathlib import Path
from typing import Annotated, List, Optional

import numpy as np
import typer
from tifffile import TiffFile, imwrite

from tiff_channels import (
    get_channel_names,
    ome_physical_size,
    read_channels,
)

UINT16_MAX = np.iinfo(np.uint16).max
COMBINED_MEMBRANE = "combined_membrane"


class CombineMethod(str, Enum):
    MAX = "max"
    NORM_MEAN = "norm_mean"


# norm_mean rescales each channel so this percentile maps to 0 and
# NORM_HIGH_PERCENTILE maps to 1, before averaging.
NORM_LOW_PERCENTILE = 50
NORM_HIGH_PERCENTILE = 99.5


def to_uint16(data: np.ndarray) -> np.ndarray:
    '''
    Cast to uint16 without wrapping or losing float resolution.

    Integer data is clipped to 0..65535. Float data has negative values set
    to 0 (e.g. from background subtraction), is scaled so its maximum is
    65535, then rounded. Rounding the raw values would leave few levels for
    small-valued float images (OPAL component data with values 0-5 would
    keep 6). The segmenters normalise each channel by its own range or
    percentiles, so only the resolution matters, not the scale.
    '''
    if data.dtype == np.uint16:
        return data
    if np.issubdtype(data.dtype, np.floating):
        data = np.clip(data, 0, None).astype(np.float32)
        peak = float(data.max())
        if peak > 0:
            data *= np.float32(UINT16_MAX / peak)
        return np.rint(data).astype(np.uint16)
    return np.clip(data, 0, UINT16_MAX).astype(np.uint16)


def combine(planes: np.ndarray, method: CombineMethod) -> np.ndarray:
    '''
    Combine a (C, Y, X) stack into one uint16 (Y, X) plane.

    max takes the per-pixel maximum and keeps the input scale.

    norm_mean rescales each channel so that its 50th percentile is 0 and its
    99.5th percentile is 1 (values below 0 set to 0, values above 1 kept),
    takes the per-pixel mean, and maps the result linearly so its maximum is
    65535. Each channel then contributes on the same scale, so a dim marker is
    not hidden by a bright one. A channel whose two percentiles are equal
    contributes 0.
    '''
    if method == CombineMethod.MAX:
        return to_uint16(np.max(planes, axis=0))
    if method == CombineMethod.NORM_MEAN:
        total = np.zeros(planes.shape[1:], dtype=np.float32)
        for plane in planes:
            low, high = np.percentile(plane, [NORM_LOW_PERCENTILE, NORM_HIGH_PERCENTILE])
            if high > low:
                scaled = plane.astype(np.float32)
                scaled -= np.float32(low)
                scaled *= np.float32(1.0 / (high - low))
                np.maximum(scaled, 0, out=scaled)
                total += scaled
        peak = float(total.max())
        if peak > 0:
            total *= np.float32(UINT16_MAX / peak)
        return to_uint16(total)
    raise ValueError(f"unknown combine method: {method}")


def membrane_copy_name(nuclear_channel: str) -> str:
    '''Output name of the membrane channel when it is the nuclear channel.'''
    return f"{nuclear_channel}_membrane"


def channel_index(name: str, channel_names: List[str]) -> int:
    '''Index of a channel name; fails with the list of available names.'''
    try:
        return channel_names.index(name)
    except ValueError:
        raise typer.BadParameter(
            f"channel '{name}' is not in the image. Available channels: {channel_names}"
        )


def build_ome_xml(width: int, height: int, channel_names: List[str],
                  physical_size: Optional[tuple]) -> str:
    '''
    OME-XML for the 2-channel output, built from scratch.

    The source XML is not edited in place: it can use any OME schema version,
    and describes N channels, per-plane TiffData and Plane elements that do
    not apply to the output. Only the pixel size is carried over.
    '''
    ome_ns = "http://www.openmicroscopy.org/Schemas/OME/2016-06"
    xsi_ns = "http://www.w3.org/2001/XMLSchema-instance"
    schema_loc = f"{ome_ns} {ome_ns}/ome.xsd"
    ET.register_namespace("", ome_ns)
    ET.register_namespace("xsi", xsi_ns)

    # Clark-notation tag prefixes, built without {{{ patterns that trip nf-core linting
    tag_ome = "{" + ome_ns + "}"
    tag_xsi = "{" + xsi_ns + "}"

    root = ET.Element(tag_ome + "OME")
    root.set(tag_xsi + "schemaLocation", schema_loc)

    image = ET.SubElement(root, tag_ome + "Image", {"ID": "Image:0", "Name": "combined"})
    pixels_attrs = {
        "ID": "Pixels:0",
        "DimensionOrder": "XYZCT",
        "Type": "uint16",
        "SizeX": str(width),
        "SizeY": str(height),
        "SizeZ": "1",
        "SizeC": str(len(channel_names)),
        "SizeT": "1",
    }
    if physical_size is not None:
        size_x, size_y, unit = physical_size
        pixels_attrs.update({
            "PhysicalSizeX": str(size_x),
            "PhysicalSizeY": str(size_y),
            "PhysicalSizeXUnit": unit,
            "PhysicalSizeYUnit": unit,
        })
    pixels = ET.SubElement(image, tag_ome + "Pixels", pixels_attrs)
    # Schema order inside Pixels: Channel elements first, then TiffData.
    for c, name in enumerate(channel_names):
        ET.SubElement(pixels, tag_ome + "Channel", {
            "ID": f"Channel:0:{c}",
            "Name": name,
            "SamplesPerPixel": "1",
        })
    for c in range(len(channel_names)):
        ET.SubElement(pixels, tag_ome + "TiffData", {
            "FirstC": str(c), "FirstT": "0", "FirstZ": "0",
            "IFD": str(c), "PlaneCount": "1",
        })

    return '<?xml version="1.0" encoding="UTF-8"?>\n' + \
        ET.tostring(root, encoding='unicode')


def source_physical_size(tiff: Path) -> Optional[tuple]:
    '''Pixel size from the source OME-XML, or None when it has none.'''
    with TiffFile(tiff) as tif:
        xml = tif.ome_metadata
        if not xml and tif.pages:
            desc = tif.pages[0].description
            xml = desc if isinstance(desc, str) and desc.lstrip().startswith("<") else None
    if not xml:
        return None
    try:
        return ome_physical_size(xml)
    except (ET.ParseError, ValueError):
        return None


def main(
    tiff: Annotated[Path, typer.Argument(
        help="Path to the TIFF input file."
    )],
    nuclear_channel: Annotated[str, typer.Option(
        help="Name of the nuclear channel."
    )],
    membrane_channel: Annotated[List[str], typer.Option(
        help="Name of a membrane channel. Repeat for several channels.")
    ],
    combine_method: Annotated[CombineMethod, typer.Option(
        help="Method used to combine several membrane channels: max or norm_mean.")
    ] = CombineMethod.MAX,
    pixel_size_microns: Annotated[Optional[float], typer.Option(
        help="Pixel size in microns. Written to OME-XML PhysicalSizeX/Y when "
             "the source image has no calibration metadata.")
    ] = None,
):
    channel_names = get_channel_names(tiff)
    nuclear_index = channel_index(nuclear_channel, channel_names)
    membrane_indices = [channel_index(name, channel_names) for name in membrane_channel]

    # Only the nuclear and membrane channels are decoded, not the whole image.
    nuclear = to_uint16(read_channels(tiff, [nuclear_index])[0])
    if len(membrane_indices) > 1:
        membrane = combine(read_channels(tiff, membrane_indices), CombineMethod(combine_method))
        membrane_name = COMBINED_MEMBRANE
    else:
        membrane = to_uint16(read_channels(tiff, membrane_indices)[0])
        membrane_name = membrane_channel[0]
        # The nuclear channel can also be the membrane channel (nuclear-only
        # input). Channel names must be unique for sopa, so rename the copy.
        if membrane_name == nuclear_channel:
            membrane_name = membrane_copy_name(nuclear_channel)

    # Report what was combined so it stays visible in the process log
    # (downstream segmenters only see the single combined membrane channel).
    typer.echo(
        f"[combine_channels] nuclear='{nuclear_channel}'; "
        f"membrane={membrane_channel} -> '{membrane_name}' "
        f"(method={CombineMethod(combine_method).value}, "
        f"n={len(membrane_channel)})",
        err=True,
    )

    output = np.stack([nuclear, membrane])
    _, height, width = output.shape

    physical_size = source_physical_size(tiff)
    if physical_size is None and pixel_size_microns is not None:
        physical_size = (pixel_size_microns, pixel_size_microns, "µm")

    ome_xml = build_ome_xml(width, height, [nuclear_channel, membrane_name], physical_size)

    imwrite(sys.stdout.buffer, output,
            photometric='minisblack',
            metadata=None,
            description=ome_xml.encode('utf-8'))


if __name__ == "__main__":
    typer.run(main)
