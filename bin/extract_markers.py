#!/usr/bin/env python
'''
Module      : extract_markers
Description : Extracts markers from COMET OME-TIFF and writes comma-delimited
              output containing markers, exposure times and background channels
              compatible with input expected by mcmicro's
              background_subtraction tool.
Copyright   : (c) WEHI SODA Hub, 2025
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import typer
from typing import Annotated, List
from tifffile import TiffFile


def _local_tag_name(tag: str) -> str:
    """XML tag name without its namespace."""
    return tag.split('}', 1)[-1]


def ome_channel_meta(ome_xml: str) -> List[dict]:
    """
    Marker name, exposure time and background channel for every OME Channel.

    Channels are numbered across all Image/Pixels elements in document order,
    so a file that stores one channel per Image (TheC is 0 in each) still gets
    each exposure on its own channel. Background channels are matched to
    channels by ID, not by position.
    """
    root = ET.fromstring(ome_xml)
    channel_meta: List[dict] = []
    by_id = {}

    for pixels in (e for e in root.iter() if _local_tag_name(e.tag) == 'Pixels'):
        offset = len(channel_meta)
        for ch in (e for e in pixels if _local_tag_name(e.tag) == 'Channel'):
            cm = {'index': ch.attrib.get('ID'), 'marker_name': ch.attrib.get('Name')}
            channel_meta.append(cm)
            if cm['index'] is not None:
                by_id[cm['index']] = cm
        n_here = len(channel_meta) - offset
        for plane in (e for e in pixels if _local_tag_name(e.tag) == 'Plane'):
            exposure_time = plane.attrib.get('ExposureTime')
            the_c = int(plane.attrib.get('TheC', 0))
            if exposure_time is not None and the_c < n_here:
                channel_meta[offset + the_c]['exposure'] = float(exposure_time)

    # COMET writes background channels as ChannelPriv elements (custom
    # annotations) whose ID is the Channel ID.
    for priv in (e for e in root.iter() if _local_tag_name(e.tag) == 'ChannelPriv'):
        cm = by_id.get(priv.attrib.get('ID'))
        if cm is None:
            continue
        background = priv.attrib.get('FluorescenceChannel')
        cm['background'] = "" if background == cm['marker_name'] else background

    return channel_meta



def main(
    tiff_path: Annotated[Path, typer.Argument(
        help="Path to the TIFF input file."
    )],
    remove_marker: Annotated[List[str], typer.Option(
        '--remove-marker',
        '-r',
        help="List of markers to remove from the output (case sensitive)."
    )] = []
):
    """
    Extracts markers, exposure times, and background channels from COMET
    OME-TIFF
    """
    with TiffFile(tiff_path) as tiff:
        ome_metadata = tiff.ome_metadata
        imagej_metadata = getattr(tiff, 'imagej_metadata', None)

    channel_meta = []

    if ome_metadata:
        try:
            channel_meta = ome_channel_meta(ome_metadata)
        except ET.ParseError:
            channel_meta = []

    # Fallback for ImageJ metadata exports (no OME-XML)
    if not channel_meta and isinstance(imagej_metadata, dict):
        labels = imagej_metadata.get('Labels') or []
        if labels:
            channel_meta = [
                {
                    'index': str(i),
                    'marker_name': str(name),
                    'background': '',
                    'exposure': 1.0,
                }
                for i, name in enumerate(labels)
            ]
            typer.echo(
                "Warning: OME-XML metadata not found/invalid. "
                "Using ImageJ channel labels with default exposure=1.0 and blank background.",
                err=True
            )

    if not channel_meta:
        raise ValueError(
            "Could not extract channel metadata from OME-XML or ImageJ labels."
        )

    # Validate remove channel list
    markers: List[str] = [cm['marker_name'] for cm in channel_meta]
    for remove in remove_marker:
        if remove not in markers:
            raise ValueError(
                f"Channel '{remove}' not found in the OME-TIFF metadata."
            )

    # Ensure required fields exist in fallback mode
    for cm in channel_meta:
        if 'background' not in cm or cm['background'] is None:
            cm['background'] = ''
        if 'exposure' not in cm or cm['exposure'] is None:
            cm['exposure'] = 1.0

    # Format output and write to stdout
    df = pd.DataFrame(channel_meta)
    df = df[['marker_name', 'background', 'exposure']]
    df['exposure'] = df['exposure'].astype(float)

    # Handle removal of channels
    markers_to_remove: List[str] = \
        ['TRUE' if marker in remove_marker else ''
            for marker in df['marker_name']]
    df['remove'] = markers_to_remove

    df.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    typer.run(main)
