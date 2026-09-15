#!/usr/bin/env python
'''
Module      : aviti_merge_plate_geojson
Description : Merges each well's cellmeasurement GeoJSON into one plate-level
              FeatureCollection in plate pixel coordinates, using the same
              well placements aviti_assemble_plate_image.py wrote to the
              plate layout CSV -- so cell polygons and the plate image agree
              on where every well is.

              A well's GeoJSON is never fully parsed into memory: a single
              AVITI well can carry on the order of 1 GB of feature text
              (10-20 GB once parsed into Python objects), and this container
              has no streaming JSON library, so features are extracted one at
              a time by a small string/escape-aware brace scanner
              (iter_geojson_features) and re-serialised straight to the
              output as each is transformed.

              cellmeasurement's coordinates are always full-resolution,
              mask-frame pixel coordinates (see cellmeasurement-py's
              geojson_writer: image_shape comes from the label mask, not from
              --downsample-factor), so merging needs only a rigid per-well
              translation -- never a rescale.
Copyright   : (c) WEHI SODA Hub, 2026
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import csv
import gzip
import json
import sys
from pathlib import Path
from typing import Annotated, Dict, Iterable, Iterator, List, TextIO

import typer

from aviti_plate_common import PlacedWell, parse_well_id, read_layout_csv

app = typer.Typer(add_completion=False)

# cellmeasurement's hardcoded bounds annotation (geojson_writer.py); dropped by
# default since the plate has its own well-scoped rectangles instead.
WHOLE_IMAGE_ANNOTATION_ID = "annotation-whole-image"

# properties.* integer namespaces that collide across wells and must be
# offset independently -- they are NOT the same counter (a cell's nucleus and
# whole-cell mask labels are drawn from separate label spaces).
ID_NAMESPACE_KEYS = ("id", "nucleus_label", "whole_cell_label")

_READ_CHUNK = 1 << 20  # 1 MiB


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def _open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def _find_features_array_start(fh: TextIO) -> str:
    '''
    Read forward from the start of a GeoJSON FeatureCollection document until
    just past the '[' opening its top-level "features" array, and return
    whatever was read past that point.

    A plain substring search over the accumulated prefix is sufficient here:
    the document's preamble (its "type" key) is a few dozen bytes, so this
    never has to scan into the -- potentially gigabyte-scale -- features
    array itself to find where that array begins.
    '''
    prefix = ""
    while True:
        idx = prefix.find('"features"')
        if idx != -1:
            after_key = prefix[idx + len('"features"'):]
            colon = after_key.find(":")
            if colon != -1:
                bracket = after_key.find("[", colon + 1)
                if bracket != -1:
                    return after_key[bracket + 1:]
        chunk = fh.read(_READ_CHUNK)
        if not chunk:
            raise ValueError('Could not find a top-level "features" array in the input.')
        prefix += chunk


def _iter_top_level_json_values(fh: TextIO, initial_buf: str) -> Iterator[str]:
    '''
    Yield the raw text of each comma-separated JSON value in an array that is
    already open (its '[' has been consumed), reading further chunks from
    ``fh`` only as needed and string/escape aware so braces or brackets
    inside a string property never end a value early.

    Memory use is bounded by roughly one feature's text plus at most one
    extra read-ahead chunk, never by the whole array.
    '''
    buf = initial_buf

    def fill(min_len: int) -> bool:
        nonlocal buf
        while len(buf) < min_len:
            chunk = fh.read(_READ_CHUNK)
            if not chunk:
                return False
            buf += chunk
        return True

    def skip_ws() -> None:
        nonlocal buf
        pos = 0
        while True:
            while pos < len(buf) and buf[pos].isspace():
                pos += 1
            if pos < len(buf):
                buf = buf[pos:]
                return
            if not fill(pos + 1):
                buf = ""
                return

    while True:
        skip_ws()
        if not buf and not fill(1):
            raise ValueError("Unexpected end of input inside a JSON array")
        if buf[0] == "]":
            buf = buf[1:]
            return
        if buf[0] != "{":
            raise ValueError(f"Expected a feature object, found {buf[0]!r}")

        depth = 0
        in_string = False
        escape = False
        pos = 0
        while True:
            if pos >= len(buf):
                if not fill(pos + 1):
                    raise ValueError("Unexpected end of input inside a feature object")
                continue
            ch = buf[pos]
            if in_string:
                if escape:
                    escape = False
                elif ch == "\\":
                    escape = True
                elif ch == '"':
                    in_string = False
            elif ch == '"':
                in_string = True
            elif ch in "{[":
                depth += 1
            elif ch in "}]":
                depth -= 1
            pos += 1
            if not in_string and depth == 0:
                break

        feature_text = buf[:pos]
        buf = buf[pos:]
        yield feature_text

        skip_ws()
        if buf and buf[0] == ",":
            buf = buf[1:]
        elif not buf:
            fill(1)


def iter_geojson_features(fh: TextIO) -> Iterator[str]:
    '''
    Yield the raw JSON text of each element of a GeoJSON FeatureCollection's
    top-level "features" array, one at a time, without ever holding the whole
    document in memory. See the module docstring for why this exists instead
    of ``json.load``.
    '''
    leftover = _find_features_array_start(fh)
    yield from _iter_top_level_json_values(fh, leftover)


def translate_coordinates(coords, dx: float, dy: float):
    '''
    Recursively translate a GeoJSON ``coordinates`` value by (dx, dy).
    Handles Point, LineString/ring, Polygon (with holes) and MultiPolygon
    uniformly: the base case is detected as the innermost [x, y] pair, so no
    geometry-type branching is needed.
    '''
    if (
        len(coords) >= 2
        and isinstance(coords[0], (int, float))
        and isinstance(coords[1], (int, float))
    ):
        return [coords[0] + dx, coords[1] + dy]
    return [translate_coordinates(c, dx, dy) for c in coords]


def transform_feature(
    feature: dict, well: str, dx: float, dy: float, id_offsets: Dict[str, int],
) -> dict:
    '''
    Return a copy of ``feature`` translated into plate coordinates, tagged
    with its well, and with its colliding integer id namespaces
    (properties.id, nucleus_label, whole_cell_label) offset -- each
    independently, since they are separate label spaces, not one counter.
    The string feature id ("cell-N") is made unique with a well prefix
    instead of arithmetic, since it isn't purely numeric.
    '''
    out = dict(feature)

    geometry = out.get("geometry")
    if geometry:
        out["geometry"] = {**geometry, "coordinates": translate_coordinates(geometry["coordinates"], dx, dy)}

    nucleus_geometry = out.get("nucleusGeometry")
    if nucleus_geometry:
        out["nucleusGeometry"] = {
            **nucleus_geometry,
            "coordinates": translate_coordinates(nucleus_geometry["coordinates"], dx, dy),
        }

    properties = dict(out.get("properties") or {})
    properties["well"] = well
    for key in ID_NAMESPACE_KEYS:
        if properties.get(key) is not None:
            properties[key] = properties[key] + id_offsets.get(key, 0)
    out["properties"] = properties

    if out.get("id") is not None:
        out["id"] = f"{well}-{out['id']}"

    return out


def well_annotation_feature(placed: PlacedWell) -> dict:
    '''
    A locked, axis-aligned rectangle at this well's exact plate placement,
    labelled with the well id, so QuPath shows where every well sits even
    before any cells are loaded.
    '''
    x0, y0 = placed.x0, placed.y0
    x1, y1 = x0 + placed.width, y0 + placed.height
    return {
        "type": "Feature",
        "id": f"well-{placed.well}",
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]],
        },
        "properties": {
            "objectType": "annotation",
            "type": "annotation",
            "name": placed.well,
            "well": placed.well,
            "isLocked": True,
        },
    }


def read_well_rows(manifest: Path) -> List[dict]:
    with open(manifest, newline="") as fh:
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in csv.DictReader(fh)
        ]
    if not rows:
        raise ValueError(f"No rows found in geojson manifest {manifest}")
    return sorted(rows, key=lambda r: parse_well_id(r["well"]))


def merge(
    rows: Iterable[dict],
    placements_by_well: Dict[str, PlacedWell],
    output: TextIO,
    *,
    keep_whole_image_annotations: bool = False,
    classify_by_well: bool = False,
) -> int:
    '''
    Stream the merged plate FeatureCollection to ``output``. Returns the
    number of cell features written (well rectangles are not counted).

    Well rectangles are written first, so they are present even if a later
    well's file is truncated or missing. Each well's running id offsets are
    seeded from the previous well's own maxima (the same running-max idiom
    ``stitch_masks()`` uses for mask labels), tracked per namespace so
    ``nucleus_label`` and ``whole_cell_label`` never bleed into each other.
    '''
    rows = list(rows)
    output.write('{"type":"FeatureCollection","features":[')

    first = True

    def write_feature(feature: dict) -> None:
        nonlocal first
        if not first:
            output.write(",")
        output.write(json.dumps(feature))
        first = False

    for row in rows:
        write_feature(well_annotation_feature(placements_by_well[row["well"]]))

    n_cells = 0
    running_max = {key: 0 for key in ID_NAMESPACE_KEYS}
    for row in rows:
        placed = placements_by_well[row["well"]]
        well_local_max = dict(running_max)
        with _open_maybe_gzip(row["geojson"]) as fh:
            for feature_text in iter_geojson_features(fh):
                feature = json.loads(feature_text)
                if feature.get("id") == WHOLE_IMAGE_ANNOTATION_ID and not keep_whole_image_annotations:
                    continue

                transformed = transform_feature(feature, placed.well, placed.x0, placed.y0, running_max)
                if classify_by_well:
                    transformed["properties"]["classification"] = {"name": placed.well}

                for key in ID_NAMESPACE_KEYS:
                    value = transformed["properties"].get(key)
                    if value is not None:
                        well_local_max[key] = max(well_local_max[key], value)

                write_feature(transformed)
                if transformed["properties"].get("objectType") != "annotation":
                    n_cells += 1
        running_max = well_local_max

    output.write("]}")
    return n_cells


@app.command()
def main(
    manifest: Annotated[Path, typer.Argument(
        exists=True, help="CSV manifest of per-well GeoJSON files (columns: well, geojson)."
    )],
    layout: Annotated[Path, typer.Option(
        exists=True, help="Plate layout CSV written by aviti_assemble_plate_image.py."
    )],
    output: Annotated[Path, typer.Option(
        help="Path to write the merged GeoJSON. A '.gz' suffix is appended when --gzip is set."
    )],
    gzip_output: Annotated[bool, typer.Option(
        "--gzip/--no-gzip", help="Gzip the output, matching cellmeasurement's own convention."
    )] = False,
    keep_whole_image_annotations: Annotated[bool, typer.Option(
        help="Keep each well's per-well 'whole_image' bounds annotation instead of dropping it "
             "in favour of the plate's own well rectangles."
    )] = False,
    classify_by_well: Annotated[bool, typer.Option(
        help="Set properties.classification.name to the well id on every cell, so QuPath can "
             "colour objects by well."
    )] = False,
):
    '''
    Merge per-well GeoJSON annotations into one plate-level FeatureCollection.
    '''
    rows = read_well_rows(manifest)
    placements = {p.well: p for p in read_layout_csv(layout)}

    missing = [row["well"] for row in rows if row["well"] not in placements]
    if missing:
        raise ValueError(f"Well(s) {missing} have no entry in the plate layout {layout}.")

    output_path = Path(f"{output}.gz") if gzip_output else output
    open_fn = gzip.open if gzip_output else open
    with open_fn(output_path, "wt") as fh:
        n_cells = merge(
            rows, placements, fh,
            keep_whole_image_annotations=keep_whole_image_annotations,
            classify_by_well=classify_by_well,
        )

    log(f"Merged {n_cells} cell(s) across {len(rows)} well(s); written to {output_path}")


if __name__ == "__main__":
    app()
