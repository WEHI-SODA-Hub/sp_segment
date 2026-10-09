"""Unit tests for ``bin/aviti_merge_plate_geojson.py``.

Covers the streaming feature scanner (the part that lets this avoid
``json.load`` on a potentially gigabyte-scale GeoJSON), coordinate
translation, per-well id offsetting, and an end-to-end merge.
"""

import gzip
import io
import json

import pytest

from aviti_merge_plate_geojson import (
    WHOLE_IMAGE_ANNOTATION_ID,
    iter_geojson_features,
    merge,
    read_well_rows,
    transform_feature,
    translate_coordinates,
    well_annotation_feature,
)
from aviti_plate_common import PlacedWell


# --- iter_geojson_features ---------------------------------------------


def _features_of(text):
    return list(iter_geojson_features(io.StringIO(text)))


def test_iter_geojson_features_yields_each_top_level_feature():
    doc = json.dumps({
        "type": "FeatureCollection",
        "features": [{"id": "a", "properties": {}}, {"id": "b", "properties": {}}],
    })
    features = [json.loads(t) for t in _features_of(doc)]
    assert [f["id"] for f in features] == ["a", "b"]


def test_iter_geojson_features_handles_pretty_printed_input():
    doc = json.dumps({
        "type": "FeatureCollection",
        "features": [{"id": "a"}, {"id": "b"}],
    }, indent=2)
    features = [json.loads(t) for t in _features_of(doc)]
    assert [f["id"] for f in features] == ["a", "b"]


def test_iter_geojson_features_handles_minified_input():
    doc = '{"type":"FeatureCollection","features":[{"id":"a"},{"id":"b"}]}'
    features = [json.loads(t) for t in _features_of(doc)]
    assert [f["id"] for f in features] == ["a", "b"]


def test_iter_geojson_features_handles_empty_features_array():
    doc = '{"type":"FeatureCollection","features":[]}'
    assert _features_of(doc) == []


def test_iter_geojson_features_does_not_split_on_braces_inside_strings():
    # A property value containing '}', ']' and an embedded double-quote must
    # not be mistaken for the end of the feature object. Building the
    # expected value as a plain variable (rather than hand-escaping it twice)
    # keeps this test's own intent unambiguous.
    note = 'has a } and a ] and a "quoted" bit'
    doc = json.dumps({
        "type": "FeatureCollection",
        "features": [
            {"id": "a", "properties": {"note": note}},
            {"id": "b", "properties": {}},
        ],
    })
    features = [json.loads(t) for t in _features_of(doc)]
    assert [f["id"] for f in features] == ["a", "b"]
    assert features[0]["properties"]["note"] == note


def test_iter_geojson_features_handles_a_backslash_at_the_end_of_a_string():
    # The classic edge case for a naive string scanner: a value ending in a
    # backslash right before the closing quote must not be read as an
    # escaped quote.
    path = "C:\\"
    doc = json.dumps({
        "type": "FeatureCollection",
        "features": [{"id": "a", "properties": {"path": path}}],
    })
    features = [json.loads(t) for t in _features_of(doc)]
    assert features[0]["properties"]["path"] == path


def test_iter_geojson_features_ignores_a_features_like_key_before_the_real_one():
    # Guard against the substring search in _find_features_array_start
    # matching a decoy "features" token that appears before the real array.
    doc = '{"type":"has \\"features\\": [1,2] inside a string","features":[{"id":"a"}]}'
    features = [json.loads(t) for t in _features_of(doc)]
    assert [f["id"] for f in features] == ["a"]


def test_iter_geojson_features_streams_across_small_chunk_boundaries(monkeypatch):
    # Force tiny reads so multi-chunk feature objects are exercised even on a
    # small fixture, without needing a huge test document.
    import aviti_merge_plate_geojson as module

    monkeypatch.setattr(module, "_READ_CHUNK", 3)
    doc = json.dumps({
        "type": "FeatureCollection",
        "features": [{"id": "a", "properties": {"n": 12345}}, {"id": "b"}],
    })
    features = [json.loads(t) for t in _features_of(doc)]
    assert [f["id"] for f in features] == ["a", "b"]
    assert features[0]["properties"]["n"] == 12345


# --- translate_coordinates --------------------------------------------


def test_translate_coordinates_polygon():
    ring = [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]]
    out = translate_coordinates(ring, dx=10, dy=20)
    assert out == [[[10, 20], [11, 20], [11, 21], [10, 21], [10, 20]]]


def test_translate_coordinates_polygon_with_hole():
    polygon = [
        [[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]],
        [[2, 2], [4, 2], [4, 4], [2, 4], [2, 2]],
    ]
    out = translate_coordinates(polygon, dx=5, dy=5)
    assert out[0][0] == [5, 5]
    assert out[1][0] == [7, 7]


def test_translate_coordinates_multipolygon():
    multi = [
        [[[0, 0], [1, 0], [1, 1], [0, 0]]],
        [[[2, 2], [3, 2], [3, 3], [2, 2]]],
    ]
    out = translate_coordinates(multi, dx=100, dy=200)
    assert out[0][0][0] == [100, 200]
    assert out[1][0][0] == [102, 202]


def test_translate_coordinates_preserves_floats():
    out = translate_coordinates([1.5, 2.25], dx=0.5, dy=0.75)
    assert out == [2.0, 3.0]


# --- transform_feature ---------------------------------------------------


def _cell_feature(cell_id, nucleus_label, whole_cell_label, x=1.0, y=2.0):
    return {
        "type": "Feature",
        "id": f"cell-{cell_id}",
        "geometry": {"type": "Polygon", "coordinates": [[[x, y], [x + 1, y], [x + 1, y + 1], [x, y]]]},
        "nucleusGeometry": {"type": "Polygon", "coordinates": [[[x, y], [x, y], [x, y], [x, y]]]},
        "properties": {
            "objectType": "cell", "id": cell_id,
            "nucleus_label": nucleus_label, "whole_cell_label": whole_cell_label,
        },
    }


def test_transform_feature_shifts_both_geometries_by_the_same_offset():
    feature = _cell_feature(1, 10, 20, x=5.0, y=6.0)
    out = transform_feature(feature, "A1", dx=100, dy=200, id_offsets={"id": 0, "nucleus_label": 0, "whole_cell_label": 0})
    assert out["geometry"]["coordinates"][0][0] == [105.0, 206.0]
    assert out["nucleusGeometry"]["coordinates"][0][0] == [105.0, 206.0]


def test_transform_feature_offsets_the_three_id_namespaces_independently():
    feature = _cell_feature(cell_id=1, nucleus_label=50, whole_cell_label=7)
    out = transform_feature(
        feature, "B1", dx=0, dy=0,
        id_offsets={"id": 100, "nucleus_label": 9, "whole_cell_label": 3},
    )
    assert out["properties"]["id"] == 101
    assert out["properties"]["nucleus_label"] == 59
    assert out["properties"]["whole_cell_label"] == 10


def test_transform_feature_prefixes_the_string_id_with_the_well_and_tags_well():
    feature = _cell_feature(1, 1, 1)
    out = transform_feature(feature, "C3", dx=0, dy=0, id_offsets={})
    assert out["id"] == "C3-cell-1"
    assert out["properties"]["well"] == "C3"


def test_transform_feature_does_not_mutate_the_input():
    feature = _cell_feature(1, 1, 1)
    original = json.loads(json.dumps(feature))
    transform_feature(feature, "A1", dx=5, dy=5, id_offsets={"id": 10})
    assert feature == original


# --- well_annotation_feature ---------------------------------------------


def test_well_annotation_feature_matches_the_layout_rectangle_exactly():
    placed = PlacedWell(well="A1", col=0, row=0, x0=10, y0=20, width=100, height=50)
    feature = well_annotation_feature(placed)
    assert feature["properties"]["objectType"] == "annotation"
    assert feature["properties"]["name"] == "A1"
    ring = feature["geometry"]["coordinates"][0]
    assert ring[0] == [10, 20]
    assert ring[2] == [110, 70]  # opposite corner: x0+width, y0+height
    assert ring[0] == ring[-1]  # closed ring


# --- read_well_rows --------------------------------------------------------


def test_read_well_rows_orders_by_plate_position(tmp_path):
    manifest = tmp_path / "m.csv"
    manifest.write_text("well,geojson\nB1,b1.geojson\nA1,a1.geojson\nA2,a2.geojson\n")
    assert [r["well"] for r in read_well_rows(manifest)] == ["A1", "A2", "B1"]


# --- end-to-end merge ------------------------------------------------------


def _write_geojson(path, features, gzip_it=False):
    doc = json.dumps({"type": "FeatureCollection", "features": features})
    if gzip_it:
        with gzip.open(path, "wt") as fh:
            fh.write(doc)
    else:
        path.write_text(doc)


def _whole_image_annotation(w=10, h=10):
    return {
        "type": "Feature", "id": WHOLE_IMAGE_ANNOTATION_ID,
        "geometry": {"type": "Polygon", "coordinates": [[[0, 0], [w, 0], [w, h], [0, h], [0, 0]]]},
        "properties": {"objectType": "annotation", "type": "annotation", "name": "whole_image"},
    }


def test_merge_produces_unique_ids_across_wells_and_drops_whole_image_annotations(tmp_path):
    a1_path = tmp_path / "a1.geojson"
    b1_path = tmp_path / "b1.geojson"
    _write_geojson(a1_path, [
        _whole_image_annotation(),
        _cell_feature(1, nucleus_label=1, whole_cell_label=1),
        _cell_feature(2, nucleus_label=2, whole_cell_label=2),
    ])
    _write_geojson(b1_path, [
        _whole_image_annotation(),
        _cell_feature(1, nucleus_label=1, whole_cell_label=1),
    ])

    placements = {
        "A1": PlacedWell("A1", 0, 0, x0=0, y0=0, width=10, height=10),
        "B1": PlacedWell("B1", 1, 0, x0=20, y0=0, width=10, height=10),
    }
    rows = [{"well": "A1", "geojson": str(a1_path)}, {"well": "B1", "geojson": str(b1_path)}]

    out = io.StringIO()
    n_cells = merge(rows, placements, out)
    result = json.loads(out.getvalue())

    assert n_cells == 3
    cells = [f for f in result["features"] if f["properties"].get("objectType") == "cell"]
    assert len(cells) == 3
    # No duplicate properties.id, nucleus_label or whole_cell_label across wells.
    assert len({c["properties"]["id"] for c in cells}) == 3
    assert len({c["properties"]["nucleus_label"] for c in cells}) == 3
    assert len({c["properties"]["whole_cell_label"] for c in cells}) == 3
    assert len({c["id"] for c in cells}) == 3
    # Whole-image annotations were dropped; well rectangles are present instead.
    assert not any(f["id"] == WHOLE_IMAGE_ANNOTATION_ID for f in result["features"])
    well_features = [f for f in result["features"] if f["properties"].get("name") in ("A1", "B1")]
    assert len(well_features) == 2


def test_merge_writes_well_annotations_before_any_cells(tmp_path):
    a1_path = tmp_path / "a1.geojson"
    _write_geojson(a1_path, [_cell_feature(1, 1, 1)])
    placements = {"A1": PlacedWell("A1", 0, 0, x0=0, y0=0, width=10, height=10)}
    rows = [{"well": "A1", "geojson": str(a1_path)}]

    out = io.StringIO()
    merge(rows, placements, out)
    result = json.loads(out.getvalue())

    assert result["features"][0]["properties"]["objectType"] == "annotation"
    assert result["features"][-1]["properties"]["objectType"] == "cell"


def test_merge_keeps_whole_image_annotations_when_asked(tmp_path):
    a1_path = tmp_path / "a1.geojson"
    _write_geojson(a1_path, [_whole_image_annotation(), _cell_feature(1, 1, 1)])
    placements = {"A1": PlacedWell("A1", 0, 0, x0=0, y0=0, width=10, height=10)}
    rows = [{"well": "A1", "geojson": str(a1_path)}]

    out = io.StringIO()
    merge(rows, placements, out, keep_whole_image_annotations=True)
    result = json.loads(out.getvalue())
    # Like every other feature, a kept whole-image annotation still goes
    # through transform_feature and gets its id well-prefixed (so the
    # original id string no longer appears verbatim) -- check by content
    # (properties.name) instead of by the original id.
    assert any(
        f["properties"].get("name") == "whole_image" for f in result["features"]
    )
    assert not any(f["id"] == WHOLE_IMAGE_ANNOTATION_ID for f in result["features"])


def test_merge_output_is_valid_json():
    out = io.StringIO()
    merge([], {}, out)
    result = json.loads(out.getvalue())
    assert result == {"type": "FeatureCollection", "features": []}


def test_merge_reads_gzipped_well_geojson_transparently(tmp_path):
    a1_path = tmp_path / "a1.geojson.gz"
    _write_geojson(a1_path, [_cell_feature(1, 1, 1)], gzip_it=True)
    placements = {"A1": PlacedWell("A1", 0, 0, x0=0, y0=0, width=10, height=10)}
    rows = [{"well": "A1", "geojson": str(a1_path)}]

    out = io.StringIO()
    n_cells = merge(rows, placements, out)
    assert n_cells == 1


def test_merge_rejects_a_well_missing_from_the_layout():
    with pytest.raises(KeyError):
        merge(
            [{"well": "Z9", "geojson": "does-not-matter.geojson"}],
            placements_by_well={},
            output=io.StringIO(),
        )


# --- CLI: --gzip output round-trip ----------------------------------------


def test_cli_gzip_output_matches_plain_output_when_decompressed(tmp_path):
    from typer.testing import CliRunner

    from aviti_merge_plate_geojson import app

    a1 = tmp_path / "a1.geojson"
    _write_geojson(a1, [_cell_feature(1, 1, 1)])
    layout = tmp_path / "layout.csv"
    layout.write_text("well,col,row,x0,y0,width,height\nA1,0,0,0,0,10,10\n")
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(f"well,geojson\nA1,{a1}\n")

    plain_out = tmp_path / "plain.geojson"
    gz_out = tmp_path / "gz.geojson"

    runner = CliRunner()
    r1 = runner.invoke(app, [str(manifest), "--layout", str(layout), "--output", str(plain_out), "--no-gzip"])
    assert r1.exit_code == 0, r1.output
    r2 = runner.invoke(app, [str(manifest), "--layout", str(layout), "--output", str(gz_out), "--gzip"])
    assert r2.exit_code == 0, r2.output

    assert plain_out.exists()
    gz_path = tmp_path / "gz.geojson.gz"
    assert gz_path.exists()
    with gzip.open(gz_path, "rt") as fh:
        gz_content = fh.read()
    assert plain_out.read_text() == gz_content
