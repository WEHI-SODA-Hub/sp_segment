process AVITIMERGEPLATEGEOJSON {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Same tifffile/numpy/typer image already used by AVITIASSEMBLEPLATE --
    // the merge script itself only needs typer/stdlib, but reuses the
    // already-pulled container rather than adding a new one.
    container 'community.wave.seqera.io/library/tifffile_pyarrow_rasterio_shapely_pruned:9cab11ac36e81144'

    input:
    // well_rows is a list of maps ([well, geojson]) built by the subworkflow;
    // geojsons are the corresponding staged per-well cellmeasurement
    // annotation files, named exactly as referenced in well_rows so the
    // manifest built below resolves them. layout is the plate layout CSV
    // AVITIASSEMBLEPLATE wrote, so both plate artefacts place every well at
    // byte-identical coordinates.
    tuple val(meta), val(well_rows), path(geojsons), path(layout)

    output:
    tuple val(meta), path("${meta.id}_plate.geojson{,.gz}"), emit: annotations
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def manifest_lines = (
        ['well,geojson'] +
        well_rows.collect { r -> "${r.well},${r.geojson}" }
    ).join('\n    ')
    """
    cat > manifest.csv <<'AVITI_PLATE_GEOJSON_MANIFEST_EOF'
    ${manifest_lines}
    AVITI_PLATE_GEOJSON_MANIFEST_EOF

    aviti_merge_plate_geojson.py \\
        manifest.csv \\
        --layout ${layout} \\
        --output ${meta.id}_plate.geojson \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    // Always the plain filename, regardless of params.gzip_geojson: an empty
    // (touch'd) file named *.gz is not valid gzip, and nf-test's own md5
    // snapshotting of a "-stub" run tries to decompress *.gz outputs to hash
    // them, which throws on a zero-byte file. Matches CELLMEASUREMENT's own
    // stub, which has the same *.geojson{,.gz} output pattern and touches
    // only the unsuffixed name for the same reason. The output glob still
    // matches either way.
    """
    touch ${meta.id}_plate.geojson

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
