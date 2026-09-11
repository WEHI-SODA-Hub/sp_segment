process AVITIASSEMBLEPLATE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Same tifffile/numpy/typer image already used by AVITISTITCHWELL/PARQUETTOTIFF --
    // the plate assembler is pure tifffile/numpy streaming, no new dependency.
    container 'community.wave.seqera.io/library/tifffile_pyarrow_rasterio_shapely_pruned:9cab11ac36e81144'

    input:
    // well_rows is a list of maps ([well, image]) built by the subworkflow;
    // images are the corresponding staged per-well stitched intensity TIFFs,
    // named exactly as referenced in well_rows so the manifest built below
    // resolves them.
    tuple val(meta), val(well_rows), path(images)

    output:
    tuple val(meta), path("${meta.id}_plate.ome.tif")     , emit: image
    tuple val(meta), path("${meta.id}_plate_overview.tif"), emit: overview
    tuple val(meta), path("${meta.id}_plate_layout.csv")  , emit: layout
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Same last-flag-wins pattern as AVITISTITCHWELL: the run's own pixel
    // size (from RunParameters.json) must come after args so it overrides
    // the global params.pixel_size_microns in ext.args.
    def pixel_size_arg = meta.pixel_size_microns ? "--pixel-size-microns ${meta.pixel_size_microns}" : ''
    def manifest_lines = (
        ['well,image'] +
        well_rows.collect { r -> "${r.well},${r.image}" }
    ).join('\n    ')
    """
    cat > manifest.csv <<'AVITI_PLATE_MANIFEST_EOF'
    ${manifest_lines}
    AVITI_PLATE_MANIFEST_EOF

    aviti_assemble_plate_image.py \\
        manifest.csv \\
        --output-image ${meta.id}_plate.ome.tif \\
        --output-layout ${meta.id}_plate_layout.csv \\
        --output-overview ${meta.id}_plate_overview.tif \\
        ${args} \\
        ${pixel_size_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    // Layout needs a real header + at least one data row, not an empty
    // touch(): AVITIMERGEPLATEGEOJSON's stub consumes this file and has
    // nothing to work with otherwise (same reasoning as AVITIDISCOVERTILES'
    // stub emitting real manifest rows).
    def first_well = well_rows[0].well
    """
    touch ${meta.id}_plate.ome.tif
    touch ${meta.id}_plate_overview.tif
    cat <<-END_LAYOUT > ${meta.id}_plate_layout.csv
    well,col,row,x0,y0,width,height
    ${first_well},0,0,0,0,8,8
    END_LAYOUT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
