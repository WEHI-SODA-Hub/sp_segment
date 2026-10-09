process AVITIDISCOVERTILES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Reuses the tifffile/pyarrow/typer image already pulled for
    // PARQUETTOTIFF: this module's own runtime needs is just typer (json/csv
    // are stdlib), and that container already satisfies it without pulling a
    // new one.
    container 'community.wave.seqera.io/library/tifffile_pyarrow_rasterio_shapely_pruned:9cab11ac36e81144'

    input:
    tuple val(meta), path(run_dir)

    output:
    tuple val(meta), path("${meta.id}.manifest.csv"), emit: manifest
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // meta.wells follows this pipeline's colon-separated list convention
    // (matching membrane_channels in the COMET/MIBI samplesheet), but the
    // underlying script takes a comma-separated list -- translate here
    // rather than adding a second list separator to the script's own CLI.
    def wells_arg = meta.wells ? "--wells '${meta.wells.replace(':', ',')}'" : ''
    """
    aviti_discover_tiles.py \
        "${run_dir}" \
        --output "${meta.id}.manifest.csv" \
        ${wells_arg} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        typer: \$(python3 -c "import typer; print(typer.__version__)")
    END_VERSIONS
    """

    stub:
    // Emits the same two tiles a real run against tests/data/aviti's WellA1
    // (RunParameters.json + Projection/WellA1/CP01_L1R01C0{1,2}S1_*) would
    // produce, for either one well (A1) or two (A1 + a second, "A2", row
    // reusing the same real files purely for staging -- their content is
    // never read in stub mode). Branching on meta.id here (rather than on
    // run_dir's contents, which is a small fixture with only one real well)
    // is test-only scaffolding: it lets downstream stub-mode tests --
    // subworkflows/local/aviti_segment in particular -- exercise both a
    // genuinely single-well sample (plate assembly must be skipped) and a
    // multi-well one (plate assembly must run) without needing a second run
    // directory fixture. Real (non-stub) discovery always reflects whatever
    // RunParameters.json actually contains.
    def wells = (meta.id ?: '').contains('single_well') ? ['A1'] : ['A1', 'A2']
    def manifest_lines = (
        ['well,tile,x_mm,y_mm,nucleus_tif,membrane_tif,actin_tif,channel_mode,pixel_size_microns'] +
        wells.collectMany { well ->
            [
                "${well},L1R01C01S1,0.0,0.0,${run_dir}/Projection/WellA1/CP01_L1R01C01S1_Nucleus.tif,${run_dir}/Projection/WellA1/CP01_L1R01C01S1_Cell-Membrane.tif,${run_dir}/Projection/WellA1/CP01_L1R01C01S1_Actin.tif,3ch,0.48",
                "${well},L1R01C02S1,0.001,0.0,${run_dir}/Projection/WellA1/CP01_L1R01C02S1_Nucleus.tif,${run_dir}/Projection/WellA1/CP01_L1R01C02S1_Cell-Membrane.tif,${run_dir}/Projection/WellA1/CP01_L1R01C02S1_Actin.tif,3ch,0.48",
            ]
        }
    ).join('\n    ')
    """
    cat > ${meta.id}.manifest.csv <<'AVITI_STUB_MANIFEST_EOF'
    ${manifest_lines}
    AVITI_STUB_MANIFEST_EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        typer: \$(python3 -c "import typer; print(typer.__version__)")
    END_VERSIONS
    """
}
