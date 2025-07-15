process PARQUETTOTIFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/tifffile_pyarrow_rasterio_shapely_pruned:9cab11ac36e81144'

    input:
    tuple val(meta), path(parquet), path(tiff), val(compartment)

    output:
    tuple val(meta), path("*.tiff"), emit: tiff
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    parquet_to_tiff.py \\
        $args \\
        $parquet \\
        $tiff > ${prefix}_${compartment}.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${compartment}.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
