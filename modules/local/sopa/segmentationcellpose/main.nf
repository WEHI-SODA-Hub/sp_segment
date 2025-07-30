/*
 * This module uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */
process SOPA_SEGMENTATIONCELLPOSE {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.0.7-cellpose'
        : 'docker.io/quentinblampey/sopa:2.0.7-cellpose'}"

    input:
    tuple val(meta), path(zarr), val(index), val(n_patches), val(nuclear_channel), val(membrane_channels)

    output:
    tuple val(meta), path("*.zarr/.sopa_cache/cellpose_boundaries/${index}.parquet"), emit: cellpose_parquet

    script:
    def args = task.ext.args ?: ''
    def membrane_channel_args = membrane_channels && membrane_channels != "" ?
        membrane_channels.split(":").collect(
        { "--channels ${it}" }
    ).join(' ') : ''
    def channels = "--channels ${nuclear_channel} ${membrane_channel_args}"
    """
    sopa segmentation cellpose \\
        ${args} \\
        --patch-index ${index} \\
        ${channels} \\
        --diameter ${params.cellpose_diameter} \\
        --min-area ${params.cellpose_min_area} \\
        ${zarr}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.zarr/.sopa_cache/cellpose_boundaries
    touch ${prefix}.zarr/.sopa_cache/cellpose_boundaries/${index}.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """
}
