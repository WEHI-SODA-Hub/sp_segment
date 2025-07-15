/*
 * This subworkflow uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */

process SOPA_SEGMENTATIONCELLPOSE {
    label "process_single"

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
    def membrane_channel_args = membrane_channels ? membrane_channels.split(":").collect(
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
    """
}

process SOPA_RESOLVECELLPOSE {
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.0.7'
        : 'docker.io/quentinblampey/sopa:2.0.7'}"

    input:
    tuple val(meta), path(zarr), val(cellpose_parquet)

    output:
    tuple val(meta), path("*.zarr/shapes/cellpose_boundaries/*.parquet"), emit: cellpose_boundaries

    when:

    script:
    """
    sopa resolve cellpose ${zarr}
    """
}

workflow SOPA_CELLPOSE {

    take:
    ch_patches // channel: [ (meta, zarr, index, n_patches, nuclear_channel, membrane_channels) ]
    ch_spatial_data // channel: [ (meta, zarr) ]

    main:

    ch_versions = Channel.empty()

    //
    // Run SOPA segmentation with cellpose
    //
    SOPA_SEGMENTATIONCELLPOSE(
        ch_patches
    )

    // Collect cellpose segmentation boundaries into one channel per sample
    SOPA_SEGMENTATIONCELLPOSE.out.cellpose_parquet
        .groupTuple()
        .join( ch_spatial_data, by: 0 )
        .map { meta, cellpose_parquet, zarr ->
            [ meta, zarr, cellpose_parquet ]
        }
        .set { ch_resolve_cellpose }

    //
    // Resolve Cellpose segmentation boundaries
    //
    SOPA_RESOLVECELLPOSE(
        ch_resolve_cellpose
    )

    emit:
    boundaries  = SOPA_RESOLVECELLPOSE.out.cellpose_boundaries  // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]

    versions = ch_versions                                     // channel: [ versions.yml ]
}
