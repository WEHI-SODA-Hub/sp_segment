/*
 * This subworkflow uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */

process SOPASEGMENTATIONCELLPOSE {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.0.7-cellpose'
        : 'docker.io/quentinblampey/sopa:2.0.7-cellpose'}"

    input:
    tuple val(meta), path(zarr), val(index), val(n_patches), val(nuclear_channel)

    output:
    tuple val(meta), path("*.zarr/.sopa_cache/cellpose_boundaries/${index}.parquet"), emit: cellpose_parquet

    script:
    def args = task.ext.args ?: ''
    """
    sopa segmentation cellpose \\
        ${args} \\
        --patch-index ${index} \\
        --channels ${nuclear_channel} \\
        --diameter ${params.cellpose_diameter} \\
        --min-area ${params.cellpose_min_area} \\
        ${zarr}
    """
}

process SOPARESOLVECELLPOSE {
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

workflow SOPACELLPOSE {

    take:
    ch_patches // channel: [ (meta, zarr, index, n_patches, nuclear_channel) ]
    ch_spatial_data // channel: [ (meta, zarr) ]

    main:

    ch_versions = Channel.empty()

    //
    // Run SOPA segmentation with cellpose
    //
    SOPASEGMENTATIONCELLPOSE(
        ch_patches
    )

    // Collect cellpose segmentation boundaries into one channel per sample
    SOPASEGMENTATIONCELLPOSE.out.cellpose_parquet
        .groupTuple()
        .join( ch_spatial_data, by: 0 )
        .map { meta, cellpose_parquet, zarr ->
            [ meta, zarr, cellpose_parquet ]
        }
        .set { ch_resolve_cellpose }

    //
    // Resolve Cellpose segmentation boundaries
    //
    SOPARESOLVECELLPOSE(
        ch_resolve_cellpose
    )

    emit:
    boundaries  = SOPARESOLVECELLPOSE.out.cellpose_boundaries  // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]

    versions = ch_versions                                     // channel: [ versions.yml ]
}
