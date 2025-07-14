/*
 * This subworkflow uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */

process SOPACONVERT {
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.0.7'
        : 'docker.io/quentinblampey/sopa:2.0.7'}"

    input:
    tuple val(meta), path(tiff)

    output:
    tuple val(meta), path("*.zarr"), emit: spatial_data

    script:
    def args = task.ext.args ?: ''
    """
    sopa convert \\
        ${args} \\
        --sdata-path ${meta.id}.zarr \\
        --technology ${params.technology} \\
        ${tiff}
    """
}

process SOPAPATCHIFYIMAGE {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.0.7'
        : 'docker.io/quentinblampey/sopa:2.0.7'}"

    input:
    tuple val(meta), path(zarr)

    output:
    tuple val(meta), path("*.zarr/.sopa_cache/patches_file_image"), path("*.zarr/shapes/image_patches"), emit: patches

    script:
    def args = task.ext.args ?: ''
    """
    sopa patchify image \\
        ${args} \\
        ${zarr} \\
        --patch-width-pixel ${params.patch_width_pixel} \\
        --patch-overlap-pixel ${params.patch_overlap_pixel}
    """
}

process SOPASEGMENTATIONCELLPOSENUCLEAR {
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

process SOPARESOLVECELLPOSENUCLEAR {
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

workflow SOPASEGMENT {

    take:
    ch_sopa // channel: [ (meta, tiff, nuclear_channel, membrane_channels) ]

    main:

    ch_versions = Channel.empty()

    ch_sopa.map { meta, tiff, nuclear_channel, membrane_channels, skip_measurements ->
        [ meta, tiff ]
    }.set { ch_convert }


    //
    // Run SOPA convert to convert tiff to zarr format
    //
    SOPACONVERT(
        ch_convert
    )

    //
    // Run SOPA patchify to create image patches
    //
    SOPAPATCHIFYIMAGE(
        SOPACONVERT.out.spatial_data
    )

    // Create a channel for each patch
    SOPAPATCHIFYIMAGE.out.patches
        .join( SOPACONVERT.out.spatial_data, by: 0 )
        .map { meta, patches_file_image, image_patches, zarr ->
            [ meta, zarr, patches_file_image.text.trim().toInteger() ] }
        .flatMap { meta, zarr, n_patches ->
            (0..<n_patches).collect { index -> [ meta, zarr, index, n_patches ] } }
        .combine(ch_sopa, by: 0)
        .map { meta, zarr, index, n_patches, tiff, nuclear_channel, membrane_channels, skip_measurements ->
            [ meta, zarr, index, n_patches, nuclear_channel.first() ]
        }.set { ch_cellpose }

    //
    // Run SOPA segmentation with cellpose
    //
    SOPASEGMENTATIONCELLPOSENUCLEAR(
        ch_cellpose
    )

    // Collect cellpose segmentation boundaries into one channel per sample
    SOPASEGMENTATIONCELLPOSENUCLEAR.out.cellpose_parquet
        .groupTuple()
        .join( SOPACONVERT.out.spatial_data, by: 0 )
        .map { meta, cellpose_parquet, zarr ->
            [ meta, zarr, cellpose_parquet ]
        }
        .set { ch_resolve_cellpose }

    //
    // Resolve Cellpose segmentation boundaries
    //
    SOPARESOLVECELLPOSENUCLEAR(
        ch_resolve_cellpose
    )

    emit:
    zarr        = SOPACONVERT.out.spatial_data                        // channel: [ val(meta), *.zarr ]
    boundaries  = SOPARESOLVECELLPOSENUCLEAR.out.cellpose_boundaries  // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]

    versions = ch_versions                                            // channel: [ versions.yml ]
}
