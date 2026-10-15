/*
 * This subworkflow uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */

include { SOPA_SEGMENTATIONCELLPOSE } from '../../../modules/local/sopa/segmentationcellpose/main.nf'
include { SOPA_RESOLVECELLPOSE      } from '../../../modules/local/sopa/resolvecellpose/main.nf'

workflow SOPA_CELLPOSE {

    take:
    ch_patches // channel: [ (meta, zarr, index, n_patches, nuclear_channel, membrane_channels) ]
    ch_spatial_data // channel: [ (meta, zarr) ]
    ch_cellpose_models // channel: value, path to the staged Cellpose weights directory

    main:

    ch_versions = channel.empty()

    //
    // Run SOPA segmentation with cellpose
    //
    SOPA_SEGMENTATIONCELLPOSE(
        ch_patches,
        ch_cellpose_models
    )
    ch_versions = ch_versions.mix(SOPA_SEGMENTATIONCELLPOSE.out.versions.first())


    // Collect cellpose segmentation boundaries into one channel per sample.
    // groupKey carries each sample's patch count, so a sample is emitted as
    // soon as its own patches finish instead of after every sample's patches.
    ch_n_patches = ch_patches
        .map { meta, _zarr, _index, n_patches, _nuclear_channel, _membrane_channels -> [ meta, n_patches ] }
        .unique()

    SOPA_SEGMENTATIONCELLPOSE.out.cellpose_parquet
        .combine( ch_n_patches, by: 0 )
        .map { meta, parquet, n_patches -> [ groupKey(meta, n_patches), parquet ] }
        .groupTuple(sort: true)
        .map { key, parquets -> [ key.getGroupTarget(), parquets ] }
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
    ch_versions = ch_versions.mix(SOPA_RESOLVECELLPOSE.out.versions.first())

    emit:
    boundaries  = SOPA_RESOLVECELLPOSE.out.cellpose_boundaries  // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]

    versions = ch_versions                                     // channel: [ versions.yml ]
}
