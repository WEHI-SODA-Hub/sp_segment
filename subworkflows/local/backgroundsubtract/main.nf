include { EXTRACTMARKERS      } from '../../../modules/local/extractmarkers/main.nf'
include { BACKSUB             } from '../../../modules/nf-core/backsub/main.nf'

workflow BACKGROUNDSUBTRACT {

    take:
    ch_backsub // channel: background subtraction parameters (sample name and tiff)

    main:

    ch_versions = channel.empty()

    //
    // Extract markers from the input tiff file
    //
    EXTRACTMARKERS(
        ch_backsub
    )
    ch_versions = ch_versions.mix(EXTRACTMARKERS.out.versions.first())


    //
    // Run background subtraction module on tiff with extracted markers
    //
    // BACKSUB takes the image and the markers file as two inputs, which
    // Nextflow pairs by arrival order. EXTRACTMARKERS tasks finish in any
    // order, so join on meta first: each markers file holds that slide's own
    // exposure times and background channels.
    ch_backsub
        .join(EXTRACTMARKERS.out.markers, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .multiMap { meta, tiff, markers ->
            image:   [ meta, tiff ]
            markers: [ meta, markers ]
        }
        .set { ch_backsub_input }

    BACKSUB(
        ch_backsub_input.image,
        ch_backsub_input.markers
    )
    ch_versions = ch_versions.mix(BACKSUB.out.versions.first())

    emit:
    backsub_tif   = BACKSUB.out.backsub_tif    // channel: [ val(meta), *.ome.tif ]
    markers       = BACKSUB.out.markerout      // channel: [ val(meta2), markers.csv ]

    versions      = ch_versions                     // channel: [ versions.yml ]
}
