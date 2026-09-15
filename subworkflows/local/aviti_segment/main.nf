include { AVITIDISCOVERTILES     } from '../../../modules/local/aviti/discovertiles/main.nf'
include { AVITIMERGETILECHANNELS } from '../../../modules/local/aviti/mergetilechannels/main.nf'
include { AVITIWHOLECELLSEGMENT  } from '../../../modules/local/aviti/wholecellsegment/main.nf'
include { AVITIMEMBRANESEGMENT   } from '../../../modules/local/aviti/membranesegment/main.nf'
include { AVITINUCLEARSEGMENT    } from '../../../modules/local/aviti/nuclearsegment/main.nf'
include { AVITISTITCHWELL        } from '../../../modules/local/aviti/stitchwell/main.nf'
include { AVITIASSEMBLEPLATE     } from '../../../modules/local/aviti/assembleplate/main.nf'
include { AVITIMERGEPLATEGEOJSON } from '../../../modules/local/aviti/mergeplategeojson/main.nf'
include { CELLMEASUREMENT         } from '../../../modules/local/cellmeasurement/main.nf'
include { SEGMENTATIONREPORT      } from '../../../modules/local/segmentationreport/main.nf'

// Top-level functions -- need to be defined here; if they are inside the
// workflow { } block, they are treated as out-of scope in sibling
// .map/.filter calls in Nextflow's strict syntax parser (24.10+/26.x)
def plate_meta_of(well_meta) {
    [
        id                 : "${well_meta.sample}__plate",
        sample             : well_meta.sample,
        channel_mode       : well_meta.channel_mode,
        pixel_size_microns : well_meta.pixel_size_microns,
        aviti_diameter     : well_meta.aviti_diameter,
        aviti_model        : well_meta.aviti_model,
    ]
}

// Top-level for the same reason as plate_meta_of above.
def to_report_input(meta, annotations, image) {
    [
        meta,
        annotations,
        false, // run_mesmer
        true,  // run_cellpose
        false, // run_cellsam
        'Nucleus',
        meta.channel_mode == '2ch' ? 'Cell-Membrane' : 'Cell-Membrane:Actin',
        image,
    ]
}

workflow AVITI_SEGMENT {

    take:
    ch_aviti_samplesheet   // channel: [ meta, run_dir ] -- meta.id is the sample name; meta.wells (optional) restricts --wells
    ch_nuclear_model       // channel: value, path to the staged Cellpose 3.x custom nuclear model
    ch_membrane_models_dir // channel: value, path to the staged directory of Cellpose 3.x membrane model files (empty if no row uses membrane_model)

    main:

    ch_versions = channel.empty()

    //
    // Discover wells/tiles for each AVITI sample and write a flat per-tile
    // manifest (one row per (well, tile)), including the stage coordinates
    // stitching needs later.
    //
    AVITIDISCOVERTILES(
        ch_aviti_samplesheet
    )
    ch_versions = ch_versions.mix(AVITIDISCOVERTILES.out.versions.first())

    //
    // Fan the manifest out into one channel item per tile. This is the point
    // where per-sample AVITI processing becomes per-tile parallel: every tile
    // in every well of every sample becomes an independent task from here
    // until the stitching step groups them back by well.
    //
    // actin_tif is emitted as a literal "NO_FILE" placeholder when empty
    // (2-channel mode), since a Nextflow tuple element cannot be "no path" --
    // downstream modules detect this by name and omit --actin-tif.
    //
    AVITIDISCOVERTILES.out.manifest
        .splitCsv(elem: 1, header: true)
        .combine(ch_aviti_samplesheet, by: 0)
        .map { sample_meta, row, run_dir ->
            def base = file(run_dir).parent
            def meta_tile = [
                id                 : "${sample_meta.id}__Well${row.well}__${row.tile}",
                sample             : sample_meta.id,
                well               : row.well,
                tile               : row.tile,
                x_mm               : row.x_mm,
                y_mm               : row.y_mm,
                channel_mode       : row.channel_mode,
                // The run's own ImageInfo.PixelSizeUm, discovered per run rather
                // than taken from the pipeline-wide default -- it drives the
                // stitch gap, the plate image's OME PhysicalSize (QuPath's scale
                // bar) and cellmeasurement's micron-based measurements.
                pixel_size_microns : row.pixel_size_microns ?: '',
                membrane_model     : sample_meta.membrane_model ?: '',
                cell_diameter      : sample_meta.cell_diameter ?: 0,
            ]
            [
                meta_tile,
                file("${base}/${row.nucleus_tif}"),
                file("${base}/${row.membrane_tif}"),
                row.actin_tif ? file("${base}/${row.actin_tif}") : file('NO_FILE'),
            ]
        }
        .set { ch_tiles }

    //
    // Merge the raw per-channel tiles into one named multi-channel OME-TIFF.
    // This carries un-normalised pixel values (unlike the Cellpose input
    // stacks built independently by the two segmentation modules below) --
    // it is the per-tile intensity image later stitched for
    // CELLMEASUREMENT/KRONOS2EMBEDDINGS.
    //
    AVITIMERGETILECHANNELS(
        ch_tiles
    )
    ch_versions = ch_versions.mix(AVITIMERGETILECHANNELS.out.versions.first())

    //
    // Whole-cell/membrane segmentation. Use the v4 SAM path when no
    // membrane model is specified, and the Cellpose 3.x membrane model path
    // when a per-row model name is present.
    //
    ch_tiles
        .branch { meta_tile, _nucleus_tif, _membrane_tif, _actin_tif ->
            v3: (meta_tile.membrane_model ?: '').trim() != ''
            v4: true
        }
        .set { ch_tiles_by_membrane_model }

    AVITIWHOLECELLSEGMENT(
        ch_tiles_by_membrane_model.v4
    )
    ch_versions = ch_versions.mix(AVITIWHOLECELLSEGMENT.out.versions.first())

    // model_path is resolved and carried in the same tuple as meta/tiles
    // (one map, one output tuple) rather than split into two separately
    // mapped channels re-paired positionally at the process call -- the
    // v3 branch predicate above already guarantees membrane_model is
    // non-blank here, so there is no "no model" case left to filter out.
    ch_membrane_input = ch_tiles_by_membrane_model.v3
        .combine(ch_membrane_models_dir)
        .map { meta_tile, nucleus_tif, membrane_tif, actin_tif, models_dir ->
            def model_path = file("${models_dir}/${meta_tile.membrane_model.trim()}", checkIfExists: true)
            [ meta_tile, nucleus_tif, membrane_tif, actin_tif, model_path ]
        }

    AVITIMEMBRANESEGMENT(
        ch_membrane_input
    )
    ch_versions = ch_versions.mix(AVITIMEMBRANESEGMENT.out.versions.first())

    //
    // Nuclear segmentation with the Cellpose 3.x custom model, in its own
    // environment/container -- Cellpose 3.x and 4.x cannot coexist in one.
    //
    AVITINUCLEARSEGMENT(
        ch_tiles.map { meta_tile, nucleus_tif, _membrane_tif, _actin_tif -> [ meta_tile, nucleus_tif ] },
        ch_nuclear_model
    )
    ch_versions = ch_versions.mix(AVITINUCLEARSEGMENT.out.versions.first())

    //
    // Group the three per-tile outputs back by well. Joining on meta_tile
    // (by: 0) is safe here because every branch above originates from the
    // same ch_tiles, so the same meta_tile map is used as the join key by
    // all three.
    //
    AVITIWHOLECELLSEGMENT.out.cell_mask
        .mix(AVITIMEMBRANESEGMENT.out.cell_mask)
        // Stitch the instance-labeled nuclear mask, not the binary
        // Nuclear.tif -- CELLMEASUREMENT needs per-nucleus instance IDs to
        // match nuclei against whole cells by centroid.
        .join(AVITINUCLEARSEGMENT.out.nuclear_label_mask, by: 0)
        .join(AVITIMERGETILECHANNELS.out.image, by: 0)
        .map { meta_tile, cell_mask, nuclear_label_mask, image ->
            // The diameter/model actually used for this sample's whole-cell
            // segmentation, carried through for SEGMENTATIONREPORT to
            // display
            def has_membrane_model = (meta_tile.membrane_model ?: '').trim() != ''
            def cell_diameter_override = (meta_tile.cell_diameter ?: 0).toFloat()
            def aviti_diameter = cell_diameter_override > 0
                ? cell_diameter_override
                : (has_membrane_model ? params.aviti_membrane_diameter : params.aviti_wholecell_diameter)
            def well_meta = [
                id                 : "${meta_tile.sample}__Well${meta_tile.well}",
                sample             : meta_tile.sample,
                well               : meta_tile.well,
                channel_mode       : meta_tile.channel_mode,
                pixel_size_microns : meta_tile.pixel_size_microns,
                aviti_diameter     : aviti_diameter,
                aviti_model        : has_membrane_model ? meta_tile.membrane_model.trim() : '',
            ]
            def row = [
                tile        : meta_tile.tile,
                x_mm        : meta_tile.x_mm,
                y_mm        : meta_tile.y_mm,
                cell_mask   : cell_mask.name,
                nuclear_mask: nuclear_label_mask.name,
                image_tif   : image.name,
            ]
            [ well_meta, row, cell_mask, nuclear_label_mask, image ]
        }
        // Collects each well's tile rows and matching files into parallel
        // lists -- exactly the shape AVITISTITCHWELL expects.
        .groupTuple(by: 0)
        // groupTuple emits in task-completion order, which would otherwise
        // make the four parallel lists' order (and so the task hash) vary
        // run to run. Sort all four together by tile -- rather than
        // groupTuple's own `sort:`, which would sort each list independently
        // and de-pair a row from its own files -- keeps -resume stable.
        .map { well_meta, rows, cell_masks, nuclear_label_masks, images ->
            def pairs = [ rows, cell_masks, nuclear_label_masks, images ]
                .transpose()
                .sort { a, b -> a[0].tile <=> b[0].tile }
            [
                well_meta,
                pairs.collect { pair -> pair[0] },
                pairs.collect { pair -> pair[1] },
                pairs.collect { pair -> pair[2] },
                pairs.collect { pair -> pair[3] },
            ]
        }
        .set { ch_well_stitch_input }

    //
    // Stitch per-tile masks/images into per-well outputs, using the stage
    // coordinates carried in each tile row.
    //
    AVITISTITCHWELL(
        ch_well_stitch_input
    )
    ch_versions = ch_versions.mix(AVITISTITCHWELL.out.versions.first())

    //
    // Feed stitched per-well outputs into the existing, unmodified
    // CELLMEASUREMENT module. meta.id is well-scoped
    // (`<sample>__Well<well>`), keeping AVITI outputs distinct from any
    // COMET/MIBI sample sharing the sample name.
    //
    AVITISTITCHWELL.out.image
        .join(AVITISTITCHWELL.out.nuclear_mask, by: 0)
        .join(AVITISTITCHWELL.out.cell_mask, by: 0)
        .set { ch_cellmeasurement }

    CELLMEASUREMENT(
        ch_cellmeasurement
    )
    ch_versions = ch_versions.mix(CELLMEASUREMENT.out.versions.first())

    ch_annotations = CELLMEASUREMENT.out.annotations

    //
    // Assemble the KRONOS input channel: stitched well image + whole-cell
    // mask. KRONOS itself is invoked once at the top level
    // (workflows/sp_segment.nf), same as every other segmenter. This and
    // `annotations` above are deliberately left well-scoped only -- KRONOS's
    // top-level join keys on this exact meta shape, so plate-level artefacts
    // are never mixed into either.
    //
    AVITISTITCHWELL.out.image
        .join(AVITISTITCHWELL.out.cell_mask, by: 0)
        .set { ch_kronos_input }

    //
    // Group per-well outputs by sample under a plate-scoped meta, and count
    // wells per sample while we're at it. Computed unconditionally (it's a
    // cheap channel-only transform, not a process) because both the plate
    // assembly block and the report-scope logic below need the *true*
    // per-sample well count -- not just whether the samplesheet row
    // restricted `wells`, which says nothing about how many wells a run
    // actually has (an unrestricted row can still resolve to one well).
    //
    // groupTuple emits in task-completion order, which would otherwise make
    // the well_rows/images lists (and so the task hash) vary run to run.
    // Sorting the two parallel lists together by well -- rather than
    // groupTuple's own `sort:`, which would sort each list independently and
    // de-pair rows from files -- keeps -resume stable.
    AVITISTITCHWELL.out.image
        .map { well_meta, image -> [ plate_meta_of(well_meta), [ well: well_meta.well, image: image.name ], image ] }
        .groupTuple(by: 0)
        .map { plate_meta, well_rows, images ->
            def pairs = [ well_rows, images ].transpose().sort { a, b -> a[0].well <=> b[0].well }
            [ plate_meta, pairs.collect { pair -> pair[0] }, pairs.collect { pair -> pair[1] } ]
        }
        .set { ch_plate_image_input_all }

    ch_plate_image_input_all
        .map { plate_meta, well_rows, _images -> [ plate_meta.sample, well_rows.size() ] }
        .set { ch_well_count }

    //
    // Plate-level assembly (additive, optional): combine every well's
    // stitched image / cellmeasurement GeoJSON for a sample into one
    // pyramidal plate OME-TIFF and one merged plate GeoJSON, so a whole run
    // can be opened in QuPath at once. Does not touch the per-well outputs
    // above in any way -- new artefacts on new emits, keyed on a distinct
    // "<sample>__plate" meta.
    //
    ch_plate_image       = channel.empty()
    ch_plate_layout      = channel.empty()
    ch_plate_annotations = channel.empty()

    if (params.aviti_plate_assembly) {

        // A sample with only one well has nothing to combine -- the "plate"
        // would just be that well's own image again, at the cost of a full
        // pyramidal-image write/read pass. Skip it.
        ch_plate_image_input_all
            .filter { _plate_meta, well_rows, _images -> well_rows.size() > 1 }
            .set { ch_plate_image_input }

        AVITIASSEMBLEPLATE(
            ch_plate_image_input
        )
        ch_versions     = ch_versions.mix(AVITIASSEMBLEPLATE.out.versions.first())
        ch_plate_image  = AVITIASSEMBLEPLATE.out.image
        ch_plate_layout = AVITIASSEMBLEPLATE.out.layout

        // The layout comes from AVITIASSEMBLEPLATE, not recomputed here, so
        // the image and the merged GeoJSON place every well at
        // byte-identical coordinates. Joining against ch_plate_layout (by:
        // 0) also means a single-well sample -- which never gets a layout,
        // per the filter above -- is naturally excluded here too, with no
        // separate well-count filter needed on this side.
        ch_annotations
            .map { well_meta, geojson -> [ plate_meta_of(well_meta), [ well: well_meta.well, geojson: geojson.name ], geojson ] }
            .groupTuple(by: 0)
            .map { plate_meta, well_rows, geojsons ->
                def pairs = [ well_rows, geojsons ].transpose().sort { a, b -> a[0].well <=> b[0].well }
                [ plate_meta, pairs.collect { pair -> pair[0] }, pairs.collect { pair -> pair[1] } ]
            }
            .join(ch_plate_layout, by: 0)
            .set { ch_plate_geojson_input }

        AVITIMERGEPLATEGEOJSON(
            ch_plate_geojson_input
        )
        ch_versions          = ch_versions.mix(AVITIMERGEPLATEGEOJSON.out.versions.first())
        ch_plate_annotations = AVITIMERGEPLATEGEOJSON.out.annotations
    }

    //
    // Optional SEGMENTATIONREPORT module. AVITI's channel names are fixed
    // ("Nucleus", "Cell-Membrane"[, "Actin"]), unlike the free-form
    // per-sample channel mapping the COMET/MIBI samplesheet supplies.
    // run_cellpose is reported as true since both AVITI segmenters are
    // Cellpose-family models; run_mesmer/run_cellsam are false.
    //
    // Report scope is per sample, based on the *actual* discovered well
    // count (ch_well_count) rather than whether the samplesheet row
    // restricted `wells`: a sample with only one well never gets plate
    // assembly (see the >1-well filter above), so it must fall back to a
    // per-well report regardless of what `wells` said -- an unrestricted
    // row can still resolve to one well. A sample with more than one well
    // gets one plate-level report instead, IF plate assembly/reporting is
    // enabled -- otherwise every sample falls back to per-well reports, so
    // disabling the plate path never silently drops a sample's report.
    //
    ch_report = channel.empty()
    if (params.generate_report) {
        def plate_reports_enabled = params.aviti_plate_assembly && params.aviti_plate_report
        ch_report_scope = ch_well_count
            .map { sample, well_count ->
                def scope = (plate_reports_enabled && well_count > 1) ? 'plate' : 'well'
                [ sample, scope ]
            }

        ch_well_report_input = AVITISTITCHWELL.out.image
            .join(ch_annotations, by: 0)
            .map { meta, image, annotations -> [ meta.sample, meta, image, annotations ] }
            .combine(ch_report_scope, by: 0)
            .filter { _sample, _meta, _image, _annotations, scope -> scope == 'well' }
            .map { _sample, meta, image, annotations, _scope -> to_report_input(meta, annotations, image) }

        // Fed the full-resolution plate image, not the overview: spatialVis
        // crops its background with terra::crop(img, ext(<bbox from
        // GeoJSON>)), so image and GeoJSON must share one pixel frame, and
        // terra/GDAL open a pyramidal TIFF lazily and read only the small QC
        // windows it actually plots.
        ch_plate_report_input = channel.empty()
        if (plate_reports_enabled) {
            ch_plate_report_input = ch_plate_image
                .join(ch_plate_annotations, by: 0)
                .map { meta, image, annotations -> [ meta.sample, meta, image, annotations ] }
                .combine(ch_report_scope, by: 0)
                .filter { _sample, _meta, _image, _annotations, scope -> scope == 'plate' }
                .map { _sample, meta, image, annotations, _scope -> to_report_input(meta, annotations, image) }
        }

        SEGMENTATIONREPORT(
            ch_well_report_input.mix(ch_plate_report_input)
        )
        ch_versions = ch_versions.mix(SEGMENTATIONREPORT.out.versions.first())
        ch_report = SEGMENTATIONREPORT.out.report
    }

    emit:
    nuclear_segmentation_mask   = AVITISTITCHWELL.out.nuclear_mask // channel: [ val(meta), *.tif ] -- stitched instance-labeled nuclear mask
    wholecell_segmentation_mask = AVITISTITCHWELL.out.cell_mask    // channel: [ val(meta), *.tif ]
    annotations                 = ch_annotations                   // channel: [ val(meta), *.geojson ] -- well-scoped only, feeds the top-level KRONOS join
    kronos_input                = ch_kronos_input                  // channel: [ val(meta), tiff, whole_cell_mask ] -- well-scoped only, see above
    report                      = ch_report                        // channel: [ val(meta), *.html ] -- mix of per-well and plate-level reports
    plate_image                 = ch_plate_image                   // channel: [ val(plate_meta), *_plate.ome.tif ]
    plate_layout                = ch_plate_layout                  // channel: [ val(plate_meta), *_plate_layout.csv ]
    plate_annotations           = ch_plate_annotations             // channel: [ val(plate_meta), *_plate.geojson{,.gz} ]

    versions = ch_versions                                         // channel: [ versions.yml ]
}
