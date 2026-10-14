//
// Subworkflow with functionality specific to the WEHI-SODA-Hub/sp_segment pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version            // boolean: Display version and exit
    validate_params    // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args  //   array: List of positional nextflow CLI args
    outdir             //  string: The output directory where the results will be saved
    _input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Fail before an hour is spent pulling a container into a full home quota.
    //
    // An Apptainer pull writes two things. The finished .img goes to
    // apptainer.cacheDir, which --container_cache_dir sets. The OCI layer blobs
    // it is converted from go to $APPTAINER_CACHEDIR -- ~/.apptainer/cache by
    // default, and several GB per image, which is what actually exhausts a home
    // quota. Nextflow runs the pull as a plain subprocess of its own JVM, so
    // that variable comes from the shell that launched the pipeline and no
    // amount of nextflow.config can redirect it. Checking it is all we can do.
    //
    // Only the finished image is retried on failure, so this surfaces as
    // "Failed to pull singularity image ... disk quota exceeded" after the full
    // download has already run.
    //
    if (workflow.containerEngine in ['apptainer', 'singularity']) {
        def home_dir   = System.getenv('HOME')
        def blob_cache = System.getenv('APPTAINER_CACHEDIR') ?: System.getenv('SINGULARITY_CACHEDIR')
        def under_home = home_dir && (
            !blob_cache || file(blob_cache).toAbsolutePath().normalize().startsWith(file(home_dir).toAbsolutePath().normalize())
        )
        def suggested = params.container_cache_dir
            ? file(params.container_cache_dir).toAbsolutePath().normalize().parent.resolve('apptainer_cache')
            : '/path/to/scratch/apptainer_cache'
        def advice = (
            "Apptainer unpacks OCI layers into \$APPTAINER_CACHEDIR (currently " +
            "${blob_cache ?: home_dir + '/.apptainer/cache'}), which is tens of GB for this pipeline's images. " +
            "Export it to the same filesystem as --container_cache_dir before launching:\n" +
            "    export APPTAINER_CACHEDIR=${suggested}\n" +
            "    export SINGULARITY_CACHEDIR=\$APPTAINER_CACHEDIR\n" +
            "This cannot be set from nextflow.config: Nextflow pulls images from its own process, " +
            "not from a task."
        )
        if (params.container_cache_dir && under_home) {
            // Opt-in: setting container_cache_dir says the home directory
            // cannot hold these images, and half of the pull would still land
            // there. Erroring beats a partial fix that fails identically.
            error("container_cache_dir is set, but the container layer cache is still under \$HOME.\n" + advice)
        }
        else if (under_home && !System.getenv('CI')) {
            // Not on CI. A GitHub runner's $HOME has no quota and tens of GB
            // free, so this fires on every run there and says nothing useful --
            // and a warning that is usually noise is one nobody reads when it
            // is not. CI=true is set by GitHub Actions, GitLab CI and the rest.
            log.warn("Container images will be cached under \$HOME. " + advice)
        }
    }

    //
    // Fail before any compute on a model that cannot be resolved.
    //
    // Built-in cellpose models, as reported by `cellpose.models.MODEL_NAMES` in
    // the pinned container (cellpose 4.2.1.1). Anything else given to
    // cellpose_pretrained_model is treated as a path to a custom model.
    //
    // cellpose_pretrained_model takes either a built-in model name or a path to
    // a custom model, so only validate as a path when it is not a known name.
    // Cellpose does not error on a --pretrained-model path that does not exist:
    // it falls back to its built-in weights, so a typo produced a full,
    // plausible, silently-wrong run.
    //
    def CELLPOSE_BUILTIN_MODELS = ['cpsam_v2', 'cpsam', 'cpdino', 'cpdino-vitb']
    if (params.cellpose_pretrained_model
        && !CELLPOSE_BUILTIN_MODELS.contains(params.cellpose_pretrained_model)
        && !file(params.cellpose_pretrained_model).exists()) {
        error(
            "cellpose_pretrained_model is neither a built-in model name nor an existing path: " +
            "${params.cellpose_pretrained_model}\n" +
            "Built-in models: ${CELLPOSE_BUILTIN_MODELS.join(', ')}"
        )
    }
    if (params.cellpose_models_dir && !file(params.cellpose_models_dir).exists()) {
        error("cellpose_models_dir does not exist: ${params.cellpose_models_dir}")
    }

    //
    // A model cache under /opt breaks the container, confusingly.
    //
    // Nextflow bind-mounts the host directory holding a staged input, so
    // --cellpose_models_dir /opt/... mounts the host's /opt over the
    // container's. The Cellpose image installs into /opt/conda, so that
    // shadows the interpreter and every task dies with
    // "sopa: command not found" -- an error that says nothing about the real
    // cause. /opt is a natural place to put a shared cache, so catch it here.
    //
    if (params.cellpose_models_dir && file(params.cellpose_models_dir).toAbsolutePath().toString().startsWith('/opt')) {
        error(
            "cellpose_models_dir must not be under /opt: ${params.cellpose_models_dir}\n" +
            "Nextflow bind-mounts the host directory containing a staged input, which would " +
            "mount the host's /opt over the container's and hide /opt/conda, where the Cellpose " +
            "image is installed. Tasks would fail with \"sopa: command not found\".\n" +
            "Put the cache somewhere else, e.g. /shared/cellpose_models."
        )
    }

    //
    // cellpose_model_type is retired rather than quietly ignored.
    //
    // Cellpose >=4.0.1 logs "model_type argument is not used in v4.0.1+" and
    // discards it, and sopa only reads it on the cellpose<4 path, so the old
    // 'cyto3' default silently produced cpsam results on every run. Failing is
    // better than repeating that.
    //
    if (params.cellpose_model_type) {
        error(
            "cellpose_model_type is no longer supported: cellpose 4 ignores --model-type.\n" +
            "Use --cellpose_pretrained_model instead " +
            "(${CELLPOSE_BUILTIN_MODELS.join(', ')}, or a path to a custom model)."
        )
    }

    //
    // deepcell_cache_dir is bind-mounted into the Mesmer and CellSAM
    // containers, and a bind source that does not exist stops the container
    // from starting ("mount source ... doesn't exist"). The modules create
    // their subdirectories inside the container, which is too late, so create
    // the directory here.
    //
    if (params.deepcell_cache_dir) {
        def cache = file(params.deepcell_cache_dir)
        if (!cache.exists() && !cache.mkdirs()) {
            error("Could not create deepcell_cache_dir: ${params.deepcell_cache_dir}")
        }
    }

    //
    // KRONOS2 inputs are checked up front so a typo fails here rather than
    // after segmentation and measurement have run.
    //
    if (params.enable_kronos) {
        if (!params.kronos_model_path) {
            error("enable_kronos is true but kronos_model_path is not set. " +
                "Set --kronos_model_path to a full KRONOS2 snapshot directory.")
        }
        def mapping = params.kronos_marker_mapping?.toString()?.trim()
        if (mapping && !mapping.startsWith('{') && !file(mapping).exists()) {
            error("kronos_marker_mapping is neither inline JSON nor an existing file: ${params.kronos_marker_mapping}")
        }
    }

    //
    // Create channel from input file provided through params.input
    //
    def samplesheet_rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
    validateInputSamplesheet(samplesheet_rows)

    Channel
        .fromList(samplesheet_rows)
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications


    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
// Each row is [ meta, run_backsub, run_mesmer, run_cellpose, run_cellsam,
// tiff, nuclear_channel, membrane_channels ]. The schema checks each value on
// its own; these rules span columns or rows, so they are checked here, before
// any task starts. Every problem is reported at once.
//
def validateInputSamplesheet(rows) {
    def errors = []

    // Outputs are published as <sample>.<ext> and channels are joined on the
    // sample, so a repeated name overwrites results and drops samples.
    rows.collect { row -> row[0].id }
        .countBy { id -> id }
        .findAll { _id, count -> count > 1 }
        .each { id, count ->
            errors << "Sample '${id}' appears ${count} times. Each row needs a unique sample name; " +
                "to run several segmenters on one image, use one row per segmenter with different sample names."
        }

    rows.each { row ->
        def (meta, run_backsub, run_mesmer, run_cellpose, run_cellsam, _tiff, nuclear_channel, membrane_channels) = row
        def segmenters = [ mesmer: run_mesmer, cellpose: run_cellpose, cellsam: run_cellsam ]
            .findAll { _name, enabled -> enabled }
            .keySet()

        if (segmenters.size() > 1) {
            errors << "Sample '${meta.id}' enables more than one segmenter (${segmenters.join(', ')}). " +
                "Set at most one of run_mesmer, run_cellpose and run_cellsam per row."
        }
        if (!run_backsub && segmenters.isEmpty()) {
            errors << "Sample '${meta.id}' has nothing to run: run_backsub, run_mesmer, run_cellpose and run_cellsam are all false."
        }
        if (!segmenters.isEmpty()) {
            if (!nuclear_channel) {
                errors << "Sample '${meta.id}' runs ${segmenters.join(', ')} but has no nuclear_channel."
            }
            if (!membrane_channels) {
                errors << "Sample '${meta.id}' runs ${segmenters.join(', ')} but has no membrane_channels."
            }
        }
    }

    if (errors) {
        error("Invalid samplesheet ${params.input}:\n  - " + errors.join("\n  - "))
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",


            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [


        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

