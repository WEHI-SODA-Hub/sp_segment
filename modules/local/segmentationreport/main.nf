process SEGMENTATIONREPORT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/wehi-soda-hub/spatialvis:0.1.0'

    input:
    tuple val(meta), path(annotations), val(run_mesmer), val(run_cellpose)

    output:
    tuple val(meta), path("*/*.html"), path("*/*_files/*"), emit: report
    tuple val(meta), path("*/*.rds"), emit: rds, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    Rscript -e "spatialVis::copy_report_template(
        template_name = 'segmentation_report_template.qmd',
        output_dir = '.'
    )"
    quarto render segmentation_report_template.qmd \\
        --to html \\
        --no-cache \\
        --output ${prefix}.html \\
        ${args} \\
        -P sample_name:${meta.id} \\
        -P geojson_file:${annotations}

    mkdir -p ${prefix}
    mv ${prefix}.html ${prefix}
    if [[ -f ${prefix}.rds ]]; then
        mv ${prefix}.rds ${prefix}
    fi
    cp -r segmentation_report_template_files ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        spatialVis: \$(Rscript -e "cat(as.character(packageVersion('spatialVis')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/segmentation_report_template_files
    touch ${prefix}/${prefix}.html
    touch ${prefix}/${prefix}.rds
    touch ${prefix}/segmentation_report_template_files/foo.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        spatialVis: \$(Rscript -e "cat(as.character(packageVersion('spatialVis')))")
    END_VERSIONS
    """
}
