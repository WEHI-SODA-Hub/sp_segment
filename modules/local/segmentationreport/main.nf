process SEGMENTATIONREPORT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/wehi-soda-hub/spatialvis:0.1.1'

    input:
    tuple val(meta),
          path(annotations),
          val(run_mesmer),
          val(run_cellpose),
          val(nuclear_channel),
          val(membrane_channels)

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
        output_dir = '.',
        overwrite = TRUE
    )"
    quarto render segmentation_report_template.qmd \\
        --to html \\
        --no-cache \\
        --output ${prefix}.html \\
        ${args} \\
        -P sample_name:${meta.id} \\
        -P geojson_file:${annotations} \\
        -P nuclear_channel:${nuclear_channel} \\
        -P membrane_channels:"${membrane_channels}" \\
        -P run_cellpose:${run_cellpose} \\
        -P cellpose_diameter:${params.cellpose_diameter} \\
        -P cellpose_min_area:${params.cellpose_min_area} \\
        -P cellpose_flow_threshold:${params.cellpose_flow_threshold} \\
        -P cellpose_cellprob_threshold:${params.cellpose_cellprob_threshold} \\
        -P cellpose_model_type:${params.cellpose_model_type} \\
        -P cellpose_pretrained_model:${params.cellpose_pretrained_model} \\
        -P run_mesmer:${run_mesmer} \\
        -P mesmer_segmentation_level:${params.mesmer_segmentation_level} \\
        -P mesmer_maxima_threshold:${params.mesmer_maxima_threshold} \\
        -P mesmer_interior_threshold:${params.mesmer_interior_threshold} \\
        -P mesmer_maxima_smooth:${params.mesmer_maxima_smooth} \\
        -P mesmer_min_nuclei_area:${params.mesmer_min_nuclei_area} \\
        -P mesmer_remove_border_cells:${params.mesmer_remove_border_cells} \\
        -P mesmer_pixel_expansion:${params.mesmer_pixel_expansion} \\
        -P mesmer_padding:${params.mesmer_padding}

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
