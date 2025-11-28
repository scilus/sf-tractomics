process STATS_METRICSINROI {
    tag "$meta.id"
    label 'process_single'

    container "scilus/scilpy:dev"

    input:
    tuple val(meta), path(metrics), path(rois), path(rois_lut)  /* optional, input = [] */

    output:
    tuple val(meta), path("*.json")                   , emit: stats_json
    tuple val(meta), path("*_desc-mean_*.{csv,tsv}")  , emit: stats_mean
    tuple val(meta), path("*_desc-std_*.{csv,tsv}")   , emit: stats_std
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.first_suffix ? "${task.ext.first_suffix}_stats" : "stats"
    def bin = task.ext.bin ? "--bin " : ""
    def normalize_weights = task.ext.normalize_weights ? "--normalize_weights " : ""
    def use_label = task.ext.use_label ? true : false
    def key_substrs_to_remove = task.ext.key_substrs_to_remove ?: []
    def value_substrs_to_remove = task.ext.value_substrs_to_remove ?: []
    def output_format = task.ext.output_format ?: 'tsv'  // 'csv' or 'tsv'

    assert output_format in ['csv', 'tsv'] : "output_format must be either 'csv' or 'tsv'"

    def sep = output_format == 'tsv' ? '\t' : ','
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1

    if $use_label;
    then
        if [[ ! -f "$rois_lut" ]];
        then
            echo "ROI LUT is missing. Will fail."
        fi

        scil_volume_stats_in_labels $rois $rois_lut \
            --metrics $metrics \
            --sort_keys > ${prefix}_${suffix}.json
    else
        scil_volume_stats_in_ROI $rois \
            --metrics $metrics \
            --sort_keys \
            --keep_unique_roi_name \
            $bin $normalize_weights > ${prefix}_${suffix}.json
    fi

    # Remove all substrings from the keys as specified
    # in the configuration via 'task.ext.key_substrs_to_remove'
    for substr in ${key_substrs_to_remove.join(' ')};
    do
        SUBSTR="\$substr" jq -r '
            with_entries(.key |= sub(env.SUBSTR; ""))
        ' ${prefix}_${suffix}.json > ${prefix}_${suffix}_tmp.json
        mv ${prefix}_${suffix}_tmp.json ${prefix}_${suffix}.json
    done

    # Remove all substrings from the values as specified
    # in the configuration via 'task.ext.value_substrs_to_remove'
    for substr in ${value_substrs_to_remove.join(' ')};
    do
        SUBSTR="\$substr" jq -r '
            with_entries(
                .value |= with_entries(
                    .key |= sub(env.SUBSTR; "")
                )
            )
        ' ${prefix}_${suffix}.json > ${prefix}_${suffix}_tmp.json
        mv ${prefix}_${suffix}_tmp.json ${prefix}_${suffix}.json
    done

    # Get all ROIs names from the JSON
    rois=\$(jq -r "keys[]" ${prefix}_${suffix}.json)

    # All ROIs have the same metrics. To get the metrics names from
    # the JSON, we can just fetch them from the first ROI.
    first_roi=\$(printf '%s\\n' \$rois | head -n 1)

    # Extract the metrics names from this first roi
    metrics=\$(FIRST_ROI="\$first_roi" jq -r ".\\"\$first_roi\\" | keys[]" ${prefix}_${suffix}.json)

    # Create the CSV/TSV headers
    # (sample, roi, metric1, metric2, ..., metricN)
    header_mean="sample${sep}roi"
    header_std="sample${sep}roi"
    for metric in \$metrics; do
        header_mean="\${header_mean}${sep}\${metric}"
        header_std="\${header_std}${sep}\${metric}"
    done
    echo "\$header_mean" > ${prefix}_desc-mean_${suffix}.${output_format}
    echo "\$header_std" > ${prefix}_desc-std_${suffix}.${output_format}

    for roi in \$rois;
    do
        # Initialize lines with sample and roi
        line_mean="${prefix}${sep}\${roi}"
        line_std="${prefix}${sep}\${roi}"

        for metric in \$metrics;
        do
            # Fetch the "mean" and "std" values from each roi/metric
            # pair from the JSON
            val_mean=\$(jq -r --arg ROI "\$roi" --arg METRIC "\$metric" '.[\$ROI].[\$METRIC].mean' ${prefix}_${suffix}.json)
            val_std=\$(jq -r --arg ROI "\$roi" --arg METRIC "\$metric" '.[\$ROI].[\$METRIC].std' ${prefix}_${suffix}.json)

            # Append values to the lines
            line_mean="\${line_mean}${sep}\${val_mean}"
            line_std="\${line_std}${sep}\${val_std}"
        done

        # Append the completed lines to the files
        echo "\$line_mean" >> ${prefix}_desc-mean_${suffix}.${output_format}
        echo "\$line_std" >> ${prefix}_desc-std_${suffix}.${output_format}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
        jq: \$(jq --version |& sed '1!d ; s/jq-//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.first_suffix ? "${task.ext.first_suffix}_stats" : "stats"
    def output_format = task.ext.output_format ?: 'tsv'  // 'csv' or 'tsv'
    assert output_format in ['csv', 'tsv'] : "output_format must be either 'csv' or 'tsv'"
    """
    scil_volume_stats_in_ROI -h
    scil_volume_stats_in_labels -h

    touch ${prefix}_${suffix}.json
    touch ${prefix}_desc-mean_${suffix}.${output_format}
    touch ${prefix}_desc-std_${suffix}.${output_format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
        jq: \$(jq --version |& sed '1!d ; s/jq-//')
    END_VERSIONS
    """
}
