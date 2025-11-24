process STATS_METRICSINROI {
    tag "$meta.id"
    label 'process_single'

    container "scilus/scilpy:dev"

    input:
    tuple val(meta), path(metrics), path(rois), path(rois_lut)  /* optional, input = [] */

    output:
    tuple val(meta), path("*_stats.json")       , emit: stats_json
    tuple val(meta), path("*_stats.{csv,tsv}")  , emit: stats_tab
    path "versions.yml"                         , emit: versions

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

    # Get all bundles names from the JSON
    bundles=\$(jq -r "keys[]" ${prefix}_${suffix}.json)

    # All bundles have the same metrics. To get the metrics names from
    # the JSON, we can just fetch them from the first bundle.
    first_bundle=\$(printf '%s\\n' \$bundles | head -n 1)

    # Extract the metrics names from this first bundle
    metrics=\$(FIRST_BUNDLE="\$first_bundle" jq -r ".\\"\$first_bundle\\" | keys[]" ${prefix}_${suffix}.json)

    # Create the CSV/TSV header
    # (SID, bundle, metric, mean, std)
    echo "sid${sep}bundle${sep}metric${sep}mean${sep}std" > ${prefix}_${suffix}.${output_format}

    for bundle in \$bundles;
    do
        for metric in \$metrics;
        do
            # Fetch the "mean" and "std" values from each bundle/metric
            # pair from the JSON
            mean=\$(jq -r --arg BUNDLE "\$bundle" --arg METRIC "\$metric" '.[\$BUNDLE].[\$METRIC].mean' ${prefix}_${suffix}.json)
            std=\$(jq -r --arg BUNDLE "\$bundle" --arg METRIC "\$metric" '.[\$BUNDLE].[\$METRIC].std' ${prefix}_${suffix}.json)

            # Build the following line and append it to the CSV/TSV file
            # (SID, bundle, metric, mean, std)
            line="${prefix}${sep}\${bundle}${sep}\${metric}${sep}\${mean}${sep}\${std}"
            echo \$line >> ${prefix}_${suffix}.${output_format}
        done
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
    """
    scil_volume_stats_in_ROI -h
    scil_volume_stats_in_labels -h

    touch ${prefix}_${suffix}.json
    touch ${prefix}_${suffix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
        jq: \$(jq --version |& sed '1!d ; s/jq-//')
    END_VERSIONS
    """
}
