process VOLUME_ROISTATS {
    tag "$meta.id"
    label 'process_high'

    container "scilus/scilpy:2.2.0_cpu"

    input:
        tuple val(meta), path(metrics_list), path(masks_list)
    output:
        tuple val(meta), path("*__stats.json"), emit: stats_json
        tuple val(meta), path("*__stats.csv"), emit: stats_csv
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    shopt -s extglob
    mkdir metrics


    # FIXME: The way we extract the metric names is fragile.
    for metric in $metrics_list;
    do
        pos=\$((\$(echo \$metric | grep -b -o __ | cut -d: -f1)+2))
        bname=\${metric:\$pos}
        bname=\$(basename \$bname .nii.gz)
        mv \$metric metrics/\${bname}.nii.gz
    done

    # FIXME: The way we extract the masks/bundles names is fragile.
    new_mask_list=""
    for mask in $masks_list;
    do
        bname=\$(echo \$mask | sed -E 's/.*__(.*)_mask__.*/\\1/')
        mv \$mask \${bname}.nii.gz
        new_mask_list="\$new_mask_list \${bname}.nii.gz"
    done

    scil_volume_stats_in_ROI \$new_mask_list --metrics_dir metrics -f > ${prefix}__stats.json

    # Convert the JSON to CSV
    bundles=\$(jq -r "keys[]" ${prefix}__stats.json)

    # Need the first bundle name to extract the metrics names from it
    first_bundle=\$(printf '%s\\n' \$bundles | head -n 1)

    # Extract the metrics names from this first bundle
    metrics=\$(FIRST_BUNDLE="\$first_bundle" jq -r ".\\"\$first_bundle\\" | keys[]" ${prefix}__stats.json)

    # Output the CSV file
    # echo "sid,bundle,\$(echo \${metrics} | sed -z 's/\\n/,/g')" > ${prefix}__stats.csv

    echo "sid,bundle,metric,mean,std" > ${prefix}__stats.csv
    for bundle in \$bundles;
    do
        for metric in \$metrics;
        do
            mean=\$(jq -r --arg BUNDLE "\$bundle" --arg METRIC "\$metric" '.[\$BUNDLE].[\$METRIC].mean' ${prefix}__stats.json)
            std=\$(jq -r --arg BUNDLE "\$bundle" --arg METRIC "\$metric" '.[\$BUNDLE].[\$METRIC].std' ${prefix}__stats.json)

            line="${prefix},\${bundle},\${metric},\${mean},\${std}"
            echo \$line >> ${prefix}__stats.csv
        done
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    scil_volume_stats_in_ROI -h

    touch ${prefix}__stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
    END_VERSIONS
    """
}
