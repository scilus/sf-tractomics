process BUNDLE_IIT {
    tag 'atlas'
    label 'process_single'

    container "scilus/scilpy:2.2.0_cpu"

    input:
    path iit_bundles
    path thresholds_txt_file

    output:
    path "*.nii.gz",        emit: bundle_masks
    path "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p /tmp
    export HOME=/tmp

    # Move all input bundles to a dedicated folder to avoid
    # name conflicts.
    mkdir -p bundles_density_maps
    for f in \$(ls *.nii.gz); do
        mv "\$f" "bundles_density_maps/\$f"
    done

    # Create masks for each bundle based on provided thresholds
    while read -r name thr; do
        scil_volume_math lower_threshold "bundles_density_maps/\${name}.nii.gz" \$thr "\${name}.nii.gz" -f
    done < ${thresholds_txt_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    scil_volume_math -h

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
    END_VERSIONS
    """
}
