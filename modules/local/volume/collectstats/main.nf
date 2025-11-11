process VOLUME_COLLECTSTATS {
    tag "$meta.id"
    label 'process_single'

    container "scilus/scilus:2.2.0"

    input:
    tuple val(meta), path(csv_files)

    output:
    tuple val(meta), path("*_grouped*.csv"),  emit: csv
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Check that all CSV files have the same header
    header=\$(head -n 1 ${csv_files[0]})
    for f in ${csv_files}; do
        if [ "\$header" != "\$(head -n 1 \$f)" ]; then
            echo "Error: CSV files have different headers!"
            exit 1
        fi
    done
    head -n 1 ${csv_files[0]} > ${prefix}_grouped.csv

    # Now concatenate all CSV files, skipping their headers
    for f in ${csv_files}; do
        tail -n +2 \$f >> ${prefix}_grouped.csv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        volume: \$(volume --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_grouped.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        volume: \$(volume --version)
    END_VERSIONS
    """
}
