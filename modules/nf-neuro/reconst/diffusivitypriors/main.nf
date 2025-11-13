process RECONST_DIFFUSIVITYPRIORS {
    tag "$meta.id"
    label 'process_single'

    container "scilus/scilpy:2.2.0_cpu"

    input:
        tuple val(meta), path(fa), path(ad), path(rd), path(md)

    output:
        tuple val(meta), path("*_para_diff.txt")        , emit: para_diff
        tuple val(meta), path("*_iso_diff.txt")         , emit: iso_diff
        tuple val(meta), path("*_perp_diff.txt")        , emit: perp_diff, optional: true
        tuple val(meta), env("para_diff")               , emit: para_diff_val
        tuple val(meta), env("iso_diff")                , emit: iso_diff_val
        tuple val(meta), env("perp_diff")               , emit: perp_diff_val, optional: true
        path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def fa_min = task.ext.fa_min ? "--fa_min " + task.ext.fa_min : ""
    def fa_max = task.ext.fa_max ? "--fa_max " + task.ext.fa_max : ""
    def md_min = task.ext.md_min ? "--md_min " + task.ext.md_min : ""
    def roi_radius = task.ext.roi_radius ? "--roi_radius " + task.ext.roi_radius : ""

    """

    scil_NODDI_priors $fa $ad $rd $md $fa_min $fa_max $md_min $roi_radius \
        --out_txt_1fiber_para ${prefix}_para_diff.txt \
        --out_txt_1fiber_perp ${prefix}_perp_diff.txt \
        --out_txt_ventricles ${prefix}_iso_diff.txt

    # Set output environment variables
    para_diff=\$(cat ${prefix}_para_diff.txt)
    iso_diff=\$(cat ${prefix}_iso_diff.txt)
    if [[ -e ${prefix}_perp_diff.txt ]]
    then
        perp_diff=\$(cat ${prefix}_perp_diff.txt)
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    scil_NODDI_priors -h

    touch ${prefix}_para_diff.txt
    touch ${prefix}_iso_diff.txt
    touch ${prefix}_perp_diff.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
    END_VERSIONS
    """
}
