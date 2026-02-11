process HARMONIZATION_CLINICALCOMBAT {
    label 'process_medium'

    container "scilus/clinical_combat:1.1.2"

    input:
    tuple path(ref_site), path(move_site)

    output:
    path("*.model.csv")                , emit: model
    path("*.harmonized.csv.gz")        , emit: harmonizedsite
    path("*bhattacharrya.txt")         , emit: bdqc
    path("*.png")                      , emit: figures
    path("*.json")                     , emit: plot_data_json
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def method = task.ext.method ? "--method " + task.ext.method : ""
    def bundles_list = task.ext.bundles ? "--bundles " + task.ext.bundles : "--bundles all"
    def regul_ref = task.ext.regul_ref ? "--regul_ref " + task.ext.regul_ref : ""
    def regul_mov = task.ext.regul_mov ? "--regul_mov " + task.ext.regul_mov : ""
    def degree = task.ext.degree ? "--degree " + task.ext.degree : ""
    def nu = task.ext.nu ? "--nu " + task.ext.nu : ""
    def tau = task.ext.tau ? "--tau " + task.ext.tau : ""
    def degree_qc = task.ext.degree_qc ? "--degree_qc " + task.ext.degree_qc : ""

    def limit_age = task.ext.limit_age_range ? "--limit_age_range " : ""
    def ignore_sex = task.ext.ignore_sex ? "--ignore_sex " : ""
    def ignore_handedness = task.ext.ignore_handedness ? "--ignore_handedness " : ""
    def no_eb = task.ext.no_empiral_bayes ? "--no_empiral_bayes " : ""

    """
    combat_pipeline $ref_site $move_site $method $bundles_list \
        $limit_age $ignore_sex $ignore_handedness \
        $regul_ref $regul_mov $degree $nu $tau $degree_qc \
        $no_eb \
        --save_curves_json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clinical-ComBAT: \$(uv pip -q -n list | grep clinical-ComBAT | tr -s ' ' | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def method = task.ext.method ?: "clinical"

    // Extract the site and metric name from the input filenames
    // We do this to avoid having file collisions when stubbing
    def ref_site_name = ref_site.getName().split("\\.")[0]
    def mov_site_name = move_site.getName().split("\\.")[0]
    def metric_name = ref_site.getName().split("\\.")[1]

    """
    combat_pipeline -h

    touch dummy_figure_1.png
    touch dummy_figure_2.png
    touch ${mov_site_name}.${metric_name}.raw.bhattacharrya.txt
    touch ${mov_site_name}.${metric_name}.${method}.harmonized.bhattacharrya.txt
    touch ${ref_site_name}_${mov_site_name}.${metric_name}.${method}.model.csv
    touch ${ref_site_name}_${mov_site_name}.${metric_name}.${method}.harmonized.csv.gz
    touch DataModels_${ref_site_name}-${mov_site_name}_${method}_${metric_name}.json
    touch AgeCurve_${ref_site_name}-${mov_site_name}_${method}_${metric_name}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clinical-ComBAT: \$(uv pip -q -n list | grep clinical-ComBAT | tr -s ' ' | cut -d' ' -f2)
    END_VERSIONS
    """
}
