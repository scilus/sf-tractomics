include { HARMONIZATION_CLINICALCOMBAT  as CLINICALCOMBAT } from '../../../modules/nf-neuro/harmonization/clinicalcombat/main'
include { HARMONIZATION_STATS2COMBAT    as STATS2COMBAT_MOVING } from '../../../modules/nf-neuro/harmonization/stats2combat/main'
include { HARMONIZATION_STATS2COMBAT    as STATS2COMBAT_REFERENCE } from '../../../modules/nf-neuro/harmonization/stats2combat/main'
include { HARMONIZATION_COMBAT2STATS    as COMBAT2STATS } from '../../../modules/nf-neuro/harmonization/combat2stats/main'

workflow HARMONIZATION {

    take:
    ch_reference_site
    ch_moving_site

    main:
    ch_versions = channel.empty()

    if (!ch_moving_site || !ch_reference_site) {
        error "HARMONIZATION workflow requires both 'ch_moving_site' and 'ch_reference_site' inputs to be provided."
    }

    // Format the input stats files
    STATS2COMBAT_REFERENCE(ch_reference_site)
    ch_versions = ch_versions.mix(STATS2COMBAT_REFERENCE.out.versions)
    ch_reference_metrics = STATS2COMBAT_REFERENCE.out.stats_for_combat
        .flatten()
        .map{ file ->
            // Extract the metric name from the filename which
            // is in the format: <site>.<metric>.*
            [[metric: file.getName().split("\\.")[1]], file]
        }

    STATS2COMBAT_MOVING(ch_moving_site)
    ch_versions = ch_versions.mix(STATS2COMBAT_MOVING.out.versions)
    ch_moving_metrics = STATS2COMBAT_MOVING.out.stats_for_combat
        .flatten()
        .map{ file ->
            // Extract the metric name from the filename which
            // is in the format: <site>.<metric>.*
            [[metric: file.getName().split("\\.")[1]], file]
        }

    // Group moving and reference metrics by metric name
    ch_grouped_metrics = ch_moving_metrics
        .join(ch_reference_metrics)
        .map{ _meta, moving_entry, reference_entry ->
            [ reference_entry, moving_entry ]
        }

    // Run the harmonization
    CLINICALCOMBAT(ch_grouped_metrics)
    ch_versions = ch_versions.mix(CLINICALCOMBAT.out.versions.first())

    // Group by site
    ch_harmonized_files = CLINICALCOMBAT.out.harmonizedsite
        .map{ file ->
            // Extract the metric name from the filename which
            // is in the format: <site>.<metric>.*
            [[site: file.getName().split("\\.")[0]], file]
        }
        .groupTuple()
        .map { _site, files -> files }

    // Combine/format the output harmonized metrics into a MultiQC friendly TSV format
    COMBAT2STATS(ch_harmonized_files)
    ch_versions = ch_versions.mix(COMBAT2STATS.out.versions)

    emit:
    harmonized_stats    = COMBAT2STATS.out.stats_for_mqc
    figures             = CLINICALCOMBAT.out.figures
    model               = CLINICALCOMBAT.out.model
    qc_reports          = CLINICALCOMBAT.out.bdqc
    qc_plot_data_json   = CLINICALCOMBAT.out.plot_data_json
    versions            = ch_versions
}
