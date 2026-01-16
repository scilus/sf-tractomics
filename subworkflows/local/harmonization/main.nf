include { HARMONIZATION_CLINICALCOMBAT } from '../../../modules/local/harmonization/clinicalcombat/main'
include { HARMONIZATION_FORMATSTATS as FORMAT_INPUT_MOVING } from '../../../modules/local/harmonization/formatstats/main'
include { HARMONIZATION_FORMATSTATS as FORMAT_INPUT_REFERENCE } from '../../../modules/local/harmonization/formatstats/main'
include { HARMONIZATION_FORMATSTATS as FORMAT_OUTPUT } from '../../../modules/local/harmonization/formatstats/main'

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
    FORMAT_INPUT_REFERENCE(ch_reference_site)
    ch_versions = ch_versions.mix(FORMAT_INPUT_REFERENCE.out.versions)
    ch_reference_metrics = FORMAT_INPUT_REFERENCE.out.raw_files
        .flatten()
        .map{ file -> extract_metric_from_filename(file) }

    FORMAT_INPUT_MOVING(ch_moving_site)
    ch_versions = ch_versions.mix(FORMAT_INPUT_MOVING.out.versions)
    ch_moving_metrics = FORMAT_INPUT_MOVING.out.raw_files
        .flatten()
        .map{ file -> extract_metric_from_filename(file) }

    // Group moving and reference metrics by metric name
    ch_grouped_metrics = ch_moving_metrics.join(ch_reference_metrics)
        .map{ _meta, moving_entry, reference_entry ->
            [ reference_entry, moving_entry ]
        }

    // Run the harmonization
    HARMONIZATION_CLINICALCOMBAT(ch_grouped_metrics)
    ch_versions = ch_versions.mix(HARMONIZATION_CLINICALCOMBAT.out.versions.first())

    // Group by site
    ch_harmonized_files = HARMONIZATION_CLINICALCOMBAT.out.harmonizedsite
        .map{ file -> extract_site_from_filename(file) }
        .groupTuple()
        .map { _site, files -> files }

    // Combine/format the output harmonized metrics into a MultiQC friendly TSV format
    FORMAT_OUTPUT(ch_harmonized_files)
    ch_versions = ch_versions.mix(FORMAT_OUTPUT.out.versions)

    emit:
    harmonized_metrics   = FORMAT_OUTPUT.out.harmonized_files
    figures              = HARMONIZATION_CLINICALCOMBAT.out.figures
    model                = HARMONIZATION_CLINICALCOMBAT.out.model
    qc_reports           = HARMONIZATION_CLINICALCOMBAT.out.bdqc

    versions = ch_versions                     // channel: [ versions.yml ]
}

def extract_metric_from_filename(file) {
    // The file will be called something like:
    // site.metric.raw.csv
    def filename = file.getName()
    def parts = filename.split("\\.")

    if (parts.length < 4) {
        throw new IllegalArgumentException(
            "Filename does not conform to expected pattern (site.metric.*.csv): " + filename)
    }
    def metric = parts[1]

    return [[metric: metric], file]
}

def extract_site_from_filename(file) {
    // The file will be called something like:
    // site.metric.raw.csv
    def filename = file.getName()
    def parts = filename.split("\\.")

    if (parts.length < 4) {
        throw new IllegalArgumentException(
            "Filename does not conform to expected pattern (site.metric.*.csv): " + filename)
    }
    def site = parts[0]

    return [[site: site], file]
}
