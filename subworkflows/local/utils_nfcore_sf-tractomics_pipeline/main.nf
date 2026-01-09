//
// Subworkflow with functionality specific to the scilus/sf-tractomics pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { IO_BIDS                   } from '../../nf-neuro/io_bids/main'
include { IO_SAFECASTINPUTS         } from '../../../modules/local/io/safecastinputs'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()
    ch_samplesheet = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        false   // Reinstate when/if we use conda/mamba :
                // workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    before_text = """
    -\033[2m------------------------------------------------------\033[0m-
                                                \033[0;32m _.--'"'.\033[0m
    \033[0;34m        ___          ___       __   __      \033[0;32m(  ( (   )\033[0m
    \033[0;34m  |\\ | |__  __ |\\ | |__  |  | |__) /  \\     \033[0;33m(o)_    ) )\033[0m
    \033[0;34m  | \\| |       | \\| |___ |__| |  \\ \\__/     \033[0;32m    (o)_.'\033[0m
                                                \033[0;32m     )/\033[0m
    \033[0;35m  scilus/sf-tractomics ${workflow.manifest.version}\033[0m
    -\033[2m------------------------------------------------------\033[0m-
    """
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { "    https://doi.org/${it.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
    * The nf-neuro project
        https://scilus.github.io/nf-neuro

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/scilus/sf-tractomics/blob/master/CITATIONS.md
    """

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    if (params.input) {
        //
        // params.input is either a BIDS compliant directory or a samplesheet
        //   - if directory, we assume it is BIDS
        //   - everything else is a samplesheet
        //
        if (file(params.input).isDirectory()) {
            IO_BIDS(
                channel.fromPath(params.input),
                channel.value(params.fsbids ?: []),
                channel.value(params.bidsignore ?: [])
            )
            ch_samplesheet = [
                t1: IO_BIDS.out.ch_t1,
                wmparc: IO_BIDS.out.ch_wmparc,
                aparc_aseg: IO_BIDS.out.ch_aparc_aseg,
                dwi_bval_bvec: IO_BIDS.out.ch_dwi_bval_bvec,
                b0: channel.empty(),
                rev_dwi_bval_bvec: IO_BIDS.out.ch_rev_dwi_bval_bvec,
                rev_b0: IO_BIDS.out.ch_rev_b0,
                lesion: channel.empty()
            ]

            if (params.participants_tsv) {
                participants_tsv_path = "${params.participants_tsv}"
            }
            else if (file("${params.input}/participants.tsv").exists()) {
                participants_tsv_path = "${params.input}/participants.tsv"
            }
            else {
                participants_tsv_path = null
                log.warn("No participants.tsv provided, covariates will not be used.")
            }
        }
        else {
            ch_input_sheets = Channel
                .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
                .map{
                    meta, dwi, bval, bvec, sbref, rev_dwi, rev_bval, rev_bvec, rev_sbref, t1, wmparc, aparc_aseg, lesion ->
                        return [
                            meta,
                            dwi,
                            bval,
                            bvec,
                            sbref ?: [],
                            rev_dwi ?: [],
                            rev_bval ?: [],
                            rev_bvec ?: [],
                            rev_sbref ?: [],
                            t1,
                            wmparc ?: [],
                            aparc_aseg ?: [],
                            lesion ?: []
                        ]
                }
                .map{ samplesheet ->
                    validateInputSamplesheet(samplesheet)
                }

            IO_SAFECASTINPUTS(ch_input_sheets)
            ch_samplesheet = IO_SAFECASTINPUTS.out.safe_inputs
                .multiMap{ meta, dwi, bval, bvec, sbref, rev_dwi, rev_bval, rev_bvec, rev_sbref, t1, wmparc, aparc_aseg, lesion ->
                    t1: [meta, t1]
                    wmparc: [meta, wmparc]
                    aparc_aseg: [meta, aparc_aseg]
                    dwi_bval_bvec: [meta, dwi, bval, bvec]
                    b0: [meta, sbref]
                    rev_dwi_bval_bvec: [meta, rev_dwi, rev_bval, rev_bvec]
                    rev_b0: [meta, rev_sbref]
                    lesion: [meta, lesion]
                }

            if (params.participants_tsv) {
                participants_tsv_path = params.participants_tsv
            } else {
                participants_tsv_path = null
                log.warn("No participants.tsv provided, covariates will not be used.")
            }
        }
    }

    // We avoid merging the covariates (i.e. the extra meta fields)
    // directly into the samplesheet's multimap meta fields, as those covariates are not used
    // in most of the pipeline steps. This means, that if the participants.tsv changes for whatever
    // reason, the entire cache of the pipeline would be invalidated, thus causing the
    // pipeline to reprocess everything from scratch. Instead, we provide the mergeCovariatesIntoMeta
    // function, which can be used to merge the covariates into the samplesheet's multimap
    // on the fly, when needed (which should be done only when the inputs requires those fields).
    ch_covariates = parseParticipantsTsv(participants_tsv_path, ch_samplesheet.t1)

    def participants_to_include = []
    def participants_to_exclude = []

    if (params.participant_label) {
        participants_to_include = params.participant_label.split(",").collect { item -> item.trim() }
        log.info "Including participants: ${participants_to_include.join(", ")}"
    }

    if (params.exclude_participant_label) {
        participants_to_exclude = params.exclude_participant_label.split(",").collect { item -> item.trim() }
        log.info "Excluding participants: ${participants_to_exclude.join(", ")}"
    }

    if (participants_to_include || participants_to_exclude) {
        ch_samplesheet = ch_samplesheet.collectEntries { key, value ->
            def filtered = value.filter { item ->
                def meta = item[0]

                // The user can provide a list of subjects to include or to exclude, separated by commas
                // in the same format as the prefix (i.e. sub-XX_ses-XX_run-XX). The implementation
                // allows the user to provide a list of subjects of different scopes to run.
                // For example, if the user provides:
                // sub-01, then all sessions and runs of sub-01 will be processed.
                // sub-01_ses-01, then all runs of sub-01_ses-01 will be processed.
                // sub-01_ses-01_run-01, then only sub-01_ses-01_run-01 will be processed.

                def sid = [meta.id].findAll { x -> x }.join("_")
                def sid_ses = [meta.id, meta.session].findAll { x -> x }.join("_")
                def sid_run = [meta.id, meta.session, meta.run].findAll { x -> x }.join("_")

                def is_included = sid in participants_to_include ||
                    sid_ses in participants_to_include ||
                    sid_run in participants_to_include

                def is_excluded = sid in participants_to_exclude ||
                    sid_ses in participants_to_exclude ||
                    sid_run in participants_to_exclude

                // If inclusion list is provided, only keep the participant if in that list
                // and not in the exclusion list.
                if (participants_to_include) {
                    return is_included && !is_excluded
                }
                // If we only have an exclusion list, filter out participants in that list
                return !is_excluded
            }

            return [key, filtered]
        }
    }

    emit:
    t1 = ch_samplesheet.t1
    wmparc = ch_samplesheet.wmparc
    aparc_aseg = ch_samplesheet.aparc_aseg
    dwi_bval_bvec = ch_samplesheet.dwi_bval_bvec
    b0 = ch_samplesheet.b0
    rev_dwi_bval_bvec = ch_samplesheet.rev_dwi_bval_bvec
    rev_b0 = ch_samplesheet.rev_b0
    lesion = ch_samplesheet.lesion
    covariates  = ch_covariates

    versions    = ch_versions
}


def parseParticipantsTsv(participants_path, ch_with_proper_meta) {

    if (participants_path == null) {
        return channel.empty()
    }

    // Define the schema for participants.tsv
    def tsv_file = file(participants_path)
    def all_tsv_headers = tsv_file.readLines()[0].split('\t')
        .collect { item -> item.trim() }
        .toList()

    // Create joining keys
    def primary_keys = ['participant_id', 'session', 'run']
    def content_keys = all_tsv_headers - primary_keys
    def default_content = content_keys.collectEntries { key -> [key, ""] }

    // Parse "${params.inputs}/participants.tsv"
    def ch_participants = channel.fromPath(tsv_file)
    def participants_content = ch_participants
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def id = row.participant_id
            def ses = row.session ? "ses-" + row.session: ""
            def run = row.run ? "run-" + row.run: ""

            def key = [id: id, session: ses, run: run]
            def content = default_content.clone()
            content_keys.each { ckey -> content[ckey] = row[ckey] }
            content = content.collectEntries { k, v -> [k.toLowerCase(), v] }
            return [key, content]
        }

    // Prepare keys
    def ch_original_meta = ch_with_proper_meta
        .map { meta, _content ->
            def key = [id: meta.id, session: meta.session ?: "", run: meta.run ?: ""]
            return [key, meta]
        }

    // Join with participants.tsv content
    def ch_covariates = ch_original_meta
        .join(participants_content, by: 0, remainder: true)
        .filter { _key, original_meta, _tsv_meta -> original_meta != null } // Remove unmatched entries from the participants.tsv
        .map { _key, original_meta, tsv_meta ->
            def extra_meta = tsv_meta ?: default_content.collectEntries { k, v -> [k.toLowerCase(), v] }
            return [original_meta, extra_meta]
        }

    return ch_covariates
}

def mergeCovariatesIntoMeta(ch_src, ch_covariates) {
    if (params.participants_tsv == null && !file(params.input).isDirectory()) {
        // The input is a samplesheet and no participants.tsv was provided.
        // So there are no covariates to parse.
        return ch_src
    }

    def ch = ch_src.join(ch_covariates, by: 0)
        .map { item ->
            def original_meta = item[0]
            def content = item[1..-2]
            def covariates = item[-1]
            // Merge original meta with covariates without overwriting existing fields
            def merged_meta = original_meta.clone()
            covariates.each { k, v ->
                if (merged_meta[k] == null || merged_meta[k] == "") {
                    merged_meta[k] = v
                }
            }
            return [merged_meta] + content
        }
    return ch
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    return input
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
