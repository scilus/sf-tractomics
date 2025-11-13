/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nf-tractoflow_pipeline'
include { TRACTOFLOW             } from '../subworkflows/nf-neuro/tractoflow'
include { RECONST_SHSIGNAL       } from '../modules/nf-neuro/reconst/shsignal'
include { RECONST_NODDI          } from '../modules/nf-neuro/reconst/noddi/main'
include { RECONST_NODDI as NODDI_KERNELS } from '../modules/nf-neuro/reconst/noddi/main'
include { RECONST_DIFFUSIVITYPRIORS } from '../modules/nf-neuro/reconst/diffusivitypriors/main'
include { RECONST_MEANDIFFUSIVITYPRIORS } from '../modules/local/reconst/meandiffusivitypriors/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NF_TRACTOFLOW {
    take:
    ch_t1
    ch_wmparc
    ch_aparc_aseg
    ch_dwi_bval_bvec
    ch_b0
    ch_rev_dwi_bval_bvec
    ch_rev_b0
    ch_lesion
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_topup_config = Channel.empty()
    ch_bet_template = Channel.empty()
    ch_bet_probability = Channel.empty()

    /* Load topup config if provided */
    if ( params.config_topup ) {
        if ( file(params.config_topup).exists()) {
            ch_topup_config = Channel.fromPath(params.config_topup, checkIfExists: true)
        }
        else {
            ch_topup_config = Channel.from( params.config_topup )
        }
    }

    /* Load bet template */
    template_directory = file(params.template_t1 ?: "$projectDir/assets/templates/mni_152_sym_09c/t1")
    if (template_directory.exists() && template_directory.isDirectory()){
        ch_bet_template = ch_t1.map{ it[0] }
            .combine(Channel.fromPath(template_directory / "t1_template.nii.gz"))
        ch_bet_probability = ch_t1.map{ it[0] }
            .combine(Channel.fromPath(template_directory / "t1_brain_probability_map.nii.gz"))
    }
    else {
        error "A T1w template is required for brain extraction. Provide its directory with params.template_t1."
    }

    TRACTOFLOW(
        ch_dwi_bval_bvec,
        ch_t1,
        ch_b0
            .filter{ it[1] },
        ch_rev_dwi_bval_bvec
            .filter{ it[1] },
        ch_rev_b0
            .filter{ it[1] },
        ch_wmparc
            .filter{ it[1] },
        ch_aparc_aseg
            .filter{ it[1] },
        ch_topup_config,
        ch_bet_template,
        ch_bet_probability,
        ch_lesion
            .filter{ it[1] }
    )
    ch_versions = ch_versions.mix(TRACTOFLOW.out.versions)

    //
    // Run RECONST/SH_METRICS
    //
    if (params.sh_fitting)
        RECONST_SHSIGNAL(
            TRACTOFLOW.out.dwi
                .map{ it + [[]] }
        )

    //
    // Run RECONST/NODDI
    //
    if (params.run_noddi) {

        // Prepare NODDI inputs. This channel will be combined/joined in the
        // lines that follow with diffusivity priors w.r.t the following 3 scenarios:
        // Option 1: The user specifies the diffusivity priors to use (via params.para_diff and params.iso_diff).
        // Option 2: The user wants to compute the mean diffusivity priors across subjects. (Recommended)
        // Option 3: The user wants to compute diffusivity priors for each subject individually.

        ch_noddi_input = TRACTOFLOW.out.dwi
            .join(TRACTOFLOW.out.b0_mask)

        if ( (params.iso_diff != null && params.para_diff == null) ||
             (params.iso_diff == null && params.para_diff != null) ) {
            error "Please provide both iso_diff and para_diff parameters to use custom diffusivity priors."
        }

        if (params.iso_diff != null && params.para_diff != null) {
            if (params.average_diff_priors) {
                log.warn "Both custom diffusivity priors and average_diff_priors parameter were provided."
                    "The specified diffusivity priors (iso_diff: ${params.iso_diff}, "
                    "para_diff: ${params.para_diff}) will be used across subjects."
            }

            ch_noddi_input = ch_noddi_input
                .combine(Channel.value(params.para_diff))
                .combine(Channel.value(params.iso_diff))
        }
        else {
            // Compute diffusivity priors for each subject.
            RECONST_DIFFUSIVITYPRIORS(
                TRACTOFLOW.out.dti_fa
                    .join(TRACTOFLOW.out.dti_ad)
                    .join(TRACTOFLOW.out.dti_rd)
                    .join(TRACTOFLOW.out.dti_md)
            )
            ch_versions = ch_versions.mix(RECONST_DIFFUSIVITYPRIORS.out.versions.first())

            // Then compute mean diffusivity priors across subjects.
            if (params.average_diff_priors) {
                RECONST_MEANDIFFUSIVITYPRIORS(
                    RECONST_DIFFUSIVITYPRIORS.out.iso_diff
                        .map{ _meta, path -> path }
                        .collect(),
                    RECONST_DIFFUSIVITYPRIORS.out.para_diff
                        .map{ _meta, path -> path }
                        .collect(),
                    RECONST_DIFFUSIVITYPRIORS.out.perp_diff
                        .map{ _meta, path -> path }
                        .collect()
                )
                ch_versions = ch_versions.mix(RECONST_MEANDIFFUSIVITYPRIORS.out.versions.first())

                ch_noddi_input = ch_noddi_input
                    .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.mean_para_diff_val)
                    .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.mean_iso_diff_val)
            }
        }

        ch_noddi_input = ch_noddi_input
            .map{ meta, dwi, bval, bvec, b0_mask, para_diff, iso_diff ->
                [meta, dwi, bval, bvec, b0_mask, [], para_diff, iso_diff] }

        RECONST_NODDI( ch_noddi_input )
        ch_versions = ch_versions.mix(RECONST_NODDI.out.versions.first())
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nf-tractoflow_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
