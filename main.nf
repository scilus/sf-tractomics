#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scilus/sf-tractomics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/scilus/sf-tractomics
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SF_TRACTOMICS  } from './workflows/sf-tractomics'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_sf-tractomics_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_sf-tractomics_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SCILUS_SF_TRACTOMICS {

    take:
    t1                  // channel: t1 read in from --input
    wmparc              // channel: wmparc read in from --input
    aparc_aseg          // channel: aparc_aseg read in from --input
    dwi_bval_bvec       // channel: dwi_bval_bvec read in from --input
    b0                  // channel: b0 read in from --input
    rev_dwi_bval_bvec   // channel: rev_dwi_bval_bvec read in from --input
    rev_b0              // channel: rev_b0 read in from --input
    lesion              // channel: lesion read in from --input
    covariates          // channel: covariates parsed from participants.tsv

    main:

    //
    // WORKFLOW: Run pipeline
    //
    SF_TRACTOMICS (
        t1,
        wmparc,
        aparc_aseg,
        dwi_bval_bvec,
        b0,
        rev_dwi_bval_bvec,
        rev_b0,
        lesion,
        covariates
    )
    emit:
    multiqc_report = SF_TRACTOMICS.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    SCILUS_SF_TRACTOMICS (
        PIPELINE_INITIALISATION.out.t1,
        PIPELINE_INITIALISATION.out.wmparc,
        PIPELINE_INITIALISATION.out.aparc_aseg,
        PIPELINE_INITIALISATION.out.dwi_bval_bvec,
        PIPELINE_INITIALISATION.out.b0,
        PIPELINE_INITIALISATION.out.rev_dwi_bval_bvec,
        PIPELINE_INITIALISATION.out.rev_b0,
        PIPELINE_INITIALISATION.out.lesion,
        PIPELINE_INITIALISATION.out.covariates
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        SCILUS_SF_TRACTOMICS.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
