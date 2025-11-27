#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scilus/nf-tractoflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/scilus/nf-tractoflow
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NF_TRACTOFLOW  } from './workflows/nf-tractoflow'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_nf-tractoflow_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_nf-tractoflow_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SCILUS_NF_TRACTOFLOW {

    take:
    t1                  // channel: t1 read in from --input
    wmparc              // channel: wmparc read in from --input
    aparc_aseg          // channel: aparc_aseg read in from --input
    dwi_bval_bvec       // channel: dwi_bval_bvec read in from --input
    b0                  // channel: b0 read in from --input
    rev_dwi_bval_bvec   // channel: rev_dwi_bval_bvec read in from --input
    rev_b0              // channel: rev_b0 read in from --input
    lesion              // channel: lesion read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    NF_TRACTOFLOW (
        t1,
        wmparc,
        aparc_aseg,
        dwi_bval_bvec,
        b0,
        rev_dwi_bval_bvec,
        rev_b0,
        lesion
    )
    //emit:
    //multiqc_report = NF_TRACTOFLOW.out.multiqc_report // channel: /path/to/multiqc_report.html
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
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    SCILUS_NF_TRACTOFLOW (
        PIPELINE_INITIALISATION.out.t1,
        PIPELINE_INITIALISATION.out.wmparc,
        PIPELINE_INITIALISATION.out.aparc_aseg,
        PIPELINE_INITIALISATION.out.dwi_bval_bvec,
        PIPELINE_INITIALISATION.out.b0,
        PIPELINE_INITIALISATION.out.rev_dwi_bval_bvec,
        PIPELINE_INITIALISATION.out.rev_b0,
        PIPELINE_INITIALISATION.out.lesion
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs
//        SCILUS_NF_TRACTOFLOW.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
