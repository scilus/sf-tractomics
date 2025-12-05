/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { QC_MULTIQC as MULTIQC  } from '../modules/nf-neuro/qc/multiqc/main'
include { QC_MULTIQC as MULTIQC_GLOBAL  } from '../modules/nf-neuro/qc/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_sf-tractomics_pipeline'
include { TRACTOFLOW             } from '../subworkflows/nf-neuro/tractoflow'
include { TRACTOGRAM_MATH as ENSEMBLE_TRACKING } from '../modules/nf-neuro/tractogram/math/main'
include { QC_TRACTOGRAM as QC_ENSEMBLE } from '../modules/nf-neuro/qc/tractogram/main'
include { ATLAS_IIT              } from '../subworkflows/nf-neuro/atlas_iit/main'
include { RECONST_SHSIGNAL       } from '../modules/nf-neuro/reconst/shsignal'
include { RECONST_FW_NODDI       } from '../subworkflows/nf-neuro/reconst_fw_noddi/main'
include { BUNDLE_SEG             } from '../subworkflows/nf-neuro/bundle_seg/main' addParams(run_easyreg: false)
include { REGISTRATION_ANTS as REGISTER_ATLAS_B0 } from '../modules/nf-neuro/registration/ants/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_ATLAS_BUNDLES } from '../modules/nf-neuro/registration/antsapplytransforms/main.nf'
include { STATS_METRICSINROI     } from '../modules/nf-neuro/stats/metricsinroi/main'
include { ATLAS_ROIMETRICS       } from '../subworkflows/nf-neuro/atlas_roimetrics/main'
include { TRACTOMETRY            } from '../subworkflows/nf-neuro/tractometry/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SF_TRACTOMICS {
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
    ch_sub_multiqc_files = Channel.empty()
    ch_global_multiqc_files = Channel.empty()
    ch_topup_config = Channel.empty()
    ch_bet_template = Channel.empty()
    ch_bet_probability = Channel.empty()

    /* Load topup config if provided */
    if ( params.config_topup ) {
        if ( file(params.config_topup).exists()) {
            ch_topup_config = Channel.fromPath(params.config_topup, checkIfExists: true)
        }
        else {
            ch_topup_config = Channel.value( params.config_topup )
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
    ch_sub_multiqc_files = ch_sub_multiqc_files.mix(TRACTOFLOW.out.mqc)
    // ch_global_multiqc_files = ch_global_multiqc_files.mix(TRACTOFLOW.out.global_mqc)

    //
    // Ensemble tracking
    //
    ch_input_tracking_qc = Channel.empty()
    if (params.run_local_tracking && params.run_pft_tracking) {
        ch_tractogram_math_input = TRACTOFLOW.out.pft_tractogram
            .join(TRACTOFLOW.out.local_tractogram)
            .map {
                meta, pft_tractogram, local_tractogram ->
                    [meta, [pft_tractogram, local_tractogram], []]}
        ENSEMBLE_TRACKING(ch_tractogram_math_input)
        ch_input_tracking_qc = ENSEMBLE_TRACKING.out.trk
    }
    else if (params.run_local_tracking || params.run_pft_tracking) {
        ch_input_tracking_qc = TRACTOFLOW.out.pft_tractogram
            .mix(TRACTOFLOW.out.local_tractogram)
    }

    QC_ENSEMBLE(ch_input_tracking_qc
        .join(TRACTOFLOW.out.wm_mask)
        .join(TRACTOFLOW.out.gm_mask))
    ch_sub_multiqc_files = ch_sub_multiqc_files.mix(QC_ENSEMBLE.out.mqc)
    ch_global_multiqc_files = ch_global_multiqc_files.mix(
        QC_ENSEMBLE.out.dice.map { _meta, dice_file -> dice_file } )
    ch_global_multiqc_files = ch_global_multiqc_files.mix(
        QC_ENSEMBLE.out.sc.map { _meta, sc_file -> sc_file } )

    //
    // Run RECONST/SH_METRICS
    //
    if (params.sh_fitting) {
        RECONST_SHSIGNAL(
            TRACTOFLOW.out.dwi
                .map{ it + [[]] }
        )
    }

    //
    // Run BundleSeg
    //
    ch_bundle_seg = Channel.empty()
    if (params.run_bundle_seg) {
        ch_input_bundle_seg = TRACTOFLOW.out.pft_tractogram
            .mix(TRACTOFLOW.out.local_tractogram)
            .groupTuple()
        BUNDLE_SEG(
            TRACTOFLOW.out.dti_fa,
            TRACTOFLOW.out.pft_tractogram
                .mix(TRACTOFLOW.out.local_tractogram)
                .groupTuple(),
            Channel.empty()
        )

        ch_versions = ch_versions.mix(BUNDLE_SEG.out.versions)
        ch_sub_multiqc_files = ch_sub_multiqc_files.mix(BUNDLE_SEG.out.mqc)
        ch_sub_multiqc_files = ch_sub_multiqc_files.mix(BUNDLE_SEG.out.bundles_mqc)
        ch_bundle_seg = BUNDLE_SEG.out.bundles
    }

    // Prepare volume ROI metric extraction
    // Start by collecting DTI metrics
    ch_input_metrics = TRACTOFLOW.out.dti_fa
        .join(TRACTOFLOW.out.dti_md)
        .join(TRACTOFLOW.out.dti_rd)
        .join(TRACTOFLOW.out.dti_ad)
        .join(TRACTOFLOW.out.afd_total)
        .join(TRACTOFLOW.out.afd_sum)
        .join(TRACTOFLOW.out.afd_max)

    //
    // Run RECONST/NODDI & RECONST/FREEWATER
    //
    if (params.run_noddi || params.run_freewater) {
        RECONST_FW_NODDI(
            TRACTOFLOW.out.dwi,
            TRACTOFLOW.out.b0_mask,
            TRACTOFLOW.out.dti_fa
                .join(TRACTOFLOW.out.dti_ad)
                .join(TRACTOFLOW.out.dti_rd)
                .join(TRACTOFLOW.out.dti_md)
        )
        ch_versions = ch_versions.mix(RECONST_FW_NODDI.out.versions)

        // Add FW/NODDI metrics to the volume
        // ROI extraction.
        ch_input_metrics = ch_input_metrics
            .join(RECONST_FW_NODDI.out.fw_fw)
            .join(RECONST_FW_NODDI.out.fw_dti_fa)
            .join(RECONST_FW_NODDI.out.fw_dti_md)
            .join(RECONST_FW_NODDI.out.fw_dti_rd)
            .join(RECONST_FW_NODDI.out.fw_dti_ad)
            .join(RECONST_FW_NODDI.out.noddi_ndi)
            .join(RECONST_FW_NODDI.out.noddi_fwf)
            .join(RECONST_FW_NODDI.out.noddi_odi)
            .join(RECONST_FW_NODDI.out.noddi_ecvf)
    }

    ch_input_metrics = ch_input_metrics.map {tuple ->
        def meta = tuple[0]
        def metrics = tuple[1..-1]
        return [meta, metrics]
    }

    if (params.run_atlas_roimetrics) {
        ATLAS_ROIMETRICS(
            TRACTOFLOW.out.b0,
            ch_input_metrics,
            [ use_atlas_iit: params.use_atlas_iit ]
        )
        ch_versions = ch_versions.mix(ATLAS_ROIMETRICS.out.versions)

        // Collect all ROI stats into a single file
        // by appending each row of the TSV/CSV files,
        // while keeping the header from the first
        // file only and skipping it in the rest.
        ch_collection_mean_input = ATLAS_ROIMETRICS.out.stats_tab_mean
            .map{ _meta, stats_tab -> stats_tab }
            .collectFile(
                storeDir: "${params.outdir}/stats/",
                name: "aggregated_atlas-iit_label-mean_desc-roi_stats.tsv",
                skip: 1,
                keepHeader: true)
        ch_global_multiqc_files = ch_global_multiqc_files.mix(ch_collection_mean_input)

        ch_collection_std_input = ATLAS_ROIMETRICS.out.stats_tab_std
            .map{ _meta, stats_tab -> stats_tab }
            .collectFile(
                storeDir: "${params.outdir}/stats/",
                name: "aggregated_atlas-iit_label-std_desc-roi_stats.tsv",
                skip: 1,
                keepHeader: true)
    }

    if ( params.run_tractometry ) {
        TRACTOMETRY(
            ch_bundle_seg,
            Channel.empty(),
            ch_input_metrics,
            Channel.empty(),
            TRACTOFLOW.out.fodf)
        ch_versions = ch_versions.mix(TRACTOMETRY.out.versions)

        ch_tractometry_mqc = TRACTOMETRY.out.mean_tsv
            .map{ _meta, stats -> stats }
            .collectFile(
                storeDir: "${params.outdir}/metrics/",
                name: "bundles_mean_stats.tsv",
                skip: 1,
                keepHeader: true
            )
        ch_global_multiqc_files = ch_global_multiqc_files.mix(ch_tractometry_mqc)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'sf-tractomics_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = Channel.empty() // To store versions, methods description, etc.

    ch_multiqc_config_subject = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_config_global = Channel.fromPath(
        "$projectDir/assets/multiqc_config_global.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.fromPath("$projectDir/assets/sf-tractomics-multiqc-logo.png", checkIfExists: true)

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

    // Reorganizing subject-specific multiqc files here.
    qc_files = ch_sub_multiqc_files
        .groupTuple()
        .map{ meta, files ->
            def f = files.flatten().findAll { it != null }
            return tuple(meta, f)
        }

    MULTIQC (
        qc_files,
        ch_multiqc_files.collect(),
        ch_multiqc_config_subject.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    ch_fd_files = ch_sub_multiqc_files
        .filter { _meta, files ->
            files.any { it.name.contains("dwi_eddy_restricted_movement_rms") }
        }
        .map{ it[1] }
    ch_global_multiqc_files = ch_global_multiqc_files.mix(ch_fd_files.flatten())
    ch_global_multiqc_files = ch_global_multiqc_files.mix(ch_multiqc_files)

    // Global multiqc
    MULTIQC_GLOBAL (
        Channel.of([meta:[id: 'global'], qc_images: []]),
        ch_global_multiqc_files.collect(),
        ch_multiqc_config_global.toList(),
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
