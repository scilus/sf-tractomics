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
include { REGISTRATION_ANTS as REGISTER_ATLAS_BUNDLES } from '../modules/nf-neuro/registration/ants/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_ATLAS_BUNDLES } from '../modules/nf-neuro/registration/antsapplytransforms/main.nf'
include { BUNDLE_IIT             } from '../modules/local/bundle/iit/main'
include { VOLUME_ROISTATS        } from '../modules/local/volume/roistats/main'
include { VOLUME_COLLECTSTATS    } from '../modules/local/volume/collectstats/main'

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

    if (params.run_atlas_based_tractometry) {
        //
        // IIT ATLAS
        //

        // Extract bundle masks from IIT atlas
        ch_iit_template_bundles = channel.fromPath( params.iit_atlas.bundle_masks_dir + "/*.nii.gz", checkIfExists: true ).collect()
        ch_iit_template_thr = channel.fromPath( params.iit_atlas.bundle_masks_thresholds, checkIfExists: true )
        BUNDLE_IIT(ch_iit_template_bundles, ch_iit_template_thr)

        // Register IIT atlas to subject space
        ch_iit_template_b0 = channel.fromPath( params.iit_atlas.template_b0 )
        ch_input_register_iit = TRACTOFLOW.out.b0
            .combine(ch_iit_template_b0)
            .map{ meta, b0, template_b0 -> [meta, b0, template_b0, []] }
        REGISTER_ATLAS_BUNDLES(ch_input_register_iit)

        // Apply the transformation to subject space to the bundles
        ch_iit_transform_bundles = TRACTOFLOW.out.b0
            .join(REGISTER_ATLAS_BUNDLES.out.forward_image_transform)
            .combine(BUNDLE_IIT.out.bundle_masks.toList())
            .map {
                meta, b0, transform, bundles ->
                    [meta, bundles, b0, transform]
            }
        TRANSFORM_ATLAS_BUNDLES(ch_iit_transform_bundles)
        ch_versions = ch_versions.mix(TRANSFORM_ATLAS_BUNDLES.out.versions)

        // Prepare volume ROI metric extraction
        // Start by collecting DTI metrics
        ch_input_volume_roistats = TRACTOFLOW.out.dti_fa
            .join(TRACTOFLOW.out.dti_md)
            .join(TRACTOFLOW.out.dti_rd)
            .join(TRACTOFLOW.out.dti_ad)

        //
        // EXTRACT ROI VOLUME STATISTICS
        //
        // Input: [meta, [metrics_list], [masks]]

        ch_input_volume_roistats = ch_input_volume_roistats
            .map {tuple ->
                def meta = tuple[0]
                def metrics = tuple[1..-1]
                return [meta, metrics]
            }
            .join(TRANSFORM_ATLAS_BUNDLES.out.warped_image)
        VOLUME_ROISTATS(ch_input_volume_roistats)

        //
        // COLLECT/GROUP ROI STATS
        //
        ch_iit_roi_stats = VOLUME_ROISTATS.out.stats_csv.collect()
        VOLUME_COLLECTSTATS(ch_iit_roi_stats)
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
