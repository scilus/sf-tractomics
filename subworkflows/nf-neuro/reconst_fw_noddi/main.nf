include { RECONST_DIFFUSIVITYPRIORS } from '../../../modules/nf-neuro/reconst/diffusivitypriors/main'
include { RECONST_MEANDIFFUSIVITYPRIORS } from '../../../modules/nf-neuro/reconst/meandiffusivitypriors/main'
include { RECONST_NODDI          } from '../../../modules/nf-neuro/reconst/noddi/main'
include { RECONST_FREEWATER      } from '../../../modules/nf-neuro/reconst/freewater/main'
include { RECONST_DTIMETRICS as FW_CORRECTED_DTIMETRICS } from '../../../modules/nf-neuro/reconst/dtimetrics/main'


workflow RECONST_FW_NODDI {

    take:
        dwi_bval_bvec   // channel: [ val(meta), dwi, bval, bvec ]
        brain_mask      // channel: [ val(meta), brain_mask ]
        fa_ad_rd_md     // channel: [ val(meta), fa, ad, rd, md ]
        diffusivities   // multiMap channel with para_diff, iso_diff, perp_diff_min, perp_diff_max
    main:
        ch_versions = channel.empty()

        // Make sure that at least one of the two reconstructions is requested
        if (!params.run_noddi && !params.run_freewater) {
            error "At least one of params.run_noddi or params.run_freewater must be true to run this subworkflow."
        }

        // Prepase base input channels
        ch_base = dwi_bval_bvec.join(brain_mask)

        // Classify priors based on their format
        ch_para = diffusivities.para_diff
            .branch { item ->
                per_subject: item instanceof List && item.size() == 2 && item[0] instanceof Map && item[0].containsKey('id')
                global: true
            }
        ch_iso = diffusivities.iso_diff
            .branch { item ->
                per_subject: item instanceof List && item.size() == 2 && item[0] instanceof Map && item[0].containsKey('id')
                global: true
            }
        ch_perp_min = diffusivities.perp_diff_min
            .branch { item ->
                per_subject: item instanceof List && item.size() == 2 && item[0] instanceof Map && item[0].containsKey('id')
                global: true
            }
        ch_perp_max = diffusivities.perp_diff_max
            .branch { item ->
                per_subject: item instanceof List && item.size() == 2 && item[0] instanceof Map && item[0].containsKey('id')
                global: true
            }

        // Combine all counts into a single validation tuple and check
        ch_para.per_subject.count()
            .combine(ch_iso.per_subject.count())
            .combine(ch_perp_min.per_subject.count())
            .combine(ch_perp_max.per_subject.count())
            .combine(ch_para.global.count())
            .combine(ch_iso.global.count())
            .combine(ch_perp_min.global.count())
            .combine(ch_perp_max.global.count())
            .map { para_s, iso_s, perp_min_s, perp_max_s, para_g, iso_g, perp_min_g, perp_max_g ->

                // NODDI validation: para and iso must be provided together
                if ( params.run_noddi ) {
                    def has_para = (para_s > 0 || para_g > 0)
                    def has_iso = (iso_s > 0 || iso_g > 0)
                    if ( has_para != has_iso ) {
                        error "For NODDI reconstruction, both para_diff and iso_diff must be provided together, " +
                            "either as per-subject values or as global values."
                    }
                }

                // Freewater validation: all four priors must be provided together
                if ( params.run_freewater ) {
                    def has_para = (para_s > 0 || para_g > 0)
                    def has_iso = (iso_s > 0 || iso_g > 0)
                    def has_perp_min = (perp_min_s > 0 || perp_min_g > 0)
                    def has_perp_max = (perp_max_s > 0 || perp_max_g > 0)
                    if ( has_para && (has_para != has_iso || has_para != has_perp_min || has_para != has_perp_max) ) {
                        error "For Freewater Elimination reconstruction, para_diff, iso_diff, perp_diff_min and " +
                            "perp_diff_max must be provided together, either as per-subject values or as global values."
                    }
                }

                // Warn if per-subject priors are provided alongside average_diff_priors
                if ( params.average_diff_priors && (para_s > 0 || iso_s > 0 || perp_min_s > 0 || perp_max_s > 0) ) {
                    log.warn "Options.average_diff_priors is set to true, but per-subject diffusivity priors " +
                            "were provided. The per-subject diffusivity priors will be ignored and the computed " +
                            "diffusivity priors will be averaged across subjects."
                }

                return true
            }
            .subscribe {}

        // Prepare NODDI inputs. This channel will be combined/joined in the
        // lines that follow with diffusivity priors w.r.t the following 3 scenarios:
        // Option 1: The user specifies the diffusivity priors to use (via params.para_diff and params.iso_diff).
        // Option 2: The user provides global diffusivity priors to be used across all subjects.
        // Option 3: The user wants to compute diffusivity priors for each subject individually or average them across subjects (recommended).

        // Branch 1: Custom diffusivity priors provided per-subject
        ch_noddi_custom_subj = ch_base
            .join( ch_para.per_subject )
            .join( ch_iso.per_subject )
        ch_freewater_custom_subj = ch_base
            .join( ch_para.per_subject )
            .join( ch_iso.per_subject )
            .join( ch_perp_min.per_subject )
            .join( ch_perp_max.per_subject )

        // Branch 2: Custom diffusivity priors provided (single value across subjects)
        ch_noddi_custom = ch_base
            .combine( ch_para.global )
            .combine( ch_iso.global )
        ch_freewater_custom = ch_base
            .combine( ch_para.global )
            .combine( ch_iso.global )
            .combine( ch_perp_min.global )
            .combine( ch_perp_max.global )

        // Branch 3: Compute diffusivity priors
        ch_subjects_with_custom_priors = ch_para.per_subject
            .map { meta, _value -> meta }
        ch_has_global_prior = ch_para.global.count()
        ch_compute_diff_priors = ch_has_global_prior
            .combine( fa_ad_rd_md )
            .filter { count, _meta, _fa, _ad, _rd, _md -> count == 0 } // Only compute diffusivity priors if no global prior is provided
            .map{ _count, meta, fa, ad, rd, md ->
                return [meta, fa, ad, rd, md]
            }

        // This should not be the case, but removing subjects that have custom priors
        ch_compute_diff_priors = ch_compute_diff_priors
            .map { meta, fa, ad, rd, md -> [meta.id, meta, fa, ad, rd, md] }
            .join(
                ch_subjects_with_custom_priors.map{ meta -> [meta.id, true] },
                remainder: true
            )
            .filter { items -> items[-1] == null } // Keep only those that do not have custom priors
            .map { _id, meta, fa, ad, rd, md, _null -> [meta, fa, ad, rd, md] }

        RECONST_DIFFUSIVITYPRIORS( ch_compute_diff_priors )
        ch_versions = ch_versions.mix(RECONST_DIFFUSIVITYPRIORS.out.versions)

        // Then compute mean diffusivity priors across subjects.
        if ( params.average_diff_priors ) {
            log.warn "Options.average_diff_priors is set to true. Averaging diffusivity priors across subjects. " +
                "This is not recommended, as it applies the same diffusivity priors to all subjects, which may not be optimal " +
                "if you have a wide age range."
            RECONST_MEANDIFFUSIVITYPRIORS(
                RECONST_DIFFUSIVITYPRIORS.out.para_diff_file
                    .map{ _meta, path -> path }
                    .collect(),
                RECONST_DIFFUSIVITYPRIORS.out.iso_diff_file
                    .map{ _meta, path -> path }
                    .collect(),
                RECONST_DIFFUSIVITYPRIORS.out.perp_diff_file
                    .map{ _meta, path -> path }
                    .collect()
            )
            ch_versions = ch_versions.mix(RECONST_MEANDIFFUSIVITYPRIORS.out.versions)

            ch_noddi_computed = ch_base
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.mean_para_diff)
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.mean_iso_diff)
            ch_freewater_computed = ch_base
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.mean_para_diff)
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.mean_iso_diff)
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.min_perp_diff)
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.max_perp_diff)
        }
        else {
            ch_noddi_computed = ch_base
                .join(RECONST_DIFFUSIVITYPRIORS.out.mean_para_diff)
                .join(RECONST_DIFFUSIVITYPRIORS.out.mean_iso_diff)
            ch_freewater_computed = ch_base
                .join(RECONST_DIFFUSIVITYPRIORS.out.mean_para_diff)
                .join(RECONST_DIFFUSIVITYPRIORS.out.mean_iso_diff)
                .join(RECONST_DIFFUSIVITYPRIORS.out.min_perp_diff)
                .join(RECONST_DIFFUSIVITYPRIORS.out.max_perp_diff)
        }

        if ( params.run_noddi ) {
            ch_noddi_input = ch_noddi_custom_subj
                .mix( ch_noddi_custom )
                .mix( ch_noddi_computed )
                .filter{ _meta, _dwi, bval, _bvec, _b0_mask, _para, _iso ->
                    def is_multi_shell = bval.text.tokenize().unique().size() > 2
                    if (!is_multi_shell && !params.noddi_silence_single_shell_warnings){
                        log.warn "Subject ${_meta.id} has single-shell data. Skipping NODDI reconstruction."
                    }
                    return is_multi_shell
                }
                .map{ meta, dwi, bval, bvec, b0_mask, para, iso ->
                    [meta, dwi, bval, bvec, b0_mask, [], para, iso] }

            RECONST_NODDI( ch_noddi_input )
            ch_versions = ch_versions.mix(RECONST_NODDI.out.versions)
        }

        if ( params.run_freewater ) {
            ch_freewater_input = ch_freewater_custom_subj
                .mix( ch_freewater_custom )
                .mix( ch_freewater_computed )
                .map{ meta, dwi, bval, bvec, b0_mask, para, iso, perp_min, perp_max ->
                    [meta, dwi, bval, bvec, b0_mask, [], para, iso, perp_min, perp_max] }

            RECONST_FREEWATER( ch_freewater_input )
            ch_versions = ch_versions.mix(RECONST_FREEWATER.out.versions)

            // -- Need to reprocess RECONST_DTIMETRICS to get
            //  FW corrected FA, MD, RD, AD, etc.
            //  using the FW corrected DWI.
            ch_fw_corrected_dti_metrics = RECONST_FREEWATER.out.dwi_fw_corrected
                .join(dwi_bval_bvec)
                .join(brain_mask)
                .map {
                    // Remove the original dwi from the join
                    meta, dwi_fw_corrected, _dwi_orig, bval, bvec, b0_mask ->
                        [meta, dwi_fw_corrected, bval, bvec, b0_mask]
                }

            FW_CORRECTED_DTIMETRICS( ch_fw_corrected_dti_metrics )
            ch_versions = ch_versions.mix(FW_CORRECTED_DTIMETRICS.out.versions)
        }

    emit:
        // NODDI
        noddi_dir           = params.run_noddi ? RECONST_NODDI.out.dir : channel.empty()
        noddi_isovf         = params.run_noddi ? RECONST_NODDI.out.isovf : channel.empty()
        noddi_icvf          = params.run_noddi ? RECONST_NODDI.out.icvf : channel.empty()
        noddi_ecvf          = params.run_noddi ? RECONST_NODDI.out.ecvf : channel.empty()
        noddi_odi           = params.run_noddi ? RECONST_NODDI.out.odi : channel.empty()

        // Freewater Elimination
        fw_dwi              = params.run_freewater ? RECONST_FREEWATER.out.dwi_fw_corrected : channel.empty()
        fw_dir              = params.run_freewater ? RECONST_FREEWATER.out.dir : channel.empty()
        fw_fibervolume      = params.run_freewater ? RECONST_FREEWATER.out.fibervolume : channel.empty()
        fw_fwf              = params.run_freewater ? RECONST_FREEWATER.out.fwf : channel.empty()
        fw_nrmse            = params.run_freewater ? RECONST_FREEWATER.out.nrmse : channel.empty()

        fw_dti_tensor       = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.tensor : channel.empty()
        fw_dti_md           = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.md : channel.empty()
        fw_dti_rd           = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.rd : channel.empty()
        fw_dti_ad           = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.ad : channel.empty()
        fw_dti_fa           = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.fa : channel.empty()
        fw_dti_rgb          = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.rgb : channel.empty()
        fw_dti_peaks        = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.evecs_v1 : channel.empty()
        fw_dti_evecs        = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.evecs : channel.empty()
        fw_dti_evals        = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.evals : channel.empty()
        fw_dti_residual     = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.residual : channel.empty()
        fw_dti_ga           = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.ga : channel.empty()
        fw_dti_mode         = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.mode : channel.empty()
        fw_dti_norm         = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.norm : channel.empty()

        versions = ch_versions
}
