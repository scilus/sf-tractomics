include { RECONST_DIFFUSIVITYPRIORS } from '../../../modules/nf-neuro/reconst/diffusivitypriors/main'
include { RECONST_MEANDIFFUSIVITYPRIORS } from '../../../modules/nf-neuro/reconst/meandiffusivitypriors/main'
include { RECONST_NODDI          } from '../../../modules/nf-neuro/reconst/noddi/main'
include { RECONST_FREEWATER      } from '../../../modules/nf-neuro/reconst/freewater/main'
include { RECONST_DTIMETRICS as FW_CORRECTED_DTIMETRICS } from '../../../modules/nf-neuro/reconst/dtimetrics/main'
workflow RECONST_FW_NODDI {

    take:
    dwi_bval_bvec
    brain_mask
    fa_ad_rd_md

    main:

    ch_versions = Channel.empty()

    // Make sure that at least one of the two reconstructions is requested
    if (!params.run_noddi && !params.run_freewater) {
        error "At least one of params.run_noddi or params.run_freewater must be true to run this subworkflow."
    }

    // Prepare NODDI inputs. This channel will be combined/joined in the
    // lines that follow with diffusivity priors w.r.t the following 3 scenarios:
    // Option 1: The user specifies the diffusivity priors to use (via params.para_diff and params.iso_diff).
    // Option 2: The user wants to compute the mean diffusivity priors across subjects. (Recommended)
    // Option 3: The user wants to compute diffusivity priors for each subject individually.

    if (params.run_noddi && ([params.iso_diff, params.para_diff].any()
        && ! [params.iso_diff, params.para_diff].every())) {
        error "Please provide both params.iso_diff and params.para_diff parameters to use custom diffusivity priors for NODDI."
    }
    else if (params.run_freewater
        && [params.iso_diff, params.para_diff, params.perp_diff_min, params.perp_diff_max].any()
        && ! [params.iso_diff, params.para_diff, params.perp_diff_min, params.perp_diff_max].every()) {
        error "Please provide all params.iso_diff, params.para_diff, params.perp_diff_min and params.perp_diff_max parameters to use custom "
            "diffusivity priors for Freewater Elimination. Otherwise, specify none and the priors will be "
            "automatically computed."
    }

    ch_noddi_input = dwi_bval_bvec
        .join(brain_mask)
    ch_freewater_input = dwi_bval_bvec
        .join(brain_mask)

    noddi_custom_priors = [params.para_diff, params.iso_diff].every()
    fw_custom_priors = [params.para_diff, params.iso_diff, params.perp_diff_min, params.perp_diff_max].every()
    if (((params.run_noddi && !params.run_freewater) && noddi_custom_priors) || (params.run_freewater && fw_custom_priors)) {
        // Use user-specified diffusivity priors across subjects.
        if (params.average_diff_priors) {
            log.warn "Both custom diffusivity priors and params.average_diff_priors parameter were provided."
                "The specified custom diffusivity priors will be used across subjects."
        }

        ch_noddi_input = ch_noddi_input
            .combine(Channel.value(params.para_diff))
            .combine(Channel.value(params.iso_diff))
        ch_freewater_input = ch_freewater_input
            .combine(Channel.value(params.para_diff))
            .combine(Channel.value(params.iso_diff))
            .combine(Channel.value(params.perp_diff_min))
            .combine(Channel.value(params.perp_diff_max))
    }
    else {
        // Compute diffusivity priors for each subject.
        RECONST_DIFFUSIVITYPRIORS(fa_ad_rd_md)
        ch_versions = ch_versions.mix(RECONST_DIFFUSIVITYPRIORS.out.versions)

        // Then compute mean diffusivity priors across subjects.
        if (params.average_diff_priors) {
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

            ch_noddi_input = ch_noddi_input
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.mean_para_diff)
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.mean_iso_diff)
            ch_freewater_input = ch_freewater_input
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.mean_para_diff)
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.mean_iso_diff)
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.min_perp_diff)
                .combine(RECONST_MEANDIFFUSIVITYPRIORS.out.max_perp_diff)
        }
        else {
            ch_noddi_input = ch_noddi_input
                .join(RECONST_DIFFUSIVITYPRIORS.out.mean_para_diff)
                .join(RECONST_DIFFUSIVITYPRIORS.out.mean_iso_diff)
            ch_freewater_input = ch_freewater_input
                .join(RECONST_DIFFUSIVITYPRIORS.out.mean_para_diff)
                .join(RECONST_DIFFUSIVITYPRIORS.out.mean_iso_diff)
                .join(RECONST_DIFFUSIVITYPRIORS.out.min_perp_diff)
                .join(RECONST_DIFFUSIVITYPRIORS.out.max_perp_diff)
        }
    }

    if (params.run_noddi) {
        ch_noddi_input = ch_noddi_input
            .map{ meta, dwi, bval, bvec, b0_mask, para, iso ->
                [meta, dwi, bval, bvec, b0_mask, [], para, iso] }

        RECONST_NODDI( ch_noddi_input )
        ch_versions = ch_versions.mix(RECONST_NODDI.out.versions)
    }

    if (params.run_freewater) {
        ch_freewater_input = ch_freewater_input
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
    noddi_dir           = params.run_noddi ? RECONST_NODDI.out.dir : Channel.empty()
    noddi_fwf           = params.run_noddi ? RECONST_NODDI.out.fwf : Channel.empty()
    noddi_ndi           = params.run_noddi ? RECONST_NODDI.out.ndi : Channel.empty()
    noddi_ecvf          = params.run_noddi ? RECONST_NODDI.out.ecvf : Channel.empty()
    noddi_odi           = params.run_noddi ? RECONST_NODDI.out.odi : Channel.empty()

    // Freewater Elimination
    fw_dwi              = params.run_freewater ? RECONST_FREEWATER.out.dwi_fw_corrected : Channel.empty()
    fw_dir              = params.run_freewater ? RECONST_FREEWATER.out.dir : Channel.empty()
    fw_fibervolume      = params.run_freewater ? RECONST_FREEWATER.out.fibervolume : Channel.empty()
    fw_fw               = params.run_freewater ? RECONST_FREEWATER.out.fw : Channel.empty()
    fw_nrmse            = params.run_freewater ? RECONST_FREEWATER.out.nrmse : Channel.empty()

    fw_dti_tensor       = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.tensor : Channel.empty()
    fw_dti_md           = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.md : Channel.empty()
    fw_dti_rd           = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.rd : Channel.empty()
    fw_dti_ad           = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.ad : Channel.empty()
    fw_dti_fa           = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.fa : Channel.empty()
    fw_dti_rgb          = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.rgb : Channel.empty()
    fw_dti_peaks        = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.evecs_v1 : Channel.empty()
    fw_dti_evecs        = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.evecs : Channel.empty()
    fw_dti_evals        = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.evals : Channel.empty()
    fw_dti_residual     = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.residual : Channel.empty()
    fw_dti_ga           = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.ga : Channel.empty()
    fw_dti_mode         = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.mode : Channel.empty()
    fw_dti_norm         = params.run_freewater ? FW_CORRECTED_DTIMETRICS.out.norm : Channel.empty()

    versions = ch_versions
}

