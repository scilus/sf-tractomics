
process IO_SAFECASTINPUTS {
    label 'process_single'

    input:
        tuple val(meta), path(dwi, stageAs: "dwi"), path(bval, stageAs: "bval"), path(bvec, stageAs: "bvec"), path(sbref, stageAs: "sbref"), path(rev_dwi, stageAs: "rev"), path(rev_bval, stageAs: "rval"), path(rev_bvec, stageAs: "rvec"), path(rev_sbref, stageAs: "rbref"), path(t1, stageAs: "t1"), path(wmparc, stageAs: "wmparc"), path(aparc_aseg, stageAs: "aparc+aseg"), path(lesion, stageAs: "lesion")
    output:
        tuple val(meta), path("$out_dwi"), path("$out_bval"), path("$out_bvec"), path("$out_sbref"), path("$out_rev_dwi"), path("$out_rev_bval"), path("$out_rev_bvec"), path("$out_rev_sbref"), path("$out_t1"), path("$out_wmparc"), path("$out_aparc_aseg"), path("$out_lesion"), emit: safe_inputs
    script:
        out_dwi = dwi ? "dwi.nii.gz" : "$dwi"
        out_bval = bval ? "dwi.bval" : "$bval"
        out_bvec = bvec ? "dwi.bvec" : "$bvec"
        out_sbref = sbref ? "sbref.nii.gz" : "$sbref"
        out_rev_dwi = rev_dwi ? "rev_dwi.nii.gz" : "$rev_dwi"
        out_rev_bval = rev_bval ? "rev_dwi.bval" : "$rev_bval"
        out_rev_bvec = rev_bvec ? "rev_dwi.bvec" : "$rev_bvec"
        out_rev_sbref = rev_sbref ? "rev_sbref.nii.gz" : "$rev_sbref"
        out_t1 = t1 ? "t1.nii.gz" : "$t1"
        out_wmparc = wmparc ? "wmparc.nii.gz" : "$wmparc"
        out_aparc_aseg = aparc_aseg ? "aparc+aseg.nii.gz" : "$aparc_aseg"
        out_lesion = lesion ? "lesion.nii.gz" : "$lesion"
    """
    [ -f "$dwi" ] && ln -sf $dwi dwi.nii.gz
    [ -f "$bval" ] && ln -sf $bval dwi.bval
    [ -f "$bvec" ] && ln -sf $bvec dwi.bvec
    [ -f "$sbref" ] && ln -sf $sbref sbref.nii.gz
    [ -f "$rev_dwi" ] && ln -sf $rev_dwi rev_dwi.nii.gz
    [ -f "$rev_bval" ] && ln -sf $rev_bval rev_dwi.bval
    [ -f "$rev_bvec" ] && ln -sf $rev_bvec rev_dwi.bvec
    [ -f "$rev_sbref" ] && ln -sf $rev_sbref rev_sbref.nii.gz
    [ -f "$t1" ] && ln -sf $t1 t1.nii.gz
    [ -f "$wmparc" ] && ln -sf $wmparc wmparc.nii.gz
    [ -f "$aparc_aseg" ] && ln -sf $aparc_aseg aparc+aseg.nii.gz
    [ -f "$lesion" ] && ln -sf $lesion lesion.nii.gz
    exit 0
    """
}
