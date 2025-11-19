def file_exists_and_not_empty(file_path) {
    if (file_path == null) {
        return false
    }
    if (file_path instanceof List && file_path.size() == 0) {
        return false
    }
    return file_path.exists()
}

process ATLAS_IIT {
    tag 'atlas'
    label 'process_single'

    container "scilus/scilus:2.2.1"

    input:
    path iit_b0
    path iit_bundles

    output:
    path "iit_bundles/*.nii.gz"      , emit: bundles, includeInputs: true
    path "*_b0.nii.gz"               , emit: b0, includeInputs: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def intermediate_dir = "maps"
    def output_dir = "iit_bundles"
    def thresholds_txt_file = "${intermediate_dir}/thresholds.txt"
    """
    #!/bin/bash

    mkdir -p "/tmp"

    export HOME="/tmp"

    # Download the B0
    ${file_exists_and_not_empty(iit_b0) ? "cp ${iit_b0} IITmean_b0.nii.gz" : "wget https://www.nitrc.org/frs/download.php/11266/IITmean_b0.nii.gz -O IITmean_b0.nii.gz" }

    # If the iit_bundles is not null and not empty, copy it locally
    echo "» Preparing IIT bundle masks... (input is ${iit_bundles})"

    if [ "${iit_bundles}" != "null" ] && [ -s "${iit_bundles}" ]; then
        echo "» Using provided IIT bundle masks from ${iit_bundles} -> ${output_dir}/"
        mv "${iit_bundles}" "${output_dir}"
    else
        mkdir -p "${intermediate_dir}"
        mkdir -p "${output_dir}"

        # Download the bundle density maps.
        wget https://www.nitrc.org/frs/download.php/11472/IIT_bundles.zip -O "IIT_bundles.zip"

        unzip -n "IIT_bundles.zip"
        mv "IIT_bundles"/* "${intermediate_dir}/"
        rm "IIT_bundles.zip" "__MACOSX" "IIT_bundles" -r

        # These thresholds are recommended by the IIT authors to create binary masks
        # from the bundle density maps.
        cat <<-END_THR > ${thresholds_txt_file}
            AC 3
            AF_L 2
            AF_R 2
            AST_L 5
            AST_R 5
            C_L 3
            C_R 3
            CC_ForcepsMajor 5
            CC_ForcepsMinor 2
            CC 1
            CCMid 4
            CST_L 0
            CST_R 0
            F_L_R 0
            FPT_L 1
            FPT_R 1
            ICP_L 5
            ICP_R 5
            IFOF_L 0
            IFOF_R 0
            ILF_L 5
            ILF_R 2
            MCP 5
            MdLF_L 5
            MdLF_R 5
            ML_L 150
            ML_R 150
            OPT_L 0
            OPT_R 0
            OR_L 5
            OR_R 10
            PPT_L 2
            PPT_R 2
            SCP 50
            SLF_L 10
            SLF_R 5
            STT_L 150
            STT_R 150
            UF_L 0
            UF_R 0
            VOF_L 15
            VOF_R 3
END_THR

        # Create binary masks using the thresholds
        # and adjust the stride for MNI compatibility.
        while read -r name thr; do
            echo "» Processing \${name} with threshold \${thr}"
            scil_volume_math lower_threshold "${intermediate_dir}/\${name}.nii.gz" "\${thr}" "${output_dir}/\${name}_mask.nii.gz" -f
            mrconvert -stride 1,2,3 "${output_dir}/\${name}_mask.nii.gz" "${output_dir}/\${name}_mask.nii.gz" -force
        done < ${thresholds_txt_file}
    fi

    echo "» All done! (masks created in ${output_dir})"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
        mrtrix: \$(mrcalc -version 2>&1 | sed -n 's/== mrcalc \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS

    """

    stub:
    def intermediate_dir = "maps"
    def output_dir = "iit_bundles"
    """
    mkdir -p ${intermediate_dir}
    mkdir -p ${output_dir}

    touch "IITmean_b0.nii.gz"

    for f in AC AF_L AF_R AST_L AST_R C_L C_R \
        CC_ForcepsMajor CC_ForcepsMinor CC CCMid \
        CST_L CST_R F_L_R FPT_L FPT_R ICP_L ICP_R \
        IFOF_L IFOF_R ILF_L ILF_R MCP MdLF_L MdLF_R \
        ML_L ML_R OPT_L OPT_R OR_L OR_R PPT_L PPT_R \
        SCP SLF_L SLF_R STT_L STT_R UF_L UF_R VOF_L \
        VOF_R; do
        touch "${output_dir}/\${f}_mask.nii.gz"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
        mrtrix: \$(mrcalc -version 2>&1 | sed -n 's/== mrcalc \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
