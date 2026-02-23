include { IMAGE_MATH as THR_BUNDLE_MASK } from '../../../modules/nf-neuro/image/math/main'
include { IMAGE_MATH as SMOOTH_MASK } from '../../../modules/nf-neuro/image/math/main'
include { IMAGE_MATH as THR_SMOOTHED_MASK } from '../../../modules/nf-neuro/image/math/main'
include { UTILS_OPTIONS } from '../utils_options/main'

def download_file(url, output_path) {
    HttpURLConnection connection = new URL(url).openConnection()
    connection.setInstanceFollowRedirects(true)
    connection.setRequestProperty("User-Agent", "nf-neuro/atlas_iit (subworkflow)")
    connection.setRequestProperty("Accept", "*/*")
    connection.setRequestProperty("Accept-Encoding", "identity")
    connection.setRequestProperty("Connection", "Keep-Alive")
    connection.connect()

    if (connection.responseCode != 200) {
        error "Failed to download file: HTTP ${connection.responseCode}"
    }

    new File(output_path).withOutputStream { out ->
        out << connection.inputStream
    }
}

// Fetch IIT Atlas Mean B0
def fetch_iit_atlas_b0(b0Url, dest) {
    def outFile = new File("$dest/IITmean_b0.nii.gz")

    download_file(b0Url, outFile.absolutePath)

    return outFile
}

// Fetch Bundles Track Density Maps
// Which are later used to generate bundle masks based on thresholds.
def fetch_iit_atlas_tdi(bundleMapsUrl, dest, thresholds) {
    def output_dir = "${workflow.workDir}/atlas_iit/bundles/tdi/"

    // If files are all there, immediately return the directory
    if (allBundleFilesExist(thresholds, new File(output_dir + "bundle_maps/IIT_bundles"))) {
        return output_dir + "bundle_maps/IIT_bundles"
    }

    def intermediate_dir = new File("$dest/intermediate")
    intermediate_dir.mkdirs()

    download_file(bundleMapsUrl, "$intermediate_dir/IIT_bundles.zip")

    def bundleMapsFile = new java.util.zip.ZipFile("$intermediate_dir/IIT_bundles.zip")
    bundleMapsFile.entries().each{ it ->
        def path = java.nio.file.Paths.get("$dest/bundle_maps/" + it.name)
        if (it.directory) {
            java.nio.file.Files.createDirectories(path)
        }
        else {
            def parentDir = path.getParent()
            if (!java.nio.file.Files.exists(parentDir)) {
                java.nio.file.Files.createDirectories(parentDir)
            }
            java.nio.file.Files.copy(
                bundleMapsFile.getInputStream(it),
                path,
                java.nio.file.StandardCopyOption.REPLACE_EXISTING)
        }
    }

    return output_dir + "bundle_maps/IIT_bundles"
}

def get_tdi_thresholds() {
    // Those thresholds were recommended by the authors of the IIT atlas
    // to generate the bundle masks from the track density images.
    // See this guide: https://www.nitrc.org/frs/download.php/12417/IIT_Atlas_v.5.0_MANUAL.pdf
    return [
        "AC": 3,
        "AF_L": 2,
        "AF_R": 2,
        "AST_L": 5,
        "AST_R": 5,
        "C_L": 3,
        "C_R": 3,
        "CC_ForcepsMajor": 5,
        "CC_ForcepsMinor": 2,
        "CC": 1,
        "CCMid": 4,
        "CST_L": 0,
        "CST_R": 0,
        "F_L_R": 0,
        "FPT_L": 1,
        "FPT_R": 1,
        "ICP_L": 5,
        "ICP_R": 5,
        "IFOF_L": 0,
        "IFOF_R": 0,
        "ILF_L": 5,
        "ILF_R": 2,
        "MCP": 5,
        "MdLF_L": 5,
        "MdLF_R": 5,
        "ML_L": 150,
        "ML_R": 150,
        "OPT_L": 0,
        "OPT_R": 0,
        "OR_L": 5,
        "OR_R": 10,
        "PPT_L": 2,
        "PPT_R": 2,
        "SCP": 50,
        "SLF_L": 10,
        "SLF_R": 5,
        "STT_L": 150,
        "STT_R": 150,
        "UF_L": 0,
        "UF_R": 0,
        "VOF_L": 15,
        "VOF_R": 3
    ]
}

// Make sure there's a bundle file for each threshold.
// This is to avoid re-downloading if files ALL the
// files are already there.
boolean allBundleFilesExist(Map thresholds, File dir) {
    thresholds.every { key, _value ->
        new File(dir, "${key}.nii.gz").exists()
    }
}

workflow ATLAS_IIT {
    take:
        options             // Map of options [ options ]

    main:
        ch_versions = channel.empty()

        // Merge options with defaults from meta.yml
        UTILS_OPTIONS("${moduleDir}/meta.yml", options, true)
        options = UTILS_OPTIONS.out.options.value

        def input_b0 = options.atlas_iit_b0 ?: null
        def input_bundle_masks_dir = options.atlas_iit_bundle_masks_dir ?: null

        // Fetch Mean B0
        if ( input_b0 ) {
            ch_b0 = channel.fromPath(input_b0, checkIfExists: true)
        }
        else {
            new File("${workflow.workDir}/atlas_iit/").mkdirs()

            if (!new File("${workflow.workDir}/atlas_iit/IITmean_b0.nii.gz").exists()) {
                fetch_iit_atlas_b0(
                    "https://www.nitrc.org/frs/download.php/11266/IITmean_b0.nii.gz",
                    "${workflow.workDir}/atlas_iit/"
                )
            }

            ch_b0 = channel.fromPath("${workflow.workDir}/atlas_iit/IITmean_b0.nii.gz", checkIfExists: true)
        }

        // Fetch and Process Bundle Masks
        if (input_bundle_masks_dir) {
            ch_bundles = channel.fromPath(input_bundle_masks_dir + "/*.nii.gz", checkIfExists: true)
                .collect(sort: { path_a, path_b ->
                    def name_a = path_a.getName()
                    def name_b = path_b.getName()
                    name_a <=> name_b
                })
        }
        else {
            def thresholds = get_tdi_thresholds()
            def atlas_tdi = fetch_iit_atlas_tdi(
                "https://www.nitrc.org/frs/download.php/11472/IIT_bundles.zip",
                "${workflow.workDir}/atlas_iit/bundles/tdi",
                thresholds
            )

            ch_bundles = channel.fromPath(atlas_tdi + "/*.nii.gz", checkIfExists: true)

            // Enabling this option will convert the track density maps to binary
            // bundle masks based on the recommended thresholds.
            //
            // One might choose to disable this thresholding step if they want to
            // use the track density maps as "soft" bundle masks instead of binary
            // masks.
            //
            if (options.threshold_bundles) {
                // Pair all bundle maps with their respective thresholds
                ch_bundle_maps_with_thresholds = ch_bundles.map { file ->
                    def file_base_name = file.baseName.replace(".nii.gz", "").replace(".nii", "")
                    def thr_find = thresholds.find { line -> line.key == file_base_name }?.value
                    def thr = thr_find != null ? thr_find : null
                    def meta = [ id: file_base_name ]
                    return [meta, file, thr]
                }

                // Threshold the track density maps to get binary bundle masks
                THR_BUNDLE_MASK(ch_bundle_maps_with_thresholds)
                ch_versions = ch_versions.mix(THR_BUNDLE_MASK.out.versions).first()

                // Smooth the mask with a gaussian kernel of sigma 1
                ch_smooth_mask_input = THR_BUNDLE_MASK.out.image
                    .map { meta, image -> [meta, image, options.smooth_sigma ?: 1.0] }
                SMOOTH_MASK(ch_smooth_mask_input)
                ch_versions = ch_versions.mix(SMOOTH_MASK.out.versions).first()

                // Threshold the smoothed mask with a threshold of 0.5 to get a binary mask again
                ch_thr_smoothed_mask_input = SMOOTH_MASK.out.image
                    .map { meta, image -> [meta, image, 0.5] }
                THR_SMOOTHED_MASK(ch_thr_smoothed_mask_input)
                ch_versions = ch_versions.mix(THR_SMOOTHED_MASK.out.versions).first()

                ch_bundles = THR_SMOOTHED_MASK.out.image
                    .map { _meta, mask -> mask }
            }

            ch_bundles = ch_bundles
                .collect(sort: { path_a, path_b ->
                    def name_a = path_a.getName()
                    def name_b = path_b.getName()
                    name_a <=> name_b
                })
        }

    emit:
        b0 = ch_b0
        bundles = ch_bundles
        versions = ch_versions
}
