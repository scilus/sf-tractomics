/**
 * Parse default options from a subworkflow meta.yml file
 * @param metaFilePath Path to the meta.yml file
 * @return Map containing default option values extracted from the options input definition
 */
def parseDefaultsFromMeta(String metaFilePath) {
    def metaFile = new File(metaFilePath)
    if (!metaFile.exists()) {
        error "Meta file not found: ${metaFilePath}"
        return [:]
    }

    def yaml = new groovy.yaml.YamlSlurper().parse(metaFile)
    def defaults = [:]

    // Extract defaults from the 'options' input definition in meta.yml
    if (yaml.input) {
        yaml.input.each { inputDef ->
            // Look for the 'options' input which contains the entries with default values
            if (inputDef.containsKey('options')) {
                def optionsInput = inputDef.options

                // The 'entries' field contains each option with its default value
                if (optionsInput.containsKey('entries')) {
                    def entries = optionsInput.entries

                    // Extract default value from each entry
                    entries.each { key, value ->
                        if (value.containsKey('default')) {
                            defaults[key] = value.default
                        }
                    }
                }
            }
        }
    }

    if (defaults.isEmpty()) {
        log.warn "Could not find default values in meta.yml file: ${metaFilePath}"
    }

    return defaults
}

/**
 * Merge provided options with defaults from meta.yml
 * @param provided Map of provided options
 * @param defaults Map of default options
 * @param strict If true, only allow keys that exist in defaults
 * @return Map containing merged options with defaults filled in
 */
def mergeWithDefaults(Object provided, Object defaults, boolean strict = false) {

    // Check to ensure options is a Map
    if (!(provided instanceof Map)) {
        throw new IllegalArgumentException("Provided options must be a Map, but received ${provided?.getClass()?.getName() ?: 'null'}")
    }
    if (!(defaults instanceof Map)) {
        throw new IllegalArgumentException("Default options must be a Map, but received ${defaults?.getClass()?.getName() ?: 'null'}")
    }

    // Create a new map and populate with defaults
    def merged = new HashMap()
    if (defaults) {
        merged.putAll(defaults)
    }

    if (strict) {
        // Only allow options that exist in defaults
        provided.each { key, value ->
            if (defaults.containsKey(key)) {
                merged[key] = value
            } else {
                log.warn "Unknown option '${key}' will be ignored (not in defaults)"
            }
        }
    } else {
        // Allow all provided options
        if (provided) {
            merged.putAll(provided)
        }
    }

    return merged
}

/**
 * Convenience function to merge options with defaults from meta.yml
 * @param options Map of provided options
 * @param metaPath Path to the meta.yml file (can use ${moduleDir}/meta.yml)
 * @return Map containing merged options with defaults filled in
 */
def getOptionsWithDefaults(Object options, String metaPath, boolean strict = true) {

    // Check to ensure options is a Map
    if (!(options instanceof Map)) {
        throw new IllegalArgumentException("Options must be a Map, but received ${options?.getClass()?.getName() ?: 'null'}")
    }

    def defaults = parseDefaultsFromMeta(metaPath)
    return mergeWithDefaults(options, defaults, strict)
}

workflow UTILS_OPTIONS {
    take:
        meta_file    // file: path(.../meta.yml)
        options      // map: val(options)
        strict       // bool: val(strict) - if true, only allow options that exist in defaults

    main:
        ch_versions = channel.empty()
        // Capture options in a local variable to avoid dataflow broadcast issues
        def provided_options = options
        def strict_mode = strict

        merged_options = getOptionsWithDefaults(provided_options, meta_file, strict_mode)
    emit:
        options = merged_options         // DataflowValue (map): [ options ]
        versions = ch_versions           // channel: [ versions.yml ]
}
