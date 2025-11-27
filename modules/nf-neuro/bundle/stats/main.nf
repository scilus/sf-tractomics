process BUNDLE_STATS {
    tag "$meta.id"
    label 'process_single'

    container "scilus/scilpy:2.2.1_cpu"

    input:
    tuple val(meta), path(bundles), path(labels_map), path(metrics), path(lesions)

    output:
    tuple val(meta), path("*__length_stats.json")               , emit: length, optional: true
    tuple val(meta), path("*__endpoints_map_raw.json")          , emit: endpoints_raw, optional: true
    tuple val(meta), path("*__endpoints_metric_stats.json")     , emit: endpoints_metric_stats, optional: true
    tuple val(meta), path("*__mean_std.json")                   , emit: mean_std, optional: true
    tuple val(meta), path("*__volume.json")                     , emit: volume, optional: true
    tuple val(meta), path("*__volume_lesions.json")             , emit: volume_lesions, optional: true
    tuple val(meta), path("*__streamline_count.json")           , emit: streamline_count, optional: true
    tuple val(meta), path("*__streamline_count_lesions.json")   , emit: streamline_count_lesions, optional: true
    tuple val(meta), path("*__volume_per_label.json")           , emit: volume_per_labels, optional: true
    tuple val(meta), path("*__volume_per_label_lesions.json")   , emit: volume_per_labels_lesions, optional: true
    tuple val(meta), path("*__mean_std_per_point.json")         , emit: mean_std_per_point, optional: true
    tuple val(meta), path("*__lesion_stats.json")               , emit: lesion_stats, optional: true
    tuple val(meta), path("*_endpoints_map_head.nii.gz")        , emit: endpoints_head, optional: true
    tuple val(meta), path("*_endpoints_map_tail.nii.gz")        , emit: endpoints_tail, optional: true
    tuple val(meta), path("*_lesion_map.nii.gz")                , emit: lesion_map, optional: true
    tuple val(meta), path("*_desc-mean_param-lesions_stats.tsv") , emit: mean_lesions_tsv, optional: true
    tuple val(meta), path("*_desc-point_param-lesions_stats.tsv"), emit: mean_per_point_lesions_tsv, optional: true
    tuple val(meta), path("*_desc-mean_stats.tsv")              , emit: mean_tsv, optional: true
    tuple val(meta), path("*_desc-point_stats.tsv")             , emit: mean_per_point_tsv, optional: true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def density_weighting = task.ext.density_weighting ? "--density_weighting" : ""
    def normalize_weights = task.ext.normalize_weights ? "--normalize_weights" : "--bin"
    def length_stats = task.ext.length_stats ?: ""
    def endpoints = task.ext.endpoints ?: ""
    def mean_std = task.ext.mean_std ?: ""
    def volume = task.ext.volume ?: ""
    def lesions_stats = task.ext.lesions_stats ?: ""
    def min_lesion_vol = task.ext.min_lesion_vol ?: ""
    def streamline_count = task.ext.streamline_count ?: ""
    def volume_per_labels = task.ext.volume_per_labels ?: ""
    def mean_std_per_point = task.ext.mean_std_per_point ?: ""

    """
    bundles=( ${bundles.join(" ")} )
    label_map=( ${labels_map.join(" ")} )
    metrics=( ${metrics.join(" ")} )

    for index in \${!bundles[@]};
    do\
        bname=\$(basename \${bundles[index]} .trk);
        bname=\${bname/${prefix}__/}
        bname=\${bname%%_labels_*}

        echo "Bundle name: \${bname}"

        # Initialize array for all relevant metrics
        b_metrics=()

        for m in \${!metrics[@]}; do
            # Include if: matches bname OR is not an afd_metric file
            if [[ "\${metrics[\$m]}" == *"\${bname}"* ]] || [[ "\${metrics[\$m]}" != *"afd_metric"* ]]; then
                b_metrics+=("\${metrics[\$m]}")
            fi
        done

        if [[ "$length_stats" ]];
        then
            scil_tractogram_print_info \${bundles[index]} > ${prefix}__\${bname}_length.json
        fi

        if [[ "$endpoints" ]];
        then
            scil_bundle_compute_endpoints_map \${bundles[index]} \
                ${prefix}__\${bname}_endpoints_map_head.nii.gz \
                ${prefix}__\${bname}_endpoints_map_tail.nii.gz --out_json \
                ${prefix}__\${bname}_endpoints_raw.json;

            scil_volume_stats_in_ROI ${prefix}__\${bname}_endpoints_map_head.nii.gz $normalize_weights\
                --metrics \${b_metrics[@]} > ${prefix}__\${bname}_head.json
            scil_volume_stats_in_ROI ${prefix}__\${bname}_endpoints_map_tail.nii.gz $normalize_weights\
                --metrics \${b_metrics[@]} > ${prefix}__\${bname}_tail.json;
            fi

        if [[ "$mean_std" ]];
        then
            scil_bundle_mean_std $density_weighting \${bundles[index]} \${b_metrics[@]} >\
                ${prefix}__\${bname}__std.json
        fi

        if [[ "$volume" ]];
        then
            scil_bundle_shape_measures \${bundles[index]} > ${prefix}__\${bname}_volume_stat.json

            if [[ "$lesions_stats" ]];
            then
                scil_lesions_info $lesions \${bname}_volume_lesions_stat.json \
                    --bundle \${bundles[index]} --out_lesion_stats ${prefix}__lesion_stats.json \
                    --out_streamlines_stats \${bname}__streamline_count_lesions_stat.json \
                    --min_lesion_vol $min_lesion_vol -f
            fi
        fi

        if [[ "$streamline_count" ]];
        then
            scil_tractogram_count_streamlines \${bundles[index]} > ${prefix}__\${bname}_streamlines.json
        fi

        if [[ "$volume_per_labels" ]];
        then
            scil_bundle_volume_per_label \${label_map[index]} \$bname --sort_keys >\
                ${prefix}__\${bname}_volume_label.json

            if [[ "$lesions_stats" ]];
            then
                scil_lesions_info $lesions ${prefix}__\${bname}_volume_per_label_lesions_stat.json \
                    --bundle_labels_map \${label_map[index]} \
                    --out_lesion_atlas "${prefix}__\${bname}_lesion_map.nii.gz" \
                    --min_lesion_vol $min_lesion_vol
            fi
        fi

        if [[ "$mean_std_per_point" ]];
        then
            scil_bundle_mean_std \${bundles[index]} \${b_metrics[@]}\
                --per_point \${label_map[index]} --sort_keys $density_weighting > ${prefix}__\${bname}_std_per_point.json
        fi
    done

    #Bundle_Length_Stats
    if [[ "$length_stats" ]];
    then
        echo "Merging Bundle_Length_Stats"
        scil_json_merge_entries *_length.json ${prefix}__length_stats.json --add_parent_key ${prefix} \
                --keep_separate
    fi

    #Bundle_Endpoints_Map
    if [[ "$endpoints" ]];
    then
        echo "Merging Bundle_Endpoints_Map"
        scil_json_merge_entries *_endpoints_raw.json ${prefix}__endpoints_map_raw.json \
            --no_list --add_parent_key ${prefix}

        #Bundle_Metrics_Stats_In_Endpoints
        scil_json_merge_entries *_tail.json *_head.json ${prefix}__endpoints_metric_stats.json \
            --no_list --add_parent_key ${prefix}
    fi

    #Bundle_Mean_Std
    if [[ "$mean_std" ]];
    then
        echo "Merging Bundle_Mean_Std"
        scil_json_merge_entries *_std.json ${prefix}__mean_std_stats.json --no_list --add_parent_key ${prefix}
    fi

    #Bundle_Volume
    if [[ "$volume" ]];
    then
        echo "Merging Bundle_Volume"
        scil_json_merge_entries *_volume_stat.json ${prefix}__volume.json --no_list --add_parent_key ${prefix} --keep_separate

        if [[ "$lesions_stats" ]];
        then
            echo "Merging Lesions Stats"
            scil_json_merge_entries *_volume_lesions_stat.json ${prefix}__volume_lesions.json --no_list --add_parent_key ${prefix}
            scil_json_merge_entries *_streamline_count_lesions_stat.json ${prefix}__streamline_count_lesions.json \
                --no_list --add_parent_key ${prefix}
            scil_json_merge_entries ${prefix}__lesion_stats.json ${prefix}__lesion_stats.json \
                --remove_parent_key --add_parent_key ${prefix} -f
        fi
    fi

    #Bundle_Streamline_Count
    if [[ "$streamline_count" ]];
    then
        echo "Merging Bundle_Streamline_Count"
        scil_json_merge_entries *_streamlines.json ${prefix}__streamline_count.json --no_list \
            --add_parent_key ${prefix}
    fi

    #Bundle_Volume_Per_Label
    if [[ "$volume_per_labels" ]];
    then
        echo "Merging Bundle_Volume_Per_Label"
        scil_json_merge_entries *_volume_label.json ${prefix}__volume_per_label.json --no_list \
            --add_parent_key ${prefix}

        if [[ "$lesions_stats" ]];
        then
            echo "Merging Bundle_Volume_Per_Label in Lesions"
            scil_json_merge_entries *_volume_per_label_lesions_stat.json ${prefix}__volume_per_label_lesions.json \
                --no_list --add_parent_key ${prefix}
        fi
    fi

    #Bundle_Mean_Std_Per_Point
    if [[ "$mean_std_per_point" ]];
    then
        echo "Merging Bundle_Mean_Std_Per_Point"
        scil_json_merge_entries *_std_per_point.json ${prefix}__mean_std_per_point_stats.json --no_list \
            --add_parent_key ${prefix}
    fi

    # -------------------
    # Mean/Std Stats TSV
    # -------------------
    # Extract mean values for each metric across bundles
    if [[ "$mean_std" ]];
    then
        f="${prefix}__mean_std_stats.json"
        out="${prefix}__mean_std_stats.tsv"
        jq -r --arg sid "${prefix}" '
            # Get sample data
            .[\$sid] as \$s

            # Extract all metric names from nested objects, removing prefix
            | ( [ \$s | to_entries[]
                | select(.value | type == "object")
                | (.value | to_entries
                    | map(select(.value | type == "object" and has("mean"))
                        | (.key
                            | sub("^" + \$sid + "__"; "")
                            | if test("desc-") then capture("desc-(?<m>[^:]+)\$").m
                            elif test("_afd_metric\$") then "afd_fixel"
                            else . end
                        )
                    )
                )
            ] ) as \$ms

            # Get unique sorted list of metrics
            | (\$ms | add // [] | unique | sort) as \$metrics

            # Create header row
            | (["sample", "bundle"] + \$metrics) | @tsv,

            # Create data rows for each bundle
            (\$s | to_entries[]
                | select(.value | type == "object")
                | . as \$b
                # Extract mean values for each metric, removing prefix from keys
                | (\$b.value | to_entries
                    | map(select(.value | type == "object" and has("mean"))
                        | {((.key
                            | sub("^" + \$sid + "__"; "")
                            | if test("desc-") then capture("desc-(?<m>[^:]+)\$").m
                            elif test("_afd_metric\$") then "afd_fixel"
                            else . end
                        )): .value.mean}
                    )
                    | add // {}
                ) as \$mm
                # Build row: sample, cleaned bundle name, then metric values
                | [\$sid, (\$b.key | sub("^" + \$sid + "__"; "") | sub("_labels_uniformized\$"; "") | gsub("_cleaned"; ""))]
                + (\$metrics | map(\$mm[.] // ""))
                | @tsv
            )
        ' "\$f" > "\$out"
    fi

    # ----------------------------
    # Mean/Std Per Point Stats TSV
    # ----------------------------
    # Extract mean values per point for each metric across bundles
    if [[ "$mean_std_per_point" ]];
    then
        f="${prefix}__mean_std_per_point_stats.json"
        out="${prefix}__mean_std_per_point_stats.tsv"
        jq -r --arg sid "${prefix}" '
            # Get sample data
            .[\$sid] as \$s

            # Extract all metric names (handling both simple and per-point metrics), removing prefix
            | ( [ \$s | to_entries[]
                | select(.value | type == "object")
                | (.value | to_entries
                    | map(
                        if .value | has("mean") then
                            # Simple metric with direct mean
                            (.key
                                | sub("^" + \$sid + "__"; "")
                                | if test("desc-") then capture("desc-(?<metric>[^:]+)\$").metric
                                elif test("_afd_metric\$") then "afd_fixel"
                                else . end)
                        elif (.value | type == "object") then
                            # Per-point metric
                            (.key
                                | sub("^" + \$sid + "__"; "")
                                | if test("desc-") then capture("desc-(?<metric>[^:]+)\$").metric
                                elif test("_afd_metric\$") then "afd_fixel"
                                else . end)
                        else empty end
                    )
                )
            ] ) as \$metric_lists

            # Get unique sorted list of metrics
            | (\$metric_lists | add // [] | unique | sort) as \$metrics

            # Create header row with "points" column
            | (["sample", "bundle", "points"] + \$metrics) | @tsv,

            # Create data rows for each bundle and point
            (\$s | to_entries[] as \$b
                | select(\$b.value | type == "object")
                | (\$b.value | to_entries) as \$me

                # Get all point identifiers (labels) for this bundle
                | (\$me
                    | map(
                        if .value | has("mean") then []
                        elif (.value | type == "object") then (.value | keys)
                        elif (.value | type == "array") then
                            # Generate point labels (001, 002, ..., 100, 101, ...)
                            [ range(1; (.value | length) + 1)
                                | (if . < 10 then ("00" + tostring)
                                elif . < 100 then ("0" + tostring)
                                else tostring end)
                            ]
                        else [] end
                    )
                    | add // []
                    | unique | sort
                ) as \$points

                # Normalize metric names, removing prefix
                | (\$me
                    | map({
                        key: (.key
                            | sub("^" + \$sid + "__"; "")
                            | if test("desc-") then capture("desc-(?<metric>[^:]+)\$").metric
                            elif test("_afd_metric\$") then "afd_fixel"
                            else . end),
                        value: .value
                    })
                ) as \$mes

                # Create one row per point (or single row if no points)
                | (\$points[]? // [""]) as \$pt

                # Build metric value map for this point
                | (\$mes
                    | map(
                        if (.value | has("mean")) then
                            # Simple mean value
                            {(.key): .value.mean}
                        else
                            # Per-point mean value
                            {(.key): ((.value[\$pt] // (if (.value | type == "array")
                                then (.value[((\$pt | tonumber) - 1)] // null)
                                else null end)) // {} | .mean // "")}
                        end
                    )
                    | add
                ) as \$rowmap

                # Build row: sample, cleaned bundle name, point, then metric values
                | [\$sid,
                    (\$b.key | sub("^" + \$sid + "__"; "") | sub("_labels_uniformized\$"; "") | gsub("_cleaned"; "")),
                    (\$pt // "")]
                + (\$metrics | map(\$rowmap[.] // ""))
                | @tsv
            )
        ' "\$f" > "\$out"
    fi

    # ----------------
    # Volume Stats TSV
    # ----------------
    # Extract volume-related metrics for each bundle
    if [[ "$volume" ]];
    then
        f="${prefix}__volume.json"
        out="${prefix}__volume.tsv"
        jq -r --arg sid "${prefix}" '
            # Get sample data
            .[\$sid] as \$s

            # Collect bundles that are objects
            | (\$s | to_entries | map(select(.value | type == "object"))) as \$bundles

            # Get union of all metric keys across bundles (no prefix removal needed - values are plain)
            | (\$bundles | map(.value | keys) | add // [] | unique | sort) as \$metrics

            # Create header row
            | (["sample", "bundle"] + \$metrics) | @tsv,

            # Create data rows for each bundle
            (\$bundles[]
                | . as \$b
                | (\$b.value) as \$vals
                # Build row: sample, cleaned bundle name, then metric values
                | ([\$sid,
                    (\$b.key | sub("^" + \$sid + "__"; "") | gsub("_cleaned"; "")
                        | gsub("(_volume_stat|_labels_uniformized)\$"; ""))]
                    + (\$metrics | map((\$vals[.] // ""))))
                | @tsv
            )
        ' "\$f" > "\$out"
    fi

    # ----------------
    # Length Stats TSV
    # ----------------
    # Extract streamline length statistics for each bundle
    if [[ "$length_stats" ]];
    then
        f="${prefix}__length_stats.json"
        out="${prefix}__length.tsv"
        jq -r --arg sid "${prefix}" '
            # Get sample data
            .[\$sid] as \$s

            # Collect bundles that are objects
            | (\$s | to_entries | map(select(.value | type == "object"))) as \$bundles

            # Get union of metric keys (excluding internal data keys), removing prefix
            | (\$bundles
                | map(.value | keys - ["data_per_point_keys", "data_per_streamline_keys"] | map(sub("^" + \$sid + "__"; "")))
                | add // []
                | unique | sort) as \$metrics

            # Create header row
            | (["sample", "bundle"] + \$metrics) | @tsv,

            # Create data rows for each bundle
            (\$bundles[]
                | . as \$b
                | (\$b.value) as \$vals
                # Build metric map with cleaned keys
                | (\$vals | to_entries | map({key: (.key | sub("^" + \$sid + "__"; "")), value: .value}) | from_entries) as \$cleaned_vals
                # Build row: sample, cleaned bundle name, then metric values
                | ([\$sid,
                    (\$b.key | sub("^" + \$sid + "__"; "") | gsub("_cleaned"; "")
                        | gsub("(_volume_stat|_labels_uniformized|_length)\$"; ""))]
                    + (\$metrics | map((\$cleaned_vals[.] // ""))))
                | @tsv
            )
        ' "\$f" > "\$out"
    fi

    # --------------------------
    # Volume Per Label Stats TSV
    # --------------------------
    # Extract volume metrics per point/label for each bundle
    if [[ "$volume_per_labels" ]];
    then
        f="${prefix}__volume_per_label.json"
        out="${prefix}__volume_per_label.tsv"
        jq -r --arg sid "${prefix}" '
            # Get sample data
            .[\$sid] as \$s

            # Collect bundles that are objects
            | (\$s | to_entries | map(select(.value | type == "object"))) as \$bundles

            # Get union of metric keys (excluding internal data keys), removing prefix
            | (\$bundles
                | map(.value | keys - ["data_per_point_keys", "data_per_streamline_keys"] | map(sub("^" + \$sid + "__"; "")))
                | add // []
                | unique | sort) as \$metrics

            # Create header row with "points" column
            | (["sample", "bundle", "points"] + \$metrics) | @tsv,

            # Create data rows for each bundle (and points if available)
            (\$bundles[]
                | . as \$b
                | (\$b.value) as \$vals

                # Get all point identifiers for this bundle
                | (\$b.value | to_entries
                    | map(select(.value | type == "object") | .value | keys)
                    | add // []
                    | unique | sort) as \$bpoints

                # Build metric map with cleaned keys
                | (\$vals | to_entries | map({key: (.key | sub("^" + \$sid + "__"; "")), value: .value}) | from_entries) as \$cleaned_vals

                # If no points, create single row; otherwise one row per point
                | ((\$bpoints | length) as \$n
                    | if \$n == 0 then
                        # No points - single row with all values
                        ([\$sid,
                            (\$b.key | sub("^" + \$sid + "__"; "") | gsub("_cleaned"; "")
                                | gsub("(_volume_stat|_labels_uniformized|_length)\$"; "")),
                            ""]
                            + (\$metrics | map((\$cleaned_vals[.] // ""))))
                    else
                        # Multiple points - one row per point
                        (\$bpoints[]
                            | . as \$pt
                            | ([\$sid,
                                (\$b.key | sub("^" + \$sid + "__"; "") | gsub("_cleaned"; "")
                                    | gsub("(_volume_stat|_labels_uniformized|_length)\$"; "")),
                                \$pt]
                                + (\$metrics | map(
                                    if (\$cleaned_vals[.] | type) == "object"
                                    then (\$cleaned_vals[.][\$pt] // "")
                                    else (\$cleaned_vals[.] // "") end
                                )))
                        )
                    end
                )
                | @tsv
            )
        ' "\$f" > "\$out"
    fi

    # ------------------------------------
    # Streamline Count Lesions Stats TSV
    # ------------------------------------
    # Extract streamline count metrics per lesion for each bundle
    if [[ "$lesions_stats" ]];
    then
        f="${prefix}__streamline_count_lesions.json"
        out="${prefix}__streamline_count_lesions.tsv"
        jq -r --arg sid "${prefix}" '
            # Get sample data
            .[\$sid] as \$s

            # Collect bundles that are objects
            | (\$s | to_entries | map(select(.value | type == "object"))) as \$bundles

            # Get union of all metric keys across bundles, removing prefix
            | (\$bundles | map(.value | keys | map(sub("^" + \$sid + "__"; ""))) | add // [] | unique | sort) as \$metrics

            # Create header row with "lesion" column
            | (["sample", "bundle", "lesion"] + \$metrics) | @tsv,

            # Create data rows for each bundle and lesion
            (\$bundles[]
                | . as \$b
                | (\$b.value) as \$vals
                # Build metric map with cleaned keys
                | (\$vals | to_entries | map({key: (.key | sub("^" + \$sid + "__"; "")), value: .value}) | from_entries) as \$cleaned_vals
                # For each metric that has nested lesion data
                | (\$metrics[] as \$metric
                    | if (\$cleaned_vals[\$metric] | type) == "object" then
                        # Extract lesion IDs for this metric
                        (\$cleaned_vals[\$metric] | keys) as \$lesion_ids
                        | \$lesion_ids[] as \$lesion_id
                        | [\$sid,
                            (\$b.key | sub("^" + \$sid + "__"; "") | gsub("_cleaned"; "")
                                | gsub("(_streamline_count_lesions_stat|_volume_stat|_labels_uniformized)\$"; "")),
                            \$lesion_id,
                            \$cleaned_vals[\$metric][\$lesion_id]]
                    else
                        empty
                    end
                )
                | @tsv
            )
        ' "\$f" > "\$out"
    fi

    # ----------------------------------------
    # Volume Per Label Lesions Stats TSV
    # ----------------------------------------
    # Extract volume metrics per lesion for bundles
    if [[ "$lesions_stats" ]];
    then
        f="${prefix}__volume_per_label_lesions.json"
        out="${prefix}__volume_per_label_lesions.tsv"
        jq -r --arg sid "${prefix}" '
            # Get sample data
            .[\$sid] as \$s

            # Collect bundles that are objects
            | (\$s | to_entries | map(select(.value | type == "object"))) as \$bundles

            # Get union of metric keys, excluding lesion_volume, removing prefix
            | (\$bundles
                | map(.value | keys - ["data_per_point_keys", "data_per_streamline_keys", "lesion_volume"] | map(sub("^" + \$sid + "__"; "")))
                | add // []
                | unique | sort) as \$metrics

            # Create header row without "lesion" column
            | (["sample", "bundle", "points"] + \$metrics) | @tsv,

            # Create data rows for each bundle and point
            (\$bundles[]
                | . as \$b
                | (\$b.value) as \$vals

                # Build metric map with cleaned keys
                | (\$vals | to_entries | map({key: (.key | sub("^" + \$sid + "__"; "")), value: .value}) | from_entries) as \$cleaned_vals

                # Get all point labels (lesion IDs) for this bundle
                | ([
                    \$cleaned_vals | to_entries[]
                    | select(.value | type == "object")
                    | .value | keys
                    ] | add // [] | unique | sort) as \$point_ids

                # For each point (lesion ID)
                | \$point_ids[] as \$pt
                | [\$sid,
                    (\$b.key | sub("^" + \$sid + "__"; "") | gsub("_cleaned|_labels_map"; "")
                        | gsub("(_volume_per_label_lesions_stat|_volume_stat|_labels_uniformized|_length)\$"; "")),
                    \$pt]
                + (\$metrics | map(
                    if (\$cleaned_vals[.] | type) == "object" then
                        # Object format: metric -> point_id -> value
                        (\$cleaned_vals[.][\$pt] // "")
                    else
                        (\$cleaned_vals[.] // "")
                    end
                ))
                | @tsv
            )
        ' "\$f" > "\$out"
    fi

    # Clean up intermediate JSON files
    # rm *.json

    # ========================================
    # Merge TSV files
    # ========================================

    # Build list of available TSV files for merging
    tsv_files=()
    [[ -f "${prefix}__mean_std_stats.tsv" ]] && tsv_files+=("${prefix}__mean_std_stats.tsv")
    [[ -f "${prefix}__volume.tsv" ]] && tsv_files+=("${prefix}__volume.tsv")
    [[ -f "${prefix}__length.tsv" ]] && tsv_files+=("${prefix}__length.tsv")

    # Only merge if we have at least one TSV file
    if [[ \${#tsv_files[@]} -gt 0 ]]; then
        # Start with the first available file
        cp "\${tsv_files[0]}" tmp.tsv

        # Merge remaining files one by one
        for ((i=1; i<\${#tsv_files[@]}; i++)); do
            f1=tmp.tsv
            f2="\${tsv_files[i]}"
            out=tmp_merged.tsv

            h1=\$(head -n1 "\$f1")
            h2=\$(head -n1 "\$f2")

            # Combine headers: first 2 cols from f1, rest from f1, rest from f2
            printf '%s\\t%s\\t%s\\n' "\$(echo "\$h1" | cut -f1-2)" "\$(echo "\$h1" | cut -f3-)" "\$(echo "\$h2" | cut -f3-)" > "\$out"

            # Combine data rows
            paste <(tail -n +2 "\$f1" | cut -f1-2) <(tail -n +2 "\$f1" | cut -f3-) <(tail -n +2 "\$f2" | cut -f3-) >> "\$out"

            mv "\$out" tmp.tsv
        done

        # Rename final merged file
        mv tmp.tsv ${prefix}_desc-mean_stats.tsv
    fi

    # Merge per-point stats if both files exist
    if [[ -f "${prefix}__mean_std_per_point_stats.tsv" ]] && [[ -f "${prefix}__volume_per_label.tsv" ]]; then
        f1=${prefix}__mean_std_per_point_stats.tsv
        f2=${prefix}__volume_per_label.tsv
        out=${prefix}_desc-point_stats.tsv

        h1=\$(head -n1 "\$f1")
        h2=\$(head -n1 "\$f2")

        # Combine headers: first 3 cols from f1, rest from f1, rest from f2
        printf '%s\\t%s\\t%s\\n' "\$(echo "\$h1" | cut -f1-3)" "\$(echo "\$h1" | cut -f4-)" "\$(echo "\$h2" | cut -f4-)" > "\$out"

        # Combine data rows
        paste <(tail -n +2 "\$f1" | cut -f1-3) <(tail -n +2 "\$f1" | cut -f4-) <(tail -n +2 "\$f2" | cut -f4-) >> "\$out"
    elif [[ -f "${prefix}__mean_std_per_point_stats.tsv" ]]; then
        # Only mean_std_per_point exists
        cp ${prefix}__mean_std_per_point_stats.tsv ${prefix}_desc-point_stats.tsv
    elif [[ -f "${prefix}__volume_per_label.tsv" ]]; then
        # Only volume_per_label exists
        cp ${prefix}__volume_per_label.tsv ${prefix}_desc-point_stats.tsv
    fi

    if [[ -f "${prefix}__volume_per_label_lesions.tsv" ]]; then
        cp ${prefix}__volume_per_label_lesions.tsv ${prefix}_desc-mean_param-lesions_stats.tsv
    fi

    if [[ -f "${prefix}__streamline_count_lesions.tsv" ]]; then
        cp ${prefix}__streamline_count_lesions.tsv ${prefix}_desc-point_param-lesions_stats.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    scil_tractogram_print_info -h
    scil_bundle_compute_endpoints_map -h
    scil_volume_stats_in_ROI -h
    scil_bundle_mean_std -h
    scil_bundle_shape_measures -h
    scil_tractogram_count_streamlines -h
    scil_bundle_volume_per_label -h
    scil_bundle_mean_std -h
    scil_json_merge_entries -h

    touch ${prefix}__length_stats.json
    touch ${prefix}__endpoints_map_raw.json
    touch ${prefix}__endpoints_metric_stats.json
    touch ${prefix}__mean_std.json
    touch ${prefix}__volume.json
    touch ${prefix}__volume_lesions.json
    touch ${prefix}__streamline_count.json
    touch ${prefix}__streamline_count_lesions.json
    touch ${prefix}__volume_per_label.json
    touch ${prefix}__volume_per_label_lesions.json
    touch ${prefix}__mean_std_per_point.json
    touch ${prefix}_endpoints_map_head.nii.gz
    touch ${prefix}_endpoints_map_tail.nii.gz
    touch ${prefix}__lesion_stats.json
    touch ${prefix}_lesion_map.nii.gz
    touch ${prefix}_desc-mean_param-lesions_stats.tsv
    touch ${prefix}_desc-point_param-lesions_stats.tsv
    touch ${prefix}_desc-mean_stats.tsv
    touch ${prefix}_desc-point_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scilpy: \$(uv pip -q -n list | grep scilpy | tr -s ' ' | cut -d' ' -f2)
    END_VERSIONS
    """
}
