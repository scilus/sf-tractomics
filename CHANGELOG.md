# scilus/sf-tractomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of scilus/sf-tractomics, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Default configuration shipped with each nf-neuro module (`modules/nf-neuro/<tool>/<sub>/nextflow.config`, selected with a module label) and each nf-neuro subworkflow (`subworkflows/nf-neuro/<name>/nextflow.config`: default params and the `ext` bindings of its modules). The pipeline includes them, then overrides the params it needs.
- Params `run_synthseg`, `run_bundleparc` and `atlas_iit_smooth_sigma`.

### `Changed`

- `ext` values are now closures (`{ params.x }`), so params set in a profile or a custom config reach the processes.
- The nf-neuro subworkflows no longer take an `options` map; they read their params, whose defaults live in their `nextflow.config`.

### `Fixed`

- `--run_t1_resampling` was ignored: T1 resampling never ran.
- Params set in a profile (e.g. `test`, `no_crop_no_mask`) did not reach the processes.
- The template registration outputs were suffixed with `null` instead of the template name.

### `Dependencies`

### `Deprecated`
