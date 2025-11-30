# scilus/sf-tractomics

[![GitHub Actions CI Status](https://github.com/scilus/sf-tractomics/actions/workflows/ci.yml/badge.svg)](https://github.com/scilus/sf-tractomics/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/scilus/sf-tractomics/actions/workflows/linting.yml/badge.svg)](https://github.com/scilus/sf-tractomics/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Deploy documentation](https://github.com/scilus/sf-tractomics/actions/workflows/deploy.yml/badge.svg?branch=main)](https://github.com/scilus/sf-tractomics/actions/workflows/deploy.yml)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.04.6-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/scilus/sf-tractomics)

## Introduction

**scilus/sf-tractomics** is a medical imaging pipeline that processes diffusion MRI images of human brains and reconstructs white matter pathways using tractography.

<!-- TODO nf-neuro : create a tube map for the pipeline -->

1. Read Inputs
2. Preprocess DWI
3. Preprocess T1
4. Segment T1
5. Reconstruct diffusion profiles
6. Run tractography

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

You can run the pipeline using:

```bash
nextflow run scilus/sf-tractomics \
   -profile <docker/singularity/.../institute> \
   --input <samplesheet.csv|/path/to/bids> \
   --outdir <OUTDIR>
```

Refer to [the usage](./docs/usage.md) to see how to format your samplesheet or for extended informations.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

scilus/sf-tractomics is a port by @AlexVCaron of tractoflow, originally written by @GuillaumeTh.

We thank the following people for their extensive assistance in the development of this pipeline:

- @jchoude
- @frheault
- @arnaudbore
- @mdesco

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
