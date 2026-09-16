[![Tests](https://github.com/oxfordmmm/tb-predict-pipeline/actions/workflows/tests.yaml/badge.svg)](https://github.com/oxfordmmm/tb-predict-pipeline/actions/workflows/tests.yaml)
# TB-Predict-Pipeline
Nextflow pipeline for producing variants, mutations and effects of a specified (minos) VCF file

## Requirements
- [Docker](https://docs.docker.com/get-docker/)
- Nextflow


## Running the NextFlow
Workflow takes parameters:
- seq_platform. `ont` or `illumina` (`illumina` by default)
- samples. path to vcf files, likely from minos
- gvcfs. Path to gvcfs
- reference. reference gbk file
- catalogue. Path to mutations catalogue. Must make "EMPTY_catalogue.csv" to not use a catalogue
- null_positions. Path to list of null positions. Must be provided but can be empty file.

To save output files need to set `--publish_dir` which will save output files to directory provided.

Example running locally:
```bash
nextflow run . -profile local --publish_dir results --seq_platform illumina \
    --samples tests/test-cases/NC_045512.2-S_E484K-minos.vcf \
    --gvcfs tests/test-cases/empty.gvcf \
    --reference tests/test-cases/NC_045512.2.gbk \
    --catalogue tests/test-cases/NC_045512.2-test-catalogue.csv \
    --null_positions tests/test-cases/no-null-positions.txt
```

### Batch locally
To run a batch of files locally though this pipeline requires input files to have the same root in order to be matched. For example:
- sample1_minos.vcf
- sample1_minos.gvcf
- sample2_minos.vcf
- sample2_minos.gvcf

Then provide globs for the samples and gvcfs.
All samples must use the same reference data.

```bash
nextflow run . -entry batch -profile local --publish_dir results --seq_platform illumina \
    --samples "tests/batch_dir/*.vcf" \
    --gvcfs "tests/batch_dir/*.gvcf" \
    --reference tests/test-cases/NC_045512.2.gbk \
    --catalogue tests/test-cases/NC_045512.2-test-catalogue.csv \
    --null_positions tests/test-cases/no-null-positions.txt
```

## Tests
As this is a very simple pipeline which just calls `gnomonicus`, the majority of testing occurs within `gnomonicus`. However, the included tests cover some edge cases of files being produced, as well as ensuring files are placed in correct places.
### Run the tests
Pipelines and evaluation are run using a bash script. This should also be automatically run on each push to the `main` branch through a CI action. Note that this repo intentionally does not use `nf-test` for these tests. This is due to a few reasons, but mostly that direct file comparisons are not reliable here due to fields such as timestamp
```bash
tests/test.sh
```

## Tags, Releases, and Committing
Use conventional commits. This is enforced with commitizen validate action and pre-commit hooks:
```bash
pre-commit install
```

This repo uses a standard gitflow approach, so changes should be first merged into develop and then released to main.
- In the develop branch semantic versioning is not used. Instead you can reference the commit hash to use it in a workflow.
- In a release branch you can create a release candidate with `cz bump a.b.c-rcX`. This also creates a tag.
- When release branch is ready for main run `cz bump a.b.c --files-only`. Manually write a human descriptive changelog. Then push these changes to main and make a release/tag there.
