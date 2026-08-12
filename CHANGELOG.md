## New

* Support gzipped input files to process

## 1.6.1
* Use image tag `v3.0.11-1` to account for image sync failure.

## 1.6.0
* Rename output `merged.vcf` to `final.vcf`.
* Rename internally used terminology from all_rows to all_calls.

## 1.5.3
* Update gnomonicus from `v3.0.10` to `v3.0.11`
    * Arbitrary frameshifts are now properly supported. e.g a gene with a range of join(1..30, 40..50)
    * Genes crossing the genome boundary are now supported
    *   Detection of a coding gene without % 3 number of bases no longer causes panic. Instead, such genes are marked as incomplete (via gene name becoming INCOMPLETE_<gene>), and have a warning printed to stderr upon parsing

No effect on TB - these fixes enable NTM genbank files to be used.

## 1.5.2
* Output `merged.vcf` (which will become `final.vcf` in a future version)

gnomonicus:
* Only merge null calls from `alternate.gvcf` (Illumina / Clockwork) or `final.full.vcf` (ONT / Rundial)
into `minos.vcf` (Illumina / Clockwork) or `final.vcf` (ONT / Rundial). Merges take place at positions
set by the catalogue.

## 1.5.1

* chore: Parameterise container prefix

## 1.5.0
Lots of edge cases which we found during the ESCMID runs

gnomonicus:
* Various fixes around always producing CSVs when empty. Not relevant here

piezo:
* Ensure a deletion starting in the promoter is interpreted as a frameshift if it is actually a frameshift. e.g gene@-2_del_10 is equivalent to gene@-2_del_2 as well as gene@1_fs,
* Ensure null calls can't hit any wildcard rules

grumpy:
* Ensure SNPs after deletions are correctly parsed from VCF rows
* Correctly place minor SNPs at positions we have major deletions (and vice versa)
* Correct filter/null behaviour

- Rename input files in runPrediction to simplify code

## 1.4.0

Use Nextflow linting via Nextflow language server.
Improve local running

### Fix

- Make seq_platform a process input rather than use params.

### Chore

- use publish_dir
- simplify nextflow configs and readme
- add examples
- allow reference data to be given as relative paths

## 1.3.3 (2024-11-13)

### Fix

- add local profile
- remove debug and errorStrategy

## 1.3.2 (2024-10-14)

### Fix

- gnomonicus 3.0.1

## 1.3.1 (2024-09-16)

### Fix

- explicity retry on error

## 1.3.0 (2024-08-29)

### Feat

- gnomonicus 3.0.0

## 1.2.13 (2024-07-25)

### Fix

- gnomonicus 2.6.13

## 1.2.12 (2024-07-24)

### Fix

- gnomonicus 2.6.12

## 1.2.11 (2024-07-24)

### Fix

- gnomonicus 2.6.11

## 1.2.10 (2024-07-23)

### Fix

- gnomonicus 2.6.10

## 1.2.9 (2024-07-19)

### Fix

- gnomonicus 2.6.9

## 1.2.8 (2024-07-10)

### Fix

- gnomonicus-2.6.8

## 1.2.7 (2024-07-10)

### Fix

- use min_dp 3 for illumina

## 1.2.6 (2024-07-09)

### Fix

- refactor input channels for running batches locally

## 1.2.5 (2024-07-08)

### Fix

- gnomonicus 2.6.7

## 1.2.4 (2024-07-05)

### Fix

- gnomonicus-2.6.6

## 1.2.3 (2024-07-03)

### Fix

- use min_dp of 5 for ont
- pass gvcf for ont

## 1.2.2 (2024-06-28)

### Fix

- gnomonicus 2.6.5

## 1.2.1 (2024-06-24)

### Fix

- use gnomonicus 2.6.4

## 1.2.0 (2024-06-21)

### Feat

- add support gvcf mediated null calls

### Fix

- track empty gvcf for tests

## 1.1.8 (2024-06-10)

### Fix

- gnomonicus 2.5.7

## 1.1.7 (2024-05-09)

### Fix

- gnom to 2.5.6

## 1.1.6 (2024-03-27)

### Fix

- use gnomonicus 2.5.3

## 1.1.5 (2024-03-18)

### Fix

- bump gnomonicus to 2.5.2

## 1.1.4 (2024-03-14)

### Fix

- use gnomonicus 2.5.1

## 1.1.3 (2024-03-11)

### Fix

- use gnomonicus 2.5.0

## 1.1.2 (2024-01-30)

### Fix

- use 16gb ram

## 1.1.1 (2024-01-22)

### Fix

- use gnomonicus v2.3.7

## 1.1.0 (2024-01-17)

### Feat

- use commitizen to manage versions and release

### Fix

- adds pre-commit hook for conventional commits
- remove integrate subworkflows action
- **test**: add docker login
- **s3fs**: ensure buckets are unmounted on exit

## v1.0.4 (2023-11-09)

## v1.0.3 (2023-11-03)

## v1.0.2 (2023-11-01)

### Feat

- added pod labels to processes

## v1.0.1 (2023-10-25)

### Fix

- bump container version

## v1.0.0 (2023-10-24)

## v0.3.3 (2023-10-16)

## v0.3.2 (2023-10-05)

## v0.3.1 (2023-09-28)

## v0.3.0 (2023-08-24)

## v0.2.4 (2023-08-16)

## v0.2.3 (2023-08-14)

## v0.2.2 (2023-07-28)

## v0.2.1 (2023-07-27)

## v0.2.0 (2023-07-24)

## 0.1.0 (2023-07-14)
