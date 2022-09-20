#!/bin/bash

#Ensure that output dirs are clear
rm tests/outputs/*/NC_045512/*

#Run the pipeline with a test case from gnomonicus unit tests for consistency
nextflow run . -profile docker --sample $(pwd)/tests/test-cases/NC_045512.2-double-minos.vcf --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv  --output_dir $(pwd)/tests/outputs/1

#Run a pipeline which should only produce a variants.csv
nextflow run . -profile docker --sample $(pwd)/tests/test-cases/NC_045512.2-just-variant.vcf --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv  --output_dir $(pwd)/tests/outputs/2

#Run a pipeline which should produce both a variants.csv and mutations.csv, but no effects.csv
nextflow run . -profile docker --sample $(pwd)/tests/test-cases/NC_045512.2-no-effects.vcf --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv  --output_dir $(pwd)/tests/outputs/3

#Run a pipeline without docker to ensure that this does not affect outputs
pip install gnomonicus
nextflow run . --sample $(pwd)/tests/test-cases/NC_045512.2-double-minos.vcf --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv  --output_dir $(pwd)/tests/outputs/4

#Install requirements and test with python
pip install pytest pandas recursive_diff
pytest -vv tests