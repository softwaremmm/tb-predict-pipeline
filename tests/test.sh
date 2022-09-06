#!/bin/bash

#Run the pipeline with a test case from gnomon unit tests for consistency
nextflow run . -profile docker --sample $(pwd)/tests/test-cases/NC_045512.2-double-minos.vcf --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv  --output_dir $(pwd)/tests/outputs/1

#Run a pipeline which should only produce a variants.csv
nextflow run . -profile docker --sample $(pwd)/tests/test-cases/NC_045512.2-just-variant.vcf --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv  --output_dir $(pwd)/tests/outputs/2



#Install requirements and test with python
pip install pytest pandas recursive_diff
pytest -vv tests