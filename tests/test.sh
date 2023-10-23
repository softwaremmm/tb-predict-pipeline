#!/bin/bash
set -e

keepOutput(){
    local guid=$1
    #Copy the output to a known dir for testing later
    mkdir -p $(pwd)/tests/outputs/$guid

    #With a clean `work` dir, the only output file is going to be the only run
    cp $(pwd)/work/*/*/resistance_prediction_report.json $(pwd)/tests/outputs/$guid

    #Keep the dir clean so we know which files are from which inputs
    sudo rm -rf $(pwd)/work
}

#Ensure that output dirs are clear
rm -rf tests/outputs/*
sudo rm -rf work

#Setup virtualenv if required
if test -d gnomonicus_venv; then
    echo "Virtualenv already existed."
else
    echo "Creating virtualenv..."
    pip install virtualenv
    python3 -m virtualenv gnomonicus_venv
fi

source gnomonicus_venv/bin/activate

gnomonicus_version=$(cat main.nf | grep -E "container\ ?=\ ?" | cut -d ":" -f 2 | tr -d \")


#Match the gnomonicus version to the version used in the container
pip install gnomonicus==$gnomonicus_version

#Do some processing...

#These are made up test cases to hit edge cases
#Most of these use COVID-19 as it is quick to process, but some use a fake genome
sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --sample $(pwd)/tests/test-cases/NC_045512.2-S_E484K-minos.vcf --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv --minor_populations $(pwd)/tests/test-cases/no-minors.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_045512.fasta --testing true
keepOutput 1

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --sample $(pwd)/tests/test-cases/NC_045512.2-S_E484K-samtools.vcf --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv --minor_populations $(pwd)/tests/test-cases/no-minors.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_045512.fasta --testing true
keepOutput 2

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --sample $(pwd)/tests/test-cases/NC_045512.2-S_F2F-minos.vcf --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv --minor_populations $(pwd)/tests/test-cases/no-minors.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_045512.fasta --testing true
keepOutput 3

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --sample $(pwd)/tests/test-cases/NC_045512.2-S_F2L-minos.vcf --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv --minor_populations $(pwd)/tests/test-cases/no-minors.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_045512.fasta --testing true
keepOutput 4

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --sample $(pwd)/tests/test-cases/NC_045512.2-S_200_indel-minos.vcf --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv --minor_populations $(pwd)/tests/test-cases/no-minors.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_045512.fasta --testing true
keepOutput 5

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --sample $(pwd)/tests/test-cases/NC_045512.2-double-minos.vcf --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv --minor_populations $(pwd)/tests/test-cases/no-minors.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_045512.fasta --testing true
keepOutput 6

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --sample "$(pwd)/tests/test-cases/NC_045512.2-S_E484K&1450_ins_a-minos.vcf" --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv --minor_populations $(pwd)/tests/test-cases/no-minors.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_045512.fasta --testing true
keepOutput 7

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --sample $(pwd)/tests/test-cases/NC_045512.2-minors.vcf --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue.csv --minor_populations $(pwd)/tests/test-cases/minor_alleles.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_045512.fasta --testing true
keepOutput 8

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_045512.2.gbk --sample $(pwd)/tests/test-cases/NC_045512.2-minors.vcf --catalogue $(pwd)/tests/test-cases/NC_045512.2-test-catalogue-COV.csv --minor_populations $(pwd)/tests/test-cases/minor_alleles.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_045512.fasta --testing true
keepOutput 9

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/TEST-DNA.gbk --sample $(pwd)/tests/test-cases/TEST-DNA-large-del.vcf --catalogue $(pwd)/tests/test-cases/TEST-DNA-catalogue.csv --minor_populations $(pwd)/tests/test-cases/minor_alleles.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_045512.fasta --testing true
keepOutput 10

#These are realistic TB cases
#As they're TB, pre-pickle the genome to save time
gbkToPkl tests/test-cases/NC_000962.3.gbk --compress

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_000962.3.gbk.pkl --sample $(pwd)/tests/test-cases/NC_000962_3_test_0001.vcf --catalogue $(pwd)/tests/test-cases/NC_000962_3_catalogue_1.csv --minor_populations $(pwd)/tests/test-cases/minor_alleles.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_000962.3.fasta --testing true
keepOutput 11

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_000962.3.gbk.pkl --sample $(pwd)/tests/test-cases/NC_000962_3_test_0002.vcf --catalogue $(pwd)/tests/test-cases/NC_000962_3_catalogue_1.csv --minor_populations $(pwd)/tests/test-cases/minor_alleles.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_000962.3.fasta --testing true
keepOutput 12

sudo nextflow run . -profile docker --reference $(pwd)/tests/test-cases/NC_000962.3.gbk.pkl --sample $(pwd)/tests/test-cases/NC_000962_3_test_0003.vcf --catalogue $(pwd)/tests/test-cases/NC_000962_3_catalogue_1.csv --minor_populations $(pwd)/tests/test-cases/minor_alleles.txt --fasta_adjudication $(pwd)/tests/test-cases/NC_000962.3.fasta --testing true
keepOutput 13


#Install requirements and test with python
pip install pytest recursive_diff

pytest --exitfirst --verbose --failed-first -vv

