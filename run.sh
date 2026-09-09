
knowledge=/mnt/block_data/knowledge_replica/knowledge
catalogues=${knowledge}/tuberculosis_amr_catalogues/catalogues/NC_000962.3/

# dir=/home/ubuntu/studies/ouh_mgit_study/data/results/illumina
# nextflow run . -entry batch -profile local --publish_dir ${dir} --seq_platform illumina -resume \
#     --samples "${dir}/*final.vcf" \
#     --gvcfs "${dir}/*alternate.gvcf" \
#     --reference ${catalogues}/NC_000962.3.gbk \
#     --catalogue ${catalogues}/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv \
#     --null_positions ${knowledge}/null_positions.txt


dir=/home/ubuntu/studies/ouh_mgit_study/data/results/ont_updated
nextflow run . -entry batch -profile local --publish_dir ${dir} --seq_platform ont -resume \
    --samples "${dir}/*variants.vcf" \
    --gvcfs "${dir}/*all_calls.vcf" \
    --reference ${catalogues}/NC_000962.3.gbk \
    --catalogue ${catalogues}/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv \
    --null_positions ${knowledge}/null_positions.txt


# dir=/home/ubuntu/studies/ont_myco_datasets/data/results/ont
# nextflow run . -entry batch -profile local --publish_dir ${dir} --seq_platform ont -resume \
#     --samples "${dir}/*final.vcf" \
#     --gvcfs "${dir}/*final.full.vcf" \
#     --reference ${catalogues}/NC_000962.3.gbk \
#     --catalogue ${catalogues}/NC_000962.3_WHO-UCN-TB-2023.5_v2.1_GARC1_RFUS.csv \
#     --null_positions ${knowledge}/null_positions.txt

###


# samples = channel
#     .fromPath("${params.samples}", checkIfExists: true, glob: true)
#     .ifEmpty { error("cannot find any reads matching ${params.samples}") }
#     .map { it -> tuple(it.getName().replaceFirst(".final.vcf", "").replaceFirst(".final.full.vcf", "").replaceFirst(".variants.vcf", ""), it) }


# gvcfs = channel
#     .fromPath("${params.gvcfs}", checkIfExists: true, glob: true)
#     .ifEmpty { error("cannot find any reads matching ${params.gvcfs}") }
#     .map { it -> tuple(it.getName().replaceFirst(".alternate.gvcf", "").replaceFirst(".final.full.vcf", "").replaceFirst(".all_calls.vcf", ""), it) }
