"""Given a fasta file, build a null gvcf file with no variants.
Used to provide gvcf files for testing as they can be too large to track with git.
"""
import argparse


def fasta_to_null_gvcf(fasta_path: str) -> list[str]:
    """Given a fasta file, produce a gvcf with a null call at each base.

    Args:
        fasta_path (str): Path to a fasta file.

    Returns:
        list[str]: List of dummy gvcf rows.
    """
    with open(fasta_path) as f:
        fasta = "".join([line.strip() for line in f if line[0] != ">"])
    
    header = [
        "##fileformat=VCFv4.2",
        "##source=clockwork merge samtools gvcf and minos vcf",
        "##fileDate=2024-06-18",
        '##ALT=<ID=*,Description="Represents allele(s) other than observed.">',
        '##FILTER=<ID=MAX_DP,Description="Maximum DP of 71.8895025378336 (= 3.0 standard deviations from the mean read depth 21.108)">',
        '##FILTER=<ID=MIN_DP,Description="Minimum DP of 2">',
        '##FILTER=<ID=MIN_FRS,Description="Minimum FRS of 0.9">',
        '##FILTER=<ID=MIN_GCP,Description="Minimum GT_CONF_PERCENTILE of 0.5">',
        '##FILTER=<ID=NO_DATA,Description="No information from minos or samtools">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FORMAT=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">',
        '##FORMAT=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">',
        '##FORMAT=<ID=AF2,Number=1,Type=Float,Description="Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)">',
        '##FORMAT=<ID=ALLELE_DP,Number=R,Type=Float,Description="Mean read depth of ref and each allele">',
        '##FORMAT=<ID=BQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)">',
        '##FORMAT=<ID=COV,Number=R,Type=Integer,Description="Number of reads on ref and alt alleles">',
        '##FORMAT=<ID=COV_TOTAL,Number=1,Type=Integer,Description="Total reads mapped at this site, from gramtools">',
        '##FORMAT=<ID=DP,Number=1,Type=Float,Description="Mean read depth of called allele (ie the ALLELE_DP of the called allele)">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">',
        '##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">',
        '##FORMAT=<ID=DP_ACGT,Number=8,Type=Integer,Description="Number of A-forward, A-reverse, C-forward, C-reverse, G-forward, G-reverse, T-forward, T-reverse bases">',
        '##FORMAT=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">',
        '##FORMAT=<ID=FRS,Number=1,Type=Float,Description="Fraction of reads that support the genotype call">',
        '##FORMAT=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias">',
        '##FORMAT=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence. Difference in log likelihood of most likely and next most likely genotype">',
        '##FORMAT=<ID=GT_CONF_PERCENTILE,Number=1,Type=Float,Description="Percentile of GT_CONF">',
        '##FORMAT=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">',
        '##FORMAT=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">',
        '##FORMAT=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">',
        '##FORMAT=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">',
        '##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">',
        '##FORMAT=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">',
        '##FORMAT=<ID=MQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)">',
        '##FORMAT=<ID=MQSBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">',
        '##FORMAT=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">',
        '##FORMAT=<ID=RPBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)">',
        '##FORMAT=<ID=SCBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)">',
        '##FORMAT=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">',
        '##FORMAT=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">',
        '##INFO=<ID=CALLER,Number=1,Description="Origin of call, one of minos, samtools, or none if there was no depth">',
        "##bcftoolsCommand=mpileup -I --output-type u -f /workspace/project/16302/myco_workflow/work/91/8116d243d1abd90d38c3ee90b6f49b/Ref_prepare/ref.fa outdir/map.bam",
        "##bcftoolsVersion=1.15.1+htslib-1.15.1",
        "##bcftools_callCommand=call -c -O v -o outdir/samtools.gvcf; Date=Tue Jun 18 12:41:58 2024",
        "##bcftools_callVersion=1.15.1+htslib-1.15.1",
        "##contig=<ID=NC_000962.3,length=4411532>",
        "##minosMeanReadDepth=21.108",
        "##minosReadDepthVariance=286.529",
        "##reference=file:///workspace/project/16302/myco_workflow/work/91/8116d243d1abd90d38c3ee90b6f49b/Ref_prepare/ref.fa",
        "##source=minos, version 0.12.5",
        "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample",
    ]
    rows = []
    for pos, char in enumerate(fasta):
        rows.append(f"NC_000962.3\t{pos+1}\t.\tC\t.\t.\tPASS\tCALLER=none\tGT:DP\t./.:0")
    return header + rows

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_path", help="Path to a fasta file.")
    parser.add_argument("--output_path", help="Path to write the gvcf file.")
    args = parser.parse_args()

    gvcf_rows = fasta_to_null_gvcf(args.fasta_path)
    with open(args.output_path, "w") as f:
        f.write("\n".join(gvcf_rows))


