#!/usr/bin/env nextflow

//Pipeline based on Lodestone tb pipeline (https://github.com/Pathogen-Genomics-Cymru/tb-pipeline)
//And SARS-CoV2_worflows (https://github.com/oxfordmmm/SARS-CoV2_workflows)

//Set DSL2 syntax
nextflow.enable.dsl=2

//Define ANSI colours for ease
ANSI_GREEN = "\033[1;32m"
ANSI_RESET = "\033[0m"


//Setup so --help triggers the help message
if (params.help) {
    helpMessage()
    exit(0)
}

def helpMessage() {
log.info """
========================================================================
M Y C O B A C T E R I A L  P R E D I C T I O N  P I P E L I N E
  
Utilises the output of the Lodestone pipeline (namely the minos VCF file) to produce variations, mutations and
drug resistance predictions based on provided a reference genome and a resistance catalogue.
	
Mandatory parameters:
------------------------------------------------------------------------
--sample            Path to the sample minos VCF
--reference         Path to the reference genome's genbank file, or a pickle dump of the corresponding gumpy Genome
--catalogue         Path to the resistance catalogue
--output_dir        Desired output path for all files produced by gnomonicus
--fasta             The kind of fasta file to generate. Defaults to fixed length (indels do not change length).
                    One of `fixed` and `variable`
"""
.stripIndent()
}


//Log pre-run info
log.info """
========================================================================
M Y C O B A C T E R I A L  P R E D I C T I O N  P I P E L I N E
Parameters used:
------------------------------------------------------------------------
--sample            ${params.sample}
--reference         ${params.reference}
--catalogue         ${params.catalogue}
--output_dir        ${params.output_dir}
--fasta             ${params.fasta}

Runtime data:
------------------------------------------------------------------------
Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
"""
.stripIndent()

//Run gnomonicus
process runPrediction {

    tag {sample_name}
    publishDir "${params.output_dir}/${sample_name}", mode: 'copy', pattern: '*', overwrite: 'true'
    input:
        path sample
        path reference
        path catalogue
        path outputDir
        val fasta
        val sample_name
    output:
        path "${sample_name}.gnomonicus.log"
        path "${sample_name}.variants.csv"
        path "${sample_name}.gnomonicus-out.json" //Always create the JSON
        path "${sample_name}.mutations.csv" optional true
        path "${sample_name}.effects.csv" optional true
        //One of these will always be created. Default is fixed length
        path "*-fixed.fasta" optional true
        path "*-variable.fasta" optional true
    script:
        """
        gnomonicus --genome_object $reference --catalogue $catalogue --vcf_file $sample --json --fasta $fasta --output_dir .
        """
}

workflow {
    main:
        //Pull out the sample name from the vcf param
        sample_name = file(params.sample).simpleName

        runPrediction(params.sample, params.reference, params.catalogue, params.output_dir, params.fasta, sample_name)
}