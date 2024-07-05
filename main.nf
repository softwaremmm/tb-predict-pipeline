#!/usr/bin/env nextflow

//Set DSL2 syntax
nextflow.enable.dsl=2

//Define ANSI colours for ease
ANSI_GREEN = "\033[1;32m"
ANSI_RESET = "\033[0m"

params.seq_platform = "illumina"

//Run gnomonicus
process runPrediction {
    container = "lhr.ocir.io/lrbvkel2wjot/oxfordmmm/gnomonicus:v2.6.6"
    cpus = 2
    maxRetries 5
    memory = { 
        params.testing=="" ? 16.GB * (0.8 + (task.attempt/5)) : "6GB"
    }

    debug true
    pod label: "name", value: "tb-predict-pipeline:runPrediction"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
        path sample
        path reference
        path catalogue
        path minor_populations
        path gvcf
        path null_positions
    output:
        path "resistance_prediction_report.json"
    script:
        """
        vcf_name=\$(basename $sample)
        guid=\${vcf_name%.vcf}

        if [ ${params.seq_platform} == 'illumina' ]
        then
            mkdir original
            mv $sample original/\$vcf_name.vcf
            merge-vcfs --minos_vcf original/\$vcf_name.vcf --gvcf $gvcf --resistant-positions $null_positions --output $sample

            gnomonicus --genome_object $reference --catalogue $catalogue --vcf_file $sample --json --output_dir . --minor_populations $minor_populations --resistance_genes --min_dp 2
        fi

        if [ ${params.seq_platform} == 'ont' ]
        then
            mkdir original
            mv $sample original/\$vcf_name.vcf
            merge-vcfs --minos_vcf original/\$vcf_name.vcf --gvcf $gvcf --resistant-positions $null_positions --output $sample

            gnomonicus --genome_object $reference --catalogue $catalogue --vcf_file $sample --json --output_dir . --minor_populations $minor_populations --resistance_genes --min_dp 5
        fi


        #Get the name of the output JSON to move it to `resistance_prediction_report.json`
        mv \$guid.gnomonicus-out.json resistance_prediction_report.json

        """
    stub:
        """
        touch resistance_prediction_report.json
        """
}

workflow gnomonicus_workflow {
    take:
        sample
        reference
        catalogue
        minor_populations
        gvcf
        null_positions

    main:
        gnomonicus_json = runPrediction(sample, reference, catalogue, minor_populations, gvcf, null_positions)

    emit:
        gnomonicus_json
}

workflow {
    main:
        //Setup so --help triggers the help message
        if (params.help) {
            log.info """
            ========================================================================
            M Y C O B A C T E R I A L  P R E D I C T I O N  P I P E L I N E
            
            Utilises a minos VCF file to produce variations, mutations and
            drug resistance predictions based on provided a reference genome and a resistance catalogue.
                
            Mandatory parameters:
            ------------------------------------------------------------------------
            --sample                Path to the sample minos VCF
            --reference             Path to the reference genome's genbank file, or a pickle dump of the corresponding gumpy Genome
            --catalogue             Path to the resistance catalogue
            --minor_populations     Path to a line separated file of genome indices to check for minor populations
            --gvcf                  Path to the non-compressed gvcf file
            --null_positions        Path to the null positions file
            --seq_platform          Sequencing platform used ('illumina' or 'ont'). Default is 'illumina'
            """
            .stripIndent()
            exit(0)
        }


        //Log pre-run info
        log.info """
        ========================================================================
        M Y C O B A C T E R I A L  P R E D I C T I O N  P I P E L I N E
        Parameters used:
        ------------------------------------------------------------------------
        --sample                ${params.sample}
        --reference             ${params.reference}
        --catalogue             ${params.catalogue}
        --minor_populations     ${params.minor_populations}
        --gvcf                  ${params.gvcf}
        --null_positions        ${params.null_positions}
        --seq_platform          ${params.seq_platform}

        Runtime data:
        ------------------------------------------------------------------------
        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        """
        .stripIndent()

        gnomonicus_workflow(params.sample, params.reference, params.catalogue, params.minor_populations, params.gvcf, params.null_positions)
}
