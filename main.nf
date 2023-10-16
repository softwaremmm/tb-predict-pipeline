#!/usr/bin/env nextflow

//Set DSL2 syntax
nextflow.enable.dsl=2

//Define ANSI colours for ease
ANSI_GREEN = "\033[1;32m"
ANSI_RESET = "\033[0m"

//Run gnomonicus
process runPrediction {
    container = "lhr.ocir.io/lrbvkel2wjot/oxfordmmm/gnomonicus:v2.2.0"
    cpus = 2
    maxRetries 5
    memory = { 
        params.testing=="" ? 8.GB * (0.8 + (task.attempt/5)) : "6GB"
    }
    input:
        path sample
        path reference
        path catalogue
        path minor_populations
    output:
        path "resistance_prediction_report.json"
    script:
        """
        if [ ${workflow.profile} == 'kubernetes' ]
        then
            echo "Running with kubernetes"
            /bin/bash ${projectDir}/lib/s3fs_setup.sh $WORKSPACE
        fi

        gnomonicus --genome_object $reference --catalogue $catalogue --vcf_file $sample --json --output_dir . --minor_populations $minor_populations --resistance_genes
        
        #Get the name of the output JSON to move it to `resistance_prediction_report.json`
        vcf_name=\$(basename $sample)
        guid=\${vcf_name%.vcf}
        mv \$guid.gnomonicus-out.json resistance_prediction_report.json

        if [ ${workflow.profile} == 'kubernetes' ]
        then
            /bin/bash ${projectDir}/lib/s3fs_teardown.sh
        fi
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

    main:
        gnomonicus_json = runPrediction(sample, reference, catalogue, minor_populations)

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
            --sample            Path to the sample minos VCF
            --reference         Path to the reference genome's genbank file, or a pickle dump of the corresponding gumpy Genome
            --catalogue         Path to the resistance catalogue
            --minor_populations Path to a line separated file of genome indices to check for minor populations
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
        --sample            ${params.sample}
        --reference         ${params.reference}
        --catalogue         ${params.catalogue}
        --minor_populations ${params.minor_populations}

        Runtime data:
        ------------------------------------------------------------------------
        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        """
        .stripIndent()

        gnomonicus_workflow(params.sample, params.reference, params.catalogue, params.minor_populations)
}
