#!/usr/bin/env nextflow

//Set DSL2 syntax
nextflow.enable.dsl=2

//Define ANSI colours for ease
ANSI_GREEN = "\033[1;32m"
ANSI_RESET = "\033[0m"

//Run gnomonicus
process runPrediction {
    container = "oxfordmmm/gnomonicus:latest"
    cpus = 2
    memory = "8GB"
    input:
        path sample
        path reference
        path catalogue
        path minor_populations
        val sample_name
    output:
        path "gnomonicus.json"
    script:
        """
        gnomonicus --genome_object $reference --catalogue $catalogue --vcf_file $sample --json --output_dir . --minor_populations $minor_populations
        
        #Get the name of the output JSON to move it to `gnomonicus.json`
        vcf_name=\$(basename $sample)
        guid=\${vcf_name%.vcf}
        mv \$guid.gnomonicus-out.json gnomonicus.json
        """
    stub:
        """
        touch gnomonicus.json
        """
}

workflow gnomonicus_workflow {
    take:
        sample
        reference
        catalogue
        minor_populations

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
        --sample            ${sample}
        --reference         ${reference}
        --catalogue         ${catalogue}
        --minor_populations ${minor_populations}

        Runtime data:
        ------------------------------------------------------------------------
        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        """
        .stripIndent()

        //Pull out the sample name from the vcf param
        sample_name = file(sample).simpleName

        gnomonicus_json = runPrediction(sample, reference, catalogue, minor_populations, sample_name)

    emit:
        gnomonicus_json
}

workflow {
    main:
        gnomonicus_workflow(params.sample, params.reference, params.catalogue, params.minor_populations)
}