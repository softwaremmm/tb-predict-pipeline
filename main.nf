#!/usr/bin/env nextflow

// default parameters
params.seq_platform = "illumina"

workflow {
    ANSI_GREEN = "\033[1;32m"
    ANSI_RESET = "\033[0m"

    //Setup so --help triggers the help message
    if (params.help) {
        log.info(
            """
            ========================================================================
            M Y C O B A C T E R I A L  P R E D I C T I O N  P I P E L I N E

            Utilises a minos VCF file to produce variations, mutations and
            drug resistance predictions based on provided a reference genome and a resistance catalogue.

            Mandatory parameters:
            ------------------------------------------------------------------------
            --sample                Path to the sample minos VCF
            --gvcf                  Path to the non-compressed gvcf file
            --reference             Path to the reference genome's genbank file, or a pickle dump of the corresponding gumpy Genome
            --catalogue             Path to the resistance catalogue
            --null_positions        Path to the null positions file
            --seq_platform          Sequencing platform used ('illumina' or 'ont'). Default is 'illumina'
            """.stripIndent()
        )
        exit(0)
    }


    //Log pre-run info
    log.info(
        """
        ========================================================================
        M Y C O B A C T E R I A L  P R E D I C T I O N  P I P E L I N E
        Parameters used:
        ------------------------------------------------------------------------
        --sample                ${params.sample}
        --gvcf                  ${params.gvcf}
        --reference             ${params.reference}
        --catalogue             ${params.catalogue}
        --null_positions        ${params.null_positions}
        --seq_platform          ${params.seq_platform}

        Runtime data:
        ------------------------------------------------------------------------
        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        """.stripIndent()
    )

    sample = Channel
        .fromPath("${params.sample}", checkIfExists: true)
        .map { it -> tuple(it.baseName, it) }
    gvcf = Channel.fromPath("${params.gvcf}", checkIfExists: true)
    input = sample.merge(gvcf)

    reference = Channel.fromPath(params.reference, checkIfExists: true).first()
    catalogue = Channel.fromPath(params.catalogue, checkIfExists: true).first()
    null_positions = Channel.fromPath(params.null_positions, checkIfExists: true).first()

    gnomonicus_workflow(input, params.seq_platform, reference, catalogue, null_positions)
}

workflow gnomonicus_workflow {
    take:
    samples // Channel of tuples with sample_id, vcf, gvcf
    seq_platform
    reference
    catalogue
    null_positions

    main:
    gnomonicus_json = runPrediction(samples, seq_platform, reference, catalogue, null_positions)

    emit:
    gnomonicus_json
}


workflow batch {
    ANSI_GREEN = "\033[1;32m"
    ANSI_RESET = "\033[0m"

    // Helper workflow for running a batch locally
    if (params.help) {
        log.info(
            """
            ========================================================================
            M Y C O B A C T E R I A L  P R E D I C T I O N  P I P E L I N E

            Utilises a minos VCF file to produce variations, mutations and
            drug resistance predictions based on provided a reference genome and a resistance catalogue.
            Expects each sample to have a corresponding directory with vcf and gvcf

            Mandatory parameters:
            ------------------------------------------------------------------------
            --samples               Path pattern for sample VCFs (e.g. 'samples/*/*.vcf')
            --gvcfs                 Path pattern for gvcf files (e.g. 'samples/*/*.gvcf')
            --reference             Path to the reference genome's genbank file, or a pickle dump of the corresponding gumpy Genome
            --catalogue             Path to the resistance catalogue
            --null_positions        Path to the null positions file
            --seq_platform          Sequencing platform used ('illumina' or 'ont'). Default is 'illumina'
            """.stripIndent()
        )
        exit(0)
    }


    //Log pre-run info
    log.info(
        """
        ========================================================================
        M Y C O B A C T E R I A L  P R E D I C T I O N  P I P E L I N E
        Parameters used:
        ------------------------------------------------------------------------
        --samples               ${params.samples}
        --gvcfs                 ${params.gvcfs}
        --reference             ${params.reference}
        --catalogue             ${params.catalogue}
        --null_positions        ${params.null_positions}
        --seq_platform          ${params.seq_platform}

        Runtime data:
        ------------------------------------------------------------------------
        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        """.stripIndent()
    )

    samples = Channel
        .fromPath("${params.samples}", checkIfExists: true, glob: true)
        .ifEmpty { error("cannot find any reads matching ${params.samples}") }
        .map { it -> tuple(it.parent.simpleName, it) }

    gvcfs = Channel
        .fromPath("${params.gvcfs}", checkIfExists: true, glob: true)
        .ifEmpty { error("cannot find any reads matching ${params.gvcfs}") }
        .map { it -> tuple(it.parent.simpleName, it) }

    input = samples.join(gvcfs)

    input.take(3).view()

    reference = Channel.fromPath(params.reference, checkIfExists: true).first()
    catalogue = Channel.fromPath(params.catalogue, checkIfExists: true).first()
    null_positions = Channel.fromPath(params.null_positions, checkIfExists: true).first()

    runPrediction(input, params.seq_platform, reference, catalogue, null_positions)
}


//Run gnomonicus
process runPrediction {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container "lhr.ocir.io/lrbvkel2wjot/oxfordmmm/gnomonicus:v3.0.9"
    cpus 2
    maxRetries 5
    memory {
        params.testing == "" ? 8.GB * (0.8 + (task.attempt / 5)) : "6GB"
    }

    pod label: "name", value: "tb-predict-pipeline:runPrediction"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path("variants.vcf"), path("all_rows.gvcf")
    val seq_platform
    path reference
    path catalogue
    path null_positions

    output:
    tuple val(sample_name), path("resistance_prediction_report.json")

    script:
    """
    merge-vcfs --minos_vcf variants.vcf --gvcf all_rows.gvcf --resistant-positions ${null_positions} --output "${sample_name}.vcf"

    if [ ${seq_platform} == 'illumina' ]
    then
        gnomonicus --genome_object ${reference} --catalogue ${catalogue} --vcf_file "${sample_name}.vcf" --json --output_dir . --resistance_genes --min_dp 3
    fi

    if [ ${seq_platform} == 'ont' ]
    then
        gnomonicus --genome_object ${reference} --catalogue ${catalogue} --vcf_file "${sample_name}.vcf" --json --output_dir . --resistance_genes --min_dp 5
    fi

    mv "${sample_name}.gnomonicus-out.json" resistance_prediction_report.json
    """

    stub:
    """
    touch resistance_prediction_report.json
    """
}
