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

    // replace is to catch both .vcf and .vcf.gz files
    sample = channel.fromPath("${params.sample}", checkIfExists: true)
        .map { it -> tuple(it.baseName.replace(".vcf", "").replace(".gvcf", ""), it) }
    gvcf = channel.fromPath("${params.gvcf}", checkIfExists: true)
    input = sample.merge(gvcf)
    input.view()

    reference = channel.fromPath(params.reference, checkIfExists: true).first()
    catalogue = channel.fromPath(params.catalogue, checkIfExists: true).first()
    null_positions = channel.fromPath(params.null_positions, checkIfExists: true).first()

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
    gnomonicus_out = runPrediction(samples, seq_platform, reference, catalogue, null_positions)

    emit:
    gnomonicus_json = gnomonicus_out.json
    gnomonicus_vcf = gnomonicus_out.vcf
    gnomonicus_variants = gnomonicus_out.variants_table
    gnomonicus_mutations = gnomonicus_out.mutations_table
    gnomonicus_effects = gnomonicus_out.effects_table
    gnomonicus_predictions = gnomonicus_out.predictions_table
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

    // replace is to catch both .vcf and .vcf.gz files
    samples = channel.fromPath("${params.samples}", checkIfExists: true, glob: true)
        .ifEmpty { error("cannot find any reads matching ${params.samples}") }
        .map { it -> tuple(it.baseName.replace(".vcf", "").replace(".gvcf", ""), it) }

    gvcfs = channel.fromPath("${params.gvcfs}", checkIfExists: true, glob: true)
        .ifEmpty { error("cannot find any reads matching ${params.gvcfs}") }
        .map { it -> tuple(it.baseName.replace(".vcf", "").replace(".gvcf", ""), it) }

    input = samples.join(gvcfs)

    input.take(3).view()

    reference = channel.fromPath(params.reference, checkIfExists: true).first()
    catalogue = channel.fromPath(params.catalogue, checkIfExists: true).first()
    null_positions = channel.fromPath(params.null_positions, checkIfExists: true).first()

    runPrediction(input, params.seq_platform, reference, catalogue, null_positions)
}


//Run gnomonicus
process runPrediction {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container params.container_prefix + "/oxfordmmm/gnomonicus:3.1.5"
    cpus 2
    maxRetries 5
    memory { params.testing == "" ? 8.GB + (4.GB * (task.attempt - 1)) : "6GB" }

    pod label: "name", value: "tb-predict-pipeline:runPrediction"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path("variants.vcf"), path("all_calls.vcf")
    val seq_platform
    path reference
    path catalogue
    path null_positions

    output:
    tuple val(sample_name), path("resistance_prediction_report.json"), emit: json
    tuple val(sample_name), path("final.vcf"), emit: vcf
    tuple val(sample_name), path("variants.parquet"), emit: variants_table, optional: true
    tuple val(sample_name), path("mutations.parquet"), emit: mutations_table, optional: true
    tuple val(sample_name), path("effects.parquet"), emit: effects_table, optional: true
    tuple val(sample_name), path("predictions.parquet"), emit: predictions_table, optional: true
    

    script:
    MIN_DP = seq_platform == 'illumina' ? 3 : 5
    """
    variants=variants.vcf
    all_calls=all_calls.vcf

    if gzip -t variants.vcf; then
        gzip -dc variants.vcf > uncompressed_variants.vcf
        variants=uncompressed_variants.vcf
    fi

    if gzip -t all_calls.vcf; then
        gzip -dc all_calls.vcf > uncompressed_all_calls.vcf
        all_calls=uncompressed_all_calls.vcf
    fi


    # Merge in GVCF rows at resistance SNPs to ensure we can detect null calls
    # at these sites (sometimes minos doesn't give us these)
    merge-vcfs --minos_vcf \${variants} \
            --gvcf \${all_calls} \
            --resistant-positions ${null_positions} \
            --output "${sample_name}.vcf" \
            --min_dp ${MIN_DP}

    gnomonicus --genome_object ${reference} \
            --catalogue ${catalogue} \
            --vcf_file "${sample_name}.vcf" \
            --json \
            --output_dir . \
            --min_dp ${MIN_DP} \
            --csvs all \
            --parquet \
            --json_resistance_genes_only


    mv "${sample_name}.gnomonicus-out.json" resistance_prediction_report.json
    mv "${sample_name}.vcf" final.vcf
    mv "${sample_name}.variants.parquet" variants.parquet || true
    mv "${sample_name}.mutations.parquet" mutations.parquet || true
    mv "${sample_name}.effects.parquet" effects.parquet || true
    mv "${sample_name}.predictions.parquet" predictions.parquet || true

    find . -type f -name "uncompressed*.vcf" -delete
    """

    stub:
    """
    touch resistance_prediction_report.json
    """
}
