#!/usr/bin/env nextflow

// default parameters
params.seq_platform = "illumina"

workflow {
    def samples_is_glob = params.samples.contains('*')
    def gvcfs_is_glob = params.gvcfs.contains('*')

    if (samples_is_glob != gvcfs_is_glob) {
        error("params.samples and params.gvcfs must either both be single files or both contain '*'")
    }

    // replace is to catch both .vcf and .vcf.gz files
    samples = channel.fromPath("${params.samples}", checkIfExists: true, glob: true)
        .ifEmpty { error("cannot find any reads matching ${params.samples}") }
        .map { it -> tuple(it.baseName.replace(".vcf", "").replace(".gvcf", ""), it) }

    gvcfs_ch = channel.fromPath("${params.gvcfs}", checkIfExists: true, glob: true)
        .ifEmpty { error("cannot find any reads matching ${params.gvcfs}") }
        .map { it -> tuple(it.baseName.replace(".vcf", "").replace(".gvcf", ""), it) }


    if (samples_is_glob) {
        // Multiple samples: match sample and gVCF by ID
        input = samples
            .join(gvcfs_ch)
            .map { sample_id, vcf, gvcf ->
                tuple(sample_id, "ref_id", vcf, gvcf)
            }
    }
    else {
        // Single sample + single gVCF: don't match IDs
        input = samples
            .combine(gvcfs_ch)
            .map { sample_id, vcf, _gvcf_id, gvcf ->
                tuple(sample_id, "ref_id", vcf, gvcf)
            }
    }



    reference = channel.fromPath(params.reference, checkIfExists: true).first()
    catalogue = channel.fromPath(params.catalogue, checkIfExists: true).first()
    null_positions = channel.fromPath(params.null_positions, checkIfExists: true).first()

    input = input.combine(reference)
        .combine(catalogue)
        .combine(null_positions)

    input.take(3).view()

    gnomonicus_workflow(
        input,
        params.seq_platform,
    )
}

workflow gnomonicus_workflow {
    take:
    samples // Channel of tuples with sample_id, ref_id, vcf, gvcf, genbank, catalogue, null positions
    seq_platform

    main:
    gnomonicus_out = runPrediction(
        samples,
        seq_platform,
    )

    emit:
    gnomonicus_json = gnomonicus_out.json
    gnomonicus_vcf = gnomonicus_out.vcf
    gnomonicus_variants = gnomonicus_out.variants_table
    gnomonicus_mutations = gnomonicus_out.mutations_table
    gnomonicus_effects = gnomonicus_out.effects_table
    gnomonicus_predictions = gnomonicus_out.predictions_table
}


//Run gnomonicus
process runPrediction {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container params.container_prefix + "/oxfordmmm/gnomonicus:3.1.6"
    cpus 2
    maxRetries 5
    memory { params.testing == "" ? 8.GB + (4.GB * (task.attempt - 1)) : "6GB" }

    pod label: "name", value: "tb-predict-pipeline:runPrediction"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), val(ref_id), path("variants.vcf"), path("all_calls.vcf"), path(reference), path(catalogue), path(null_positions)
    val seq_platform

    output:
    tuple val(sample_name), val(ref_id), path("resistance_prediction_report.json"), emit: json
    tuple val(sample_name), val(ref_id), path("final.vcf"), emit: vcf
    tuple val(sample_name), val(ref_id), path("variants.parquet"), emit: variants_table, optional: true
    tuple val(sample_name), val(ref_id), path("mutations.parquet"), emit: mutations_table, optional: true
    tuple val(sample_name), val(ref_id), path("effects.parquet"), emit: effects_table, optional: true
    tuple val(sample_name), val(ref_id), path("predictions.parquet"), emit: predictions_table, optional: true

    script:
    MIN_DP = seq_platform == 'illumina' ? 3 : 5
    CATALOGUE_CMD = catalogue.name == "EMPTY_catalogue.csv" ? "" : "--catalogue " + catalogue
    // Only merge VCFs if a catalogue is provided, otherwise just copy the variants.vcf to the output
    MERGE_VCFS_CMD = (CATALOGUE_CMD == ""
        ? "cp variants.vcf ${sample_name}.vcf"
        : "merge-vcfs --minos_vcf variants.vcf --gvcf all_calls.vcf --resistant-positions ${null_positions} --output \"${sample_name}.vcf\" --min_dp ${MIN_DP}")
    """
    ${MERGE_VCFS_CMD}
    gnomonicus --genome_object ${reference} \
        ${CATALOGUE_CMD} \
        --vcf_file "${sample_name}.vcf" \
        --json \
        --output_dir . \
        --min_dp ${MIN_DP} \
        --csvs all \
        --parquet \
        --json_resistance_genes_only

    mv "${sample_name}.gnomonicus-out.json" resistance_prediction_report.json
    mv "${sample_name}.vcf" final.vcf

    # Depending on the sample, these may or may not exist so suppress errors
    mv "${sample_name}.variants.parquet" variants.parquet || true
    mv "${sample_name}.mutations.parquet" mutations.parquet || true
    mv "${sample_name}.effects.parquet" effects.parquet || true
    mv "${sample_name}.predictions.parquet" predictions.parquet || true
    """

    stub:
    """
    touch resistance_prediction_report.json
    touch final.vcf
    """
}
