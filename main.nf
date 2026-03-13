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
            --sample                  Path to the sample minos VCF
            --species                 Species of Mycobacteria. Default is 'Mycobacterium tuberculosis'
            --genbank_reference_dir   Path to the directory containing genbank reference files
            --gvcf                    Path to the non-compressed gvcf file
            --catalogue               Path to the resistance catalogue
            --null_positions          Path to the null positions file
            --seq_platform            Sequencing platform used ('illumina' or 'ont'). Default is 'illumina'
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
        --sample                  ${params.sample}
        --species                 ${params.species}
        --genbank_reference_dir   ${params.genbank_reference_dir}
        --gvcf                    ${params.gvcf}
        --catalogue               ${params.catalogue}
        --null_positions          ${params.null_positions}
        --seq_platform            ${params.seq_platform}

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
    input = sample.merge(gvcf).map { sample_id, vcf, gvcf -> tuple(sample_id, vcf, params.species, gvcf) }

    genbank_reference_dir = Channel.fromPath(params.genbank_reference_dir, checkIfExists: true).first()
    catalogue = Channel.fromPath(params.catalogue, checkIfExists: true).first()
    null_positions = Channel.fromPath(params.null_positions, checkIfExists: true).first()

    gnomonicus_workflow(input, params.seq_platform, genbank_reference_dir, catalogue, null_positions)
}

workflow gnomonicus_workflow {
    take:
    samples // Channel of tuples with sample_id, vcf, species, gvcf
    seq_platform
    reference_dir
    catalogue
    null_positions

    main:
    reference_picked = pick_reference(samples, reference_dir)
    gnomonicus_out = runPrediction(samples, seq_platform, reference_picked.reference, catalogue, null_positions)

    emit:
    gnomonicus_json = gnomonicus_out.json
    gnomonicus_vcf = gnomonicus_out.vcf
    gnomonicus_variants_csv = gnomonicus_out.variants_csv
    gnomonicus_mutations_csv = gnomonicus_out.mutations_csv
    gnomonicus_effects_csv = gnomonicus_out.effects_csv
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
            Expects each sample to have a corresponding directory with vcf and gvcf.
            Requires all samples to be of the same species.

            Mandatory parameters:
            ------------------------------------------------------------------------
            --samples                 Path pattern for sample VCFs (e.g. 'samples/*/*.vcf')
            --gvcfs                   Path pattern for gvcf files (e.g. 'samples/*/*.gvcf')
            --species                 Species of Mycobacteria. Default is 'Mycobacterium tuberculosis'
            --genbank_reference_dir   Path to the directory containing genbank reference files
            --catalogue               Path to the resistance catalogue
            --null_positions          Path to the null positions file
            --seq_platform            Sequencing platform used ('illumina' or 'ont'). Default is 'illumina'
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
        --samples                 ${params.samples}
        --gvcfs                   ${params.gvcfs}
        --genbank_reference_dir   ${params.genbank_reference_dir}
        --catalogue               ${params.catalogue}
        --null_positions          ${params.null_positions}
        --seq_platform            ${params.seq_platform}

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

    input = samples.join(gvcfs).map { sample_id, vcf, gvcf -> tuple(sample_id, vcf, params.species, gvcf) }

    input.take(3).view()

    genbank_reference_dir = Channel.fromPath(params.genbank_reference_dir, checkIfExists: true).first()
    catalogue = Channel.fromPath(params.catalogue, checkIfExists: true).first()
    null_positions = Channel.fromPath(params.null_positions, checkIfExists: true).first()

    reference_picked = pick_reference(samples, reference_dir)
    runPrediction(input, params.seq_platform, reference_picked.reference, catalogue, null_positions)
}

process pick_reference {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    // Using viridian-utils as it has a small footprint (~20MB) and includes bash
    container params.container_prefix + '/vtap/viridian-utils:1.0.2'
    cpus 1
    maxRetries 5
    memory "1GB"

    pod label: "name", value: "tb-predict-pipeline:pick_reference"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path("variants.vcf"), val(species), path("all_rows.gvcf")
    path genbank_reference_dir

    output:
    tuple val(sample_name), path("reference.gbk"), emit: reference, optional: true

    script:
    """
    if [[ "$species" == "Mycobacterium tuberculosis" ]]; then
        cp "${genbank_reference_dir}/NC_000962.3.gbk" reference.gbk
    elif [[ "$species" == "Mycobacterium abscessus" ]]; then
        cp "${genbank_reference_dir}/CU458896.1.gbk" reference.gbk
    elif [[ "$species" == "Mycobacterium avium" ]]; then
        cp "${genbank_reference_dir}/CP018019.2.gbk" reference.gbk
    elif [[ "$species" == "Mycobacterium chelonae" ]]; then
        cp "${genbank_reference_dir}/CP031516.1.gbk" reference.gbk
    elif [[ "$species" == "Mycobacterium chelonae_A" ]]; then
        cp "${genbank_reference_dir}/CP031516.1.gbk" reference.gbk
    elif [[ "$species" == "Mycobacterium gwanakae" ]]; then
        cp "${genbank_reference_dir}/CP031516.1.gbk" reference.gbk
    elif [[ "$species" == "Mycobacterium chimaera" ]]; then
        cp "${genbank_reference_dir}/CP015278.1.gbk" reference.gbk
    elif [[ "$species" == "Mycobacterium fortuitum" ]]; then
        cp "${genbank_reference_dir}/AP025518.1.gbk" reference.gbk
    elif [[ "$species" == "Mycobacterium intracellulare" ]]; then
        cp "${genbank_reference_dir}/NZ_CP085945.1.gbk" reference.gbk
    elif [[ "$species" == "Mycobacterium kansasii" ]]; then
        cp "${genbank_reference_dir}/CP006835.1.gbk" reference.gbk

    # These are mostly for testing (with the species TEST not being real a genome)
    elif [[ "$species" == "TEST" ]]; then
        cp "${genbank_reference_dir}/TEST-DNA.gbk" reference.gbk
    elif [[ "$species" == "SARS-CoV2" ]]; then
        cp "${genbank_reference_dir}/NC_045512.2.gbk" reference.gbk
    else
        echo "Unsupported species: $species" >&2
    fi
    """

    stub:
    """
    touch reference.gbk
    """
}


//Run gnomonicus
process runPrediction {
    publishDir "${params.publish_dir}", enabled: params.publish_dir != "", mode: "copy", saveAs: { filename -> sample_name + "_" + filename }
    container params.container_prefix + "/oxfordmmm/gnomonicus:v3.1.1"
    cpus 2
    maxRetries 5
    memory {
        params.testing == "" ? 8.GB * (0.8 + (task.attempt / 5)) : "6GB"
    }

    pod label: "name", value: "tb-predict-pipeline:runPrediction"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), path("variants.vcf"), val(species), path("all_rows.gvcf")
    val seq_platform
    tuple val(sample_name), path(reference)
    path catalogue
    path null_positions

    output:
    tuple val(sample_name), path("resistance_prediction_report.json"), emit: json
    tuple val(sample_name), path("merged.vcf"), emit: vcf
    tuple val(sample_name), path("variants.csv"), emit: variants_csv, optional: true
    tuple val(sample_name), path("mutations.csv"), emit: mutations_csv, optional: true
    tuple val(sample_name), path("effects.csv"), emit: effects_csv, optional: true

    script:
    MIN_DP = seq_platform == 'illumina' ? 3 : 5
    // For now we only have a TB catalogue (and catalogues for testing), 
    // so only pass the catalogue arg if the species is one of these
    // In future this will likely need updating to dynamically select catalogues for other species
    CATALOGUE = species == "Mycobacterium tuberculosis" || species == "TEST" || species == "SARS-CoV2" ? "--catalogue " + catalogue : ""
    MERGE_VCFS = CATALOGUE == "" ? "cp variants.vcf ${sample_name}.vcf" : "merge-vcfs --minos_vcf variants.vcf --gvcf all_rows.gvcf --resistant-positions ${null_positions} --output \"${sample_name}.vcf\" --min_dp ${MIN_DP}"
    """
    ${MERGE_VCFS}
    gnomonicus --genome_object ${reference} $CATALOGUE --vcf_file "${sample_name}.vcf" --json --csvs all --output_dir . --min_dp ${MIN_DP}

    mv "${sample_name}.gnomonicus-out.json" resistance_prediction_report.json
    mv "${sample_name}.vcf" merged.vcf # This will be renamed "final.vcf" in a future update

    # Depending on the sample, these may or may not exist so suppress errors
    mv "${sample_name}.variants.csv" variants.csv 2> /dev/null || true
    mv "${sample_name}.mutations.csv" mutations.csv 2> /dev/null || true
    mv "${sample_name}.effects.csv" effects.csv 2> /dev/null || true
    """

    stub:
    """
    touch resistance_prediction_report.json
    touch merged.vcf
    """
}
