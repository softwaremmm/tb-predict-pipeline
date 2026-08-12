#!/usr/bin/env nextflow

// default parameters
params.seq_platform = "illumina"

workflow workflow {
    ANSI_GREEN = "\033[1;32m"
    ANSI_RESET = "\033[0m"

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
            --reference_data_dir      Path to the directory containing reference files (genbank, catalogue, null positions)
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
        --reference_data_dir      ${params.reference_data_dir}
        --seq_platform            ${params.seq_platform}

        Runtime data:
        ------------------------------------------------------------------------
        Running with profile  ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
        Running as user       ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
        Launch directory      ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
        """.stripIndent()
    )

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
                tuple(sample_id, params.species, vcf, gvcf)
            }
    }
    else {
        // Single sample + single gVCF: don't match IDs
        input = samples
            .combine(gvcfs_ch)
            .map { sample_id, vcf, _gvcf_id, gvcf ->
                tuple(sample_id, params.species, vcf, gvcf)
            }
    }


    input.take(3).view()

    reference_data_dir = channel.fromPath(params.reference_data_dir, checkIfExists: true).first()

    gnomonicus_workflow(
        input,
        reference_data_dir,
        params.seq_platform,
    )
}

workflow gnomonicus_workflow {
    take:
    samples // Channel of tuples with sample_id, species, vcf, gvcf
    reference_data_dir // Path to the directory containing reference files (genbank, catalogue, null positions)
    seq_platform

    main:
    pick_reference(samples, reference_data_dir)
    gnomonicus_out = runPrediction(
        samples.join(
            pick_reference,
            by: [0, 1]
        ),
        seq_platform,
    )

    emit:
    gnomonicus_json = gnomonicus_out.json
    gnomonicus_vcf = gnomonicus_out.vcf
    gnomonicus_variants_csv = gnomonicus_out.variants_csv
    gnomonicus_mutations_csv = gnomonicus_out.mutations_csv
    gnomonicus_effects_csv = gnomonicus_out.effects_csv
    gnomonicus_predictions_csv = gnomonicus_out.predictions_csv
}


// expects ref data directory to contain subfolder for each species
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
    tuple val(sample_name), val(species), path("variants.vcf"), path("all_rows.gvcf")
    path reference_data_dir

    output:
    tuple val(sample_name), val(species), path("*reference.gbk"), path("*catalogue.csv"), path("*null_positions.txt"), emit: reference_data

    script:
    """
    if [[ "${species}" == "Mycobacterium tuberculosis" ]]; then
        cp "${reference_data_dir}/NC_000962.3.gbk" reference.gbk
        cp "${reference_data_dir}/tb_catalogue.csv" catalogue.csv
        cp "${reference_data_dir}/null_positions.txt" null_positions.txt
    elif [[ "${species}" == "Mycobacterium abscessus" ]]; then
        cp "${reference_data_dir}/CU458896.1.gbk" reference.gbk
        touch EMPTY_catalogue.csv
        touch EMPTY_null_positions.txt
    elif [[ "${species}" == "Mycobacterium avium" ]]; then
        cp "${reference_data_dir}/CP018019.2.gbk" reference.gbk
        touch EMPTY_catalogue.csv
        touch EMPTY_null_positions.txt
    elif [[ "${species}" == "Mycobacterium chelonae" ]]; then
        cp "${reference_data_dir}/CP031516.1.gbk" reference.gbk
        touch EMPTY_catalogue.csv
        touch EMPTY_null_positions.txt
    elif [[ "${species}" == "Mycobacterium chelonae_A" ]]; then
        cp "${reference_data_dir}/CP031516.1.gbk" reference.gbk
        touch EMPTY_catalogue.csv
        touch EMPTY_null_positions.txt
    elif [[ "${species}" == "Mycobacterium gwanakae" ]]; then
        cp "${reference_data_dir}/CP031516.1.gbk" reference.gbk
        touch EMPTY_catalogue.csv
        touch EMPTY_null_positions.txt
    elif [[ "${species}" == "Mycobacterium chimaera" ]]; then
        cp "${reference_data_dir}/CP015278.1.gbk" reference.gbk
        touch EMPTY_catalogue.csv
        touch EMPTY_null_positions.txt
    elif [[ "${species}" == "Mycobacterium fortuitum" ]]; then
        cp "${reference_data_dir}/AP025518.1.gbk" reference.gbk
        touch EMPTY_catalogue.csv
        touch EMPTY_null_positions.txt
    elif [[ "${species}" == "Mycobacterium intracellulare" ]]; then
        cp "${reference_data_dir}/NZ_CP085945.1.gbk" reference.gbk
        touch EMPTY_catalogue.csv
        touch EMPTY_null_positions.txt
    elif [[ "${species}" == "Mycobacterium kansasii" ]]; then
        cp "${reference_data_dir}/CP006835.1.gbk" reference.gbk
        touch EMPTY_catalogue.csv
        touch EMPTY_null_positions.txt

    # These are mostly for testing (with the species TEST not being real a genome)
    elif [[ "${species}" == "TEST" ]]; then
        cp "${reference_data_dir}/TEST-DNA.gbk" reference.gbk
        cp "${reference_data_dir}/TEST-DNA-catalogue.csv" catalogue.csv
        touch EMPTY_null_positions.txt
    elif [[ "${species}" == "SARS-CoV2" ]]; then
        cp "${reference_data_dir}/NC_045512.2.gbk" reference.gbk
        cp "${reference_data_dir}/NC_045512.2-test-catalogue.csv" catalogue.csv
        touch EMPTY_null_positions.txt
    else
        touch EMPTY_reference.gbk
        touch EMPTY_catalogue.csv
        touch EMPTY_null_positions.txt
        echo "Unsupported species: ${species}" >&2
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
    memory { params.testing == "" ? 8.GB + (4.GB * (task.attempt - 1)) : "6GB" }

    pod label: "name", value: "tb-predict-pipeline:runPrediction"
    pod label: "sample_id", value: "${params.sample_id}"
    pod label: "run_id", value: "${params.run_id}"

    input:
    tuple val(sample_name), val(species), path("variants.vcf"), path("all_calls.vcf"), path(reference), path(catalogue), path(null_positions)
    val seq_platform

    output:
    tuple val(sample_name), path("resistance_prediction_report.json"), val(species), emit: json
    tuple val(sample_name), path("final.vcf"), val(species), emit: vcf
    tuple val(sample_name), path("variants.csv"), val(species), emit: variants_csv, optional: true
    tuple val(sample_name), path("mutations.csv"), val(species), emit: mutations_csv, optional: true
    tuple val(sample_name), path("effects.csv"), val(species), emit: effects_csv, optional: true
    tuple val(sample_name), path("predictions.csv"), val(species), emit: predictions_csv, optional: true

    script:
    MIN_DP = seq_platform == 'illumina' ? 3 : 5
    // For now we only have a TB catalogue (and catalogues for testing),
    // so only pass the catalogue arg if the species is one of these
    // In future this will likely need updating to dynamically select catalogues for other species
    CATALOGUE_CMD = catalogue == "EMPTY_catalogue.csv" ? "" : "--catalogue " + catalogue
    MERGE_VCFS_CMD = (CATALOGUE_CMD == ""
        ? "cp variants.vcf ${sample_name}.vcf"
        : "merge-vcfs --minos_vcf variants.vcf --gvcf all_rows.gvcf --resistant-positions ${null_positions} --output \"${sample_name}.vcf\" --min_dp ${MIN_DP}")
    """
    ${MERGE_VCFS_CMD}
    gnomonicus --genome_object ${reference} ${CATALOGUE_CMD} --vcf_file "${sample_name}.vcf" --json --csvs all --output_dir . --min_dp ${MIN_DP}

    mv "${sample_name}.gnomonicus-out.json" resistance_prediction_report.json
    mv "${sample_name}.vcf" final.vcf

    # Depending on the sample, these may or may not exist so suppress errors
    mv "${sample_name}.variants.csv" variants.csv 2> /dev/null || true
    mv "${sample_name}.mutations.csv" mutations.csv 2> /dev/null || true
    mv "${sample_name}.effects.csv" effects.csv 2> /dev/null || true
    mv "${sample_name}.predictions.csv" predictions.csv 2> /dev/null || true
    """

    stub:
    """
    touch resistance_prediction_report.json
    touch final.vcf
    """
}
