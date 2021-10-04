#!/usr/bin/env nextflow

/*
========================================================================================
    DEFAULT PARAMETERS
========================================================================================
*/

nextflow.enable.dsl=2
params.synid = "syn26240435.txt"
params.outdir = "./"
params.results_dir = "samplesheet"
params.strandedness = "reverse"

/*
========================================================================================
    DEFINE PROCESSES
========================================================================================
*/

// Generate list of Synapse IDs and corresponding files to be downloaded
process SYNAPSE_LIST {
    publishDir "${params.outdir}/synapse/"

    input:
    val synid

    output:
    path "*.synlist.csv", emit: synlist_csv

    script:
    """
    synapse list -l $synid | cut -c-11 > ${synid}.synlist.csv
    """
}

// Download FASTQ files by Synapse ID
process SYNAPSE_GET {
    publishDir "${params.outdir}/fastq/"

    input:
    val(synapseID)

    output:
    path "*", emit: fastq

    script:
    """
    synapse get $synapseID
    """
}

// Download metadata by Synapse ID
process SYNAPSE_SHOW {
    publishDir "${params.outdir}/synapse/"

    input:
    val(synapseID)

    output:
    path "*.metadata.txt", emit: metadata

    script:
    """
    synapse show $synapseID | sed -n '1,3p;15,16p;20p;23p' > ${synapseID}.metadata.txt
    """
}

process READ_PAIRS_TO_SAMPLESHEET {
    publishDir "${params.outdir}/samplesheet/"

    input:
    tuple val(id), val(files)
    val strandedness

    output:
    path("*samplesheet.csv"), emit: samplesheet

    exec:
    // Add relevant fields to the beginning of the map
    pipeline_map = [
        sample  : "${id}",
        fastq_1 : "${params.outdir}/${params.results_dir}/${files[0].getBaseName()}",
        fastq_2 : "${params.outdir}/${params.results_dir}/${files[1].getBaseName()}"
    ]
    // Add Strandedness
    pipeline_map << [ strandedness: "${strandedness}" ]
    
    // Create Samplesheet
    samplesheet  = pipeline_map.keySet().collect{ '"' + it + '"'}.join(",") + '\n'
    samplesheet += pipeline_map.values().collect{ '"' + it + '"'}.join(",")

    def samplesheet_file2 = task.workDir.resolve("${pipeline_map.sample}.samplesheet.csv")
    samplesheet_file2.text = samplesheet
}

process MERGE_SAMPLESHEET {
    publishDir "${params.outdir}"

    input:
    path ('samplesheets/*')

    output:
    path "samplesheet.csv", emit: samplesheet

    script:
    """
    head -n 1 `ls ./samplesheets/* | head -n 1` > samplesheet.csv
    for fileid in `ls ./samplesheets/*`; do
        awk 'NR>1' \$fileid >> samplesheet.csv
    done
    """
}

/*
========================================================================================
    VALIDATE INPUT SYNAPSE IDS
========================================================================================
*/

// Create channel for input Synapse directory IDs
Channel
    .from(file(params.synid, checkIfExists: true))
    .splitCsv(header:false, sep:'', strip:true)
    .map { it[0] }
    .unique()
    .set { ch_ids }

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    SYNAPSE_LIST (
        ch_ids,
    )

    // Create channel for FQ SynapseIDs
    SYNAPSE_LIST
        .out
        .synlist_csv
        .splitCsv(header:false, strip:true).flatten()
        .set { ch_samples }

    // Download FQ Files by SynapseID
    SYNAPSE_GET(
        ch_samples,
    )

    // Create Read Pairs Channel
    SYNAPSE_GET
        .out
        .collect().flatten()
        .toSortedList().flatten()
        .map { meta ->  
            def sampleId = meta.name.toString().tokenize('_').get(0)
            [sampleId, meta]
        }
        .groupTuple()
        .set{ ch_read_pairs }

    // Download FQ Metadata by SynapseID
    SYNAPSE_SHOW(
        ch_samples,
    )

    // Clean Metadata
    SYNAPSE_SHOW
        .out
        .metadata
        .splitCsv()

    // Create Samplesheet
    READ_PAIRS_TO_SAMPLESHEET(
        ch_read_pairs,
        params.strandedness
    )

    // Merge Samplesheets
    MERGE_SAMPLESHEET (
        READ_PAIRS_TO_SAMPLESHEET.out.samplesheet.collect()
    )
}
