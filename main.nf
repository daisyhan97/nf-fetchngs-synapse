#!/usr/bin/env nextflow

/*
========================================================================================
    DEFAULT PARAMETERS
========================================================================================
*/

params.synid = "syn26240435"
params.synapse_config = "~/.synapseConfig"
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
    file synapseconfig

    output:
    path "*.synlist.csv", emit synlist_csv

    script:
    """
    synapse -c $synapseconfig list $synid > ${synid}.synlist.csv
    """
}

// Download FASTQ files by Synapse ID
process DOWNLOAD_SYNID {
    publishDir "${params.outdir}/fastq/"

    input:
    tuple (synapseID, fileName)
    file synapseconfig

    output:
    path "*", emit fastq

    script:
    """
    synapse -c $synapseconfig get $synapseID
    """
}

// Generate Samplesheet from Read Pairs
process READ_PAIRS_TO_SAMPLESHEET {

    input:
    set id, files
    val strandedness

    exec:
    // Split read channel to corresponding files
    (id, fastq_1, fastq_2) = [id, files[0].getBaseName(), files[1].getBaseName()]

    // Generate samplesheet
    File file = new File("${params.outdir}/${params.results_dir}/samplesheet.txt")
    file.write "Sample, Read_1, Read_2, Standedness\n"
    file.append "$id, ${params.outdir}/${params.results_dir}/${fastq_1}, ${params.outdir}/${params.results_dir}/${fastq_2}, $strandedness\n"
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

    // Download CSV of Synapse IDs and FQ File Names
    SYNAPSE_LIST (
        ch_ids,
        ch_synapse_config
    )

    // Create channel for FQ SynapseIDs
    SYNAPSE_LIST
        .out
        .synlist_csv
        .splitCsv(header:false)
        .map{ row-> tuple(row.synapseId, row.fileName) }
        .set { ch_samples }

    // Download FQ Files by SynapseID
    DOWNLOAD_SYNID(
        ch_samples,
        ch_synapse_config
    )

    // Create Read Pairs Channel
    DOWNLOAD_SYNID
        .out
        .fastq
        .collect().flatten()
        .toSortedList().flatten()
        .map { meta ->  
            def sampleId = meta.name.toString().tokenize('_').get(0)
            [sampleId, meta]
        }
        .groupTuple()
        .set{ ch_read_pairs }

    // Create Samplesheet
    READ_PAIRS_TO_SAMPLESHEET(
        ch_read_pairs,
        params.strandedness
    )
}