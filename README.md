# nf-fetchngs-synapse

This is a Nextflow tool that uses the SynapseClient CLI to download FASTQ files and generate a samplesheet for use in the nf-core rnaseq pipeline.

### Inputs
* `synID` - the SynapseID of the directory containing all the required FASTQ files
* `.synapseConfig` - Synapse configuration file required to log in to Synapse. For more information, see [Synapse Docs - Client Configuration](https://help.synapse.org/docs/Client-Configuration.1985446156.html)
* (Optional) `outdir` - Directory to which all outputs are saved. Defaults to current working directory.
* (Optional) `results_dir` - Subdirectory to which the final samplesheet is saved. Defaults to `./samplesheet`
* (Optional) `strandedness` - Strandedness of FASTQ files. Defaults to `reverse`.

