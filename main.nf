#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

def help() {
	log.info"""
		nextflow <nextflow options> run main.nf <options>

		Required options:
		-profile        Comma-separated list of run profiles to use: local,slurm,docker,singularity
		--input         Input directory with fastq files
		--outdir        Directory for results
		--ticket        Batch name

		Additional options:
		--primers       Primers used during amplification [Default:?]
		--adapters      Sequencing adapters [Default:?]
		--initDir       Shiver initialization directory for configurations
		--config        Shiver configuration file [Default:?]
		--min_cov_eti   Coverage threshold in time of infection calculations [Default: 300]
		--large         Use more computing resources
		--help          Display this help message and exit

		Initialisation of Shiver directory:
		nextflow <nextflow options> run main.nf --init <options>

        Required options:
        -profile        Comma-separated list of run profiles to use: local,slurm,docker,singularity
        --primers       Primers used during amplification [Default:?]
		--adapters      Sequencing adapters [Default:?]
		--config        Shiver configuration file [Default:?]
		--references    HIV reference dataset

		Nextflow options:
		Specify config file with Nextflow options:
		-c        Path to additional config file [Default: Nextflow uses nextflow.config in current and script
		          directory in addition to configurations in ~/.nextflow/config]
		-C        Path to config file. Other config files are ignored.
		"""
}

if (params.help) {
	help()
		exit 0
}

// Verify input parameters
if (!params.primers){
    println("Please specify primers with --primers")
    System.exit(1)
}
if (!params.adapters){
    println("Please specify adapters with --adapters")
    System.exit(1)
}
if (!params.config){
    println("Please specify Shiver config file with --config")
    System.exit(1)
}
if (!params.outdir){
    println("Please specify directory for results with --outdir")
    System.exit(1)
}
if ( params.init ) {
    if (!params.references){
        println("Please specify reference fasta file with --references")
        System.exit(1)
    }
}
else {
    if (!params.input){
        println("Please specify input directory with --input")
        System.exit(1)
    }
    if (!params.ticket){
        println("Please specify a name for the batch --ticket")
        System.exit(1)
    }
    if (!params.initDir){
        println("Please specify Shiver intialisation directory with --initDir")
        System.exit(1)
    }
}

// include workflows
include {timeAnalysis} from './workflows/time.nf'
include {initialisation} from './workflows/initialisation.nf'

workflow {

    // Set general workflow input
    Channel.fromPath(params.primers)
        .set{ ch_primersFile }
    Channel.fromPath(params.adapters)
        .set{ ch_adaptersFile }
    Channel.fromPath(params.config)
        .set{ ch_shiverConfigFile }

    if ( params.init ) {
        // Set initialisation input
        Channel.fromPath(params.references)
            .set{ ch_referenceFile }
        // Run initialisation
        main:
            initialisation(ch_primersFile, ch_adaptersFile, ch_shiverConfigFile, ch_referenceFile)
    }
    else {
        // Set analysis input
        Channel.fromFilePairs( params.input, size: 2 )
            .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
            .ifEmpty { exit 1, "Cannot find fastq-files in input directory" }
            .set{ch_fastq}
        Channel.fromPath(params.initDir)
            .set{ ch_shiverInitDir }
        Channel.fromPath(params.hostGenome)
            .set{ ch_hostGenome }

        // run analysis
        main:
            timeAnalysis(ch_fastq, ch_primersFile, ch_adaptersFile, ch_shiverConfigFile, ch_shiverInitDir, ch_hostGenome)
    }

}

