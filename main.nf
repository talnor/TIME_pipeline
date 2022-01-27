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
		--config        Shiver configuration file
		--min_cov_eti   Coverage threshold in time of infection calculations [Default: 300]
		--help          Display this help message and exit

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
if (!params.input){
    println("Please specify input directory with --input")
    System.exit(1)
}
if (!params.outdir){
    println("Please specify directory for results with --outdir")
    System.exit(1)
}
if (!params.ticket){
    println("Please specify a name for the batch --ticket")
    System.exit(1)
}

// include workflows
include {timeAnalysis} from './workflows/time.nf'

workflow {
    // Find fastq files
    Channel.fromFilePairs( params.input, size: 2 )
    .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
    .ifEmpty { exit 1, "Cannot find fastq-files in input directory" }
    .set{ch_fastq}

    Channel.fromPath(params.primers)
        .set{ ch_primersFile }
    Channel.fromPath(params.adapters)
        .set{ ch_adaptersFile }
    Channel.fromPath(params.config)
        .set{ ch_shiverConfigFile }
    Channel.fromPath(params.initDir)
        .set{ ch_shiverInitDir }
    Channel.fromPath(params.hostGenome)
        .set{ ch_hostGenome }

    // run analysis
    main:
    timeAnalysis(ch_fastq, ch_primersFile, ch_adaptersFile, ch_shiverConfigFile, ch_shiverInitDir, ch_hostGenome)
}

