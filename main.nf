#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

def help() {
	log.info"""
	    Run the analysis:
		nextflow <nextflow options> run main.nf -profile <profiles> --input <fastqs> --outdir <dir> --ticket <batch> <options>

		Initialize the Shiver input:
		nextflow <nextflow options> run main.nf --init -profile <profiles> --outdir <dir> <options>

		Create database for host genome:
		nextflow <nextflow options> run main.nf --setup -profile <profiles> --outdir <dir> <options>

		Nextflow options:
		-c                    Path to additional config file [Default: Nextflow uses nextflow.config in current and script
		                      directory in addition to configurations in ~/.nextflow/config]
		-C                    Path to config file. Other config files are ignored.

		Options:
		-profile            Comma-separated list of run profiles to use: local,slurm,docker,singularity
		--input             Input directory with fastq files
		--outdir            Directory for results
		--ticket            Batch name
		--primers           Primers used during amplification [Default:?]
		--adapters          Sequencing adapters [Default:?]
		--initDir           Shiver initialization directory for configurations
		--config            Shiver configuration file [Default:?]
		--references        HIV reference dataset
        --hostGenome        Directory with reference database for host genome
        --hostGenomeBase    Name of reference database for host genome
        --hostURL           URL to download host reference genome from
		--min_cov_eti       Coverage threshold in time of infection calculations [Default: 300]
		--large             Use more computing resources
		--help              Display this help message and exit

        Initialization options:
        --primers           Primers used during amplification [Default:?]
		--adapters          Sequencing adapters [Default:?]
		--config            Shiver configuration file [Default:?]
		--references        HIV reference dataset

        Setup options:
        --hostGenome        Directory to create reference database for host genome in
        --hostGenomeBase    Name of reference database for host genome
        --hostFasta         Host reference genome path
		"""
}

if (params.help) {
	help()
		exit 0
}

// Verify input parameters

if (!params.outdir){
    println("Please specify directory for output with --outdir")
    System.exit(1)
}

if ( params.setup ) {
    if (!params.hostGenome){
        println("Please specify reference database directory with --hostGenome")
        System.exit(1)
    }
    if (!params.hostGenomeBase){
        println("Please specify reference database name with --hostGenomeBase")
        System.exit(1)
    }
    if (!params.hostFasta){
        println("Please specify url for host reference fasta with --hostFasta")
        System.exit(1)
    }
}
else if ( params.init ) {
    if (!params.references){
        println("Please specify reference fasta file with --references")
        System.exit(1)
    }
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
}

// include workflows and modules
include {timeAnalysis} from './workflows/time.nf'
include {initialisation} from './workflows/initialisation.nf'
include {setup} from './workflows/pipeline_setup.nf'
include {getVersion} from './modules/version.nf'

workflow {
    getVersion()

    if ( params.setup ) {
        // Set initialisation input
        Channel.fromPath(params.hostFasta)
            .set{ ch_hostFasta }

        // Run setup
        main:
            setup(ch_hostFasta)
    }

    else if ( params.init ) {
        // Set initialisation input
        Channel.fromPath(params.references)
            .set{ ch_referenceFile }
        Channel.fromPath(params.primers)
            .set{ ch_primersFile }
        Channel.fromPath(params.adapters)
            .set{ ch_adaptersFile }
        Channel.fromPath(params.config)
            .set{ ch_shiverConfigFile }

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
        Channel.fromPath(params.primers)
            .set{ ch_primersFile }
        Channel.fromPath(params.adapters)
            .set{ ch_adaptersFile }
        Channel.fromPath(params.config)
            .set{ ch_shiverConfigFile }

        // run analysis
        main:
            timeAnalysis(ch_fastq, ch_primersFile, ch_adaptersFile, ch_shiverConfigFile, ch_shiverInitDir, ch_hostGenome)
    }

}

