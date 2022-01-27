# TIME_pipeline
Determining time of infection for human immunodeficiency viruses

## Installation and set up

### Install nextflow
`conda create -n time_analysis nextflow`

### Set up container image
A container image is available from Docker at `talnor/hiv_time_analysis`.
By default, this image will be used when using the docker or singularity profiles.

Alternatively, the container can be manually downloaded, or rebuilt from the Dockerfile in this 
repo. If so update the settings in the `nextflow.config`.
* manually pull singularity image: `singularity pull path/to/hiv_time_analysis.sif docker://talnor/hiv_time_analysis:<version>`
* manually pull docker image: `docker pull talnor/hiv_time_analysis:<version>`

## Usage
`conda activate time_analysis`

`nextflow run main.nf -profile slurm,singularity --input 'path/to/*_R{1,2}.fastq.gz' --outdir path/to/results/ --ticket <batch_name>`

## Options