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

### Configure pipeline options in the nextflow config file
Default parameters and settings for running the pipeline are specified in `nextflow.config`. 
Important settings to change include:

| Settings                  | Description                                                     | 
| ------------------------- | --------------------------------------------------------------- |
| hostGenome                | Directory with human reference genome database, indexed for bwa |
| hostGenomeBase            | The name (base) of the database files                           |
| cache                     | Directory for cache files                                       |
| process.clusterOptions    | Cluster options if using slurm for execution                    |

In addition, default settings for primers, adapters and similar configurations are described in more 
detail [here](data/README.md).

### Run Shiver initialisation
Shiver initilisation directories are included in this repository. Information on these are 
available [here](data/README.md). To create your own initilisation directory, run the following command:
```
nextflow run main.nf --init -profile slurm,singularity --primers <primers.fasta> --adapters <adapters.fasta> --config <shiver_config.sh> --references <references.fasta>
```

## Usage
`conda activate time_analysis`

Basic usage:

```
nextflow run main.nf -profile slurm,singularity --input 'path/to/*_R{1,2}.fastq.gz' --outdir path/to/results/ --ticket <batch_name>
```

The pipeline can be executed on your **local** computer or with a **slurm** resource manager. The container can be run using 
**docker** or **singularity**. The above command would run the pipeline using slurm and singularity. 

## Options
Check the command help for more info and options.
```
nextflow run main.nf --help
```