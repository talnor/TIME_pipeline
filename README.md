# TIME_pipeline
Determining time of infection for human immunodeficiency viruses.

The analysis builds upon the [Shiver tool](https://github.com/ChrisHIV/shiver) for mapping paired-end short reads to a 
custom reference sequence constructed using do novo assembled contigs. Base frequencies from the alignment are then
used to calculate the time of infection based on the accumulated mutations in the pol gene as described in the 
publication by Puller et al.

> Neher R, Albert J (2017), Estimating time of HIV-1 infection from next-generation sequence diversity.
PLoS Comput Biol 13(10): e1005775. https://doi.org/10.1371/journal.pcbi.1005775

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
| hostGenome                | Directory with human reference genome database (bwa)            |
| hostGenomeBase            | The name (base) of the database files                           |
| cache                     | Directory for cache files                                       |
| process.clusterOptions    | Cluster options if using slurm for execution                    |

In addition, default settings for primers, adapters and similar configurations are described in more 
detail [here](data/README.md).

### Download host reference genome
Set `nextflow.config` parameter `hostFasta` to the path of the host reference genome of your choice. 
Then set up the host database with the following command. The database will be placed in the `hostGenome` directory
and will be namned as `hostGenomeBase`.
```
nextflow run main.nf --setup -profile slurm,singularity --outdir <outdir>
```

### Run Shiver initialisation
Shiver initilisation directories are included in this repository. Information on these are 
available [here](data/README.md). To create your own initilisation directory, run the following command:
```
nextflow run main.nf --init -profile slurm,singularity --primers <primers.fasta> --adapters <adapters.fasta> --config <shiver_config.sh> --references <references.fasta> --outdir <outdir>
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