# TIME_pipeline
Determining time of infection for human immunodeficiency viruses.

The analysis builds upon the [Shiver tool](https://github.com/ChrisHIV/shiver) for mapping paired-end short reads to a 
custom reference sequence constructed using do novo assembled contigs. Base frequencies from the alignment are then
used to calculate the time of infection based on the accumulated mutations in the pol gene as described in the 
publication by Puller et al.

> Neher R, Albert J (2017), Estimating time of HIV-1 infection from next-generation sequence diversity.
PLoS Comput Biol 13(10): e1005775. https://doi.org/10.1371/journal.pcbi.1005775

## Installation and set up

### Install required software

Make sure the following is installed, or install them:
- Singularity or Docker
- Nextflow or conda (install nextflow using conda as described below)

```
conda create -n time_analysis nextflow
conda activate time_analysis
```

### Install TIME_pipeline

```
git clone https://github.com/talnor/TIME_pipeline.git
```

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

For example, the following command can be run to download the human reference genome:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
```
 
Then set up the host database with the following command. The database will be placed in the `hostGenome` directory
and will be namned as `hostGenomeBase`.
```
nextflow run main.nf --setup -profile slurm,singularity --hostFasta <path_to_genome> --outdir <outdir>
```

### Run Shiver initialisation
The Shiver initialization directory includes the set of primers used during the amplification of the samples as well as a 
reference dataset to be used in the analysis. Several options are included in this repository. Information on these are 
available [here](data/README.md). To create your own initialization directory, run the following command:
```
nextflow run main.nf --init -profile slurm,singularity --primers <primers.fasta> --adapters <adapters.fasta> --config <shiver_config.sh> --references <references.fasta> --outdir <outdir>
```

## Usage
- Ensure the settings in the `nextflow.config` are correct for your samples. Importantly, the primer set and the 
  initialization directory needs to match the primers used during amplification of the samples. Or override the
  default values by supplying them as parameters in the command below.

Basic usage:
```
conda activate time_analysis
nextflow run main.nf -profile slurm,singularity --input 'path/to/*_R{1,2}.fastq.gz' --outdir path/to/results/ --ticket <batch_name>
```

The pipeline can be executed on your **local** computer or with a **slurm** resource manager. The container can be run using 
**docker** or **singularity**. The above command would run the pipeline using slurm and singularity. 

## Options
Check the command help for more info and options.
```
nextflow run main.nf --help
```

## Optional installation steps

### Set up container image
A container image is available from Docker at `talnor/hiv_time_analysis`.
By default, this image will be used when using the docker or singularity profiles. If this doesn't work, the container can be manually downloaded, or rebuilt from the Dockerfile in this 
repo. If so update the settings in the `nextflow.config`.
* manually pull singularity image: `singularity pull path/to/hiv_time_analysis.sif docker://talnor/hiv_time_analysis:<version>`
* manually pull docker image: `docker pull talnor/hiv_time_analysis:<version>`