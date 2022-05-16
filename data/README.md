# Configurations for the analysis

## Parameters to set

To run the analysis four file paths and directories are required as input.
The default files are listed in the 'nextflow.config' and can be changed with the following parameters added to the 
run command:

| File        | Command     |
| ----------- | ----------- |
| Primer      | --primers "file path" |
| Adapter      | --adapters "file path" | 
| InitDir      | --initDir "file path" |
| Configurations      | --config "directory path" | 

## Description of files

Below follows a description of the files that are made available in this directory.

### Primers
| Version     | File  | Primer type   | Description | 
| ----------- | ----------- | ----------- | ----------- |
| 1 | A_primers_6AMP_PCR1-2.fasta | A | 6 primer pairs, PCR1+PCR2 | 
| 2 | A_primers_6AMP_PCR1.fasta | A | 6 primer pairs, PCR1 | 
| 3 | A_primers_6AMP_PCR1_F4-Bprimers.fasta | A+B | 6 primer pairs, PCR1, B primers used for fragment 4 |
| 4 | B_primers_1AMP_PCR1-2.fasta | B | 1 primer pair, PCR1+PCR2 | 
| 5 | AB_primers_1AMP_PCR1.fasta | A, B | 1 primer pair, PCR1, A and B primers are identical |
| 6 | B_primers_6AMP_PCR1-2.fasta | B | 6 primer pairs, PCR1+PCR2 |
| 7 | B_primers_6AMP_PCR1.fasta | B | 6 primer pairs, PCR1 |

### Adapters
| Version     | Adapters    | Description           | 
| ----------- | ----------- | --------------------- |
| 1 | NexteraPE-PE.fa | Nextera paired-end adapters | 

### Shiver configuration file
| Version     | Configurations | Description | 
| ----------- | ----------- | ----------- |
| 1 | original_config.sh | Default settings used in Shiver |
| 2 | shiver_config_BQ20_notrimming.sh | TIME-study settings |
| 3 | config_BQ30.sh | Older settings |

#### Configuration file 2

The following options in Shiver are altered. For the full list of options see the default config.

| Parameter   | Value       | Default     | Description |
| ----------- | ----------- | ----------- | ----------- |
| TrimReadsForAdaptersAndQual      | false | true | Trim adapaters and low quality bases from reads using trimmomatic? |
| TrimReadsForPrimers      | false | true | Trim exact matches to PCR primers from the end of reads using fastaq? |
| mpileupOptions      | --min-BQ 20 | --min-BQ 5 | Higher quality threshold for individual bases |
| deduplicate      | true | false | Remove read pairs marked as duplicates? This can cause loss of diversity in the reads due to true biological variation as well sequencing error. |

#### Configuration file 3

The following options in Shiver are altered. For the full list of options see the default config.

| Parameter   | Value       | Default     | Description |
| ----------- | ----------- | ----------- | ----------- |
| mpileupOptions      | --min-BQ 30 | --min-BQ 5 | Higher quality threshold for individual bases |
| deduplicate      | true | false | Remove read pairs marked as duplicates? This can cause loss of diversity in the reads due to true biological variation as well sequencing error. |


### Shiver init directory
| Version     | InitDir     | Description | 
| ----------- | ----------- | ----------- |
| 1 | InitDirShiver220223_BQ30_1amp | 1 amplicon primers, 2020 references (no UTRs) |
| 2 | InitDirShiver220128_BQ30_1amp | 1 amplicon primers, 2020 references |
| 3 | InitDirShiver190405_BQ30 |  |
| 4 | InitDirShiver191022_BQ30_PANHIV | 1 amplicon primers, 2018 references |
| 5 | InitDirShiver220516_BQ20_1AMP_Bprimers_PCR1 | 1 amplicon primers, 2020 references (no UTRs) |
| 6 | InitDirShiver220516_BQ20_1AMP_Bprimers_PCR2 | 1 amplicon primers, 2020 references (no UTRs) |
| 7 | InitDirShiver220516_BQ20_6AMP_ABprimers_PCR1 | 6 amplicon primers, 2020 references (no UTRs) |
| 8 | InitDirShiver220516_BQ20_6AMP_Aprimers_PCR1-2 | 6 amplicon primers, 2020 references (no UTRs) |
| 9 | InitDirShiver220516_BQ20_6AMP_Aprimers_PCR1 | 6 amplicon primers, 2020 references (no UTRs) |

#### Shiver init directory 1

**Name**: InitDirShiver220223_BQ30_1amp 
**Created**: 2022-02-23

| Content     | Description     |
| ----------- | ----------- |
| Primer      | primers_1_amplicon_PCR1_190620.fasta | 
| Adapter      | NexteraPE-PE.fa | 
| Configurations      | shiver_config_BQ30_notrimming.sh |
| References      | HIV1_COM_2020_547-9592_DNA.fasta |

#### Shiver init directory 2

**Name**: InitDirShiver220128_BQ30_1amp 
**Created**: 2022-01-28

| Content     | Description     |
| ----------- | ----------- |
| Primer      | primers_1_amplicon_PCR1_190620.fasta | 
| Adapter      | NexteraPE-PE.fa | 
| Configurations      | shiver_config_BQ30_notrimming.sh |
| References      | HIV1_COM_2020_genome_DNA.fasta |

#### Shiver init directory 3

**Name**: InitDirShiver190405_BQ30
**Created**: 2019-04-05

#### Shiver init directory 4

**Name**: InitDirShiver191022_BQ30_PANHIV 
**Created**: 2019-10-22

| Content     | Description     |
| ----------- | ----------- |
| Primer      | primers_1_amplicon_PCR1-2_190620.fasta | 
| Adapter      | NexteraPE-PE.fa | 
| Configurations      | config_BQ30.sh |
| References      | HIV1_COM_2017_547-9592_DNA_2018Compendium.fasta |

### References to use in Shiver alignments
Reference compendiums with representative genomes can be downloaded from the 
[LANL HIV database](http://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html).

| Reference  file                | Description                                                          |
| ------------------------------ | -------------------------------------------------------------------- |
| HIV1_COM_2020_genome_DNA.fasta | Represenative genome alignment with references from 2020 and earlier. |
| HIV1_COM_2020_547-9592_DNA.fasta | Represenative genome alignment with references from 2020 and earlier. Genomic positions 547-9592 included. |
| HIV1_COM_2017_547-9592_DNA_2018Compendium.fasta | Represenative genome alignment with references from 2018 and earlier. Genomic positions 547-9592 included. |

