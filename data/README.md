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
| Version     | Primers     | Description | 
| ----------- | ----------- | ----------- |
| 1 | Primers_A_elife-11282-supp2-v2_PCR1-2_primers_A_primers_RC.fasta |        |
| 2 | Primers_A_elife-11282-supp2-v2_PCR1_primers_A_primers_RC.fasta |        |
| 3 | Primers_A_elife-11282-supp2-v2_PCR1_primers_A_primers_RC_Bprimers_fragment4.fasta |        |
| 4 | Primers_A_elife-11282-supp2-v2_PCR2_primers_A_primers.fasta |        |
| 5 | primers_1_amplicon_PCR1-2_190620.fasta |        |
| 6 | primers_1_amplicon_PCR1_190620.fasta | Full genome amplified with 1 primer pair |
| 7 | primers_1_amplicon_PCR2_190620.fasta |        |
| 8 | primers_B1_180119.fasta |        |
| 9 | primers_B1_201203.fasta |        |
| 10 | primers_B1_201203_PCR1.fasta |        |

### Adapters
| Version     | Adapters    | Description           | 
| ----------- | ----------- | --------------------- |
| 1 | NexteraPE-PE.fa | Nextera paired-end adapters | 

### Shiver configuration file
| Version     | Configurations | Description | 
| ----------- | ----------- | ----------- |
| 1 | original_config.sh | Default settings used in Shiver |
| 2 | shiver_config_BQ30_notrimming.sh | TIME-study settings |
| 3 | config_BQ30.sh | Older settings |

#### Configuration file 2

The following options in Shiver are altered. For the full list of options see the default config.

| Parameter   | Value       | Default     | Description |
| ----------- | ----------- | ----------- | ----------- |
| TrimReadsForAdaptersAndQual      | false | true | Trim adapaters and low quality bases from reads using trimmomatic? |
| TrimReadsForPrimers      | false | true | Trim exact matches to PCR primers from the end of reads using fastaq? |
| mpileupOptions      | --min-BQ 30 | --min-BQ 5 | Higher quality threshold for individual bases |
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
| 1 | InitDirShiver220128_BQ30_1amp | 1 amplicon primers, 2020 references |
| 2 | InitDirShiver190405_BQ30 |  |
| 3 | InitDirShiver191022_BQ30_PANHIV | 1 amplicon primers, 2018 references |

#### Shiver init directory 1

**Name**: InitDirShiver220128_BQ30_1amp 
**Created**: 2022-01-28

| Content     | Description     |
| ----------- | ----------- |
| Primer      | primers_1_amplicon_PCR1_190620.fasta | 
| Adapter      | NexteraPE-PE.fa | 
| Configurations      | shiver_config_BQ30_notrimming.sh |
| References      | HIV1_COM_2020_genome_DNA.fasta' |

#### Shiver init directory 2

**Name**: InitDirShiver190405_BQ30
**Created**: 2019-04-05

#### Shiver init directory 3

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
| HIV1_COM_2017_547-9592_DNA_2018Compendium.fasta | Represenative genome alignment with references from 2018 and earlier. Genomic positions 547-9592 included. |

