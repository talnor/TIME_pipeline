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
Primers used when?
| Version     | Primers     | Description | 
| ----------- | ----------- | ----------- |
| 1 | Primers_A_elife-11282-supp2-v2_PCR1-2_primers_A_primers_RC.fasta | Description       |
| 2 | Primers_A_elife-11282-supp2-v2_PCR1_primers_A_primers_RC.fasta | Description       |
| 3 | Primers_A_elife-11282-supp2-v2_PCR1_primers_A_primers_RC_Bprimers_fragment4.fasta | Description       |
| 4 | Primers_A_elife-11282-supp2-v2_PCR2_primers_A_primers.fasta | Description       |
| 5 | primers_1_amplicon_PCR1-2_190620.fasta | Description       |
| 6 | primers_1_amplicon_PCR1_190620.fasta | Full genome amplified with 1 primer pair |
| 7 | primers_1_amplicon_PCR2_190620.fasta | Description       |
| 8 | primers_B1_180119.fasta | Description       |
| 9 | primers_B1_201203.fasta | Description       |
| 10 | primers_B1_201203_PCR1.fasta | Description       |

### Adapters
Adaptors used when?
| Version     | Adapters    | Description | 
| ----------- | ----------- | ----------- |
| 1 | adapters.fasta | Description | 

### Shiver configuration file
| Version     | Configurations | Description | 
| ----------- | ----------- | ----------- |
| 1 | original_config.sh | Default settings used in Shiver |
| 2 | shiver_config_BQ30_notrimming.sh | Description |

#### Configuration file 1

The following options in Shiver are altered. For the full list of options see the default config.

| Parameter   | Value       | Default     | Description |
| ----------- | ----------- | ----------- | ----------- |
| TrimReadsForAdaptersAndQual      | false | true | Trim adapaters and low quality bases from reads using trimmomatic? |
| TrimReadsForPrimers      | false | true | Trim exact matches to PCR primers from the end of reads using fastaq? |
| mpileupOptions      | --min-BQ 30 | --min-BQ 5 | Higher quality threshold for individual bases |
| deduplicate      | true | false | Remove read pairs marked as duplicates? This can cause loss of diversity in the reads due to true biological variation as well sequencing error. |

### Shiver init directory
| Version     | InitDir     | Description | 
| ----------- | ----------- | ----------- |
| 1 | InitDirShiver190405_BQ30 | Description |

#### Shiver init directory 1

**Name**: InitDirShiver190405_BQ30
**Created**: 2019-04-05

| Content     | Version     | Description | 
| ----------- | ----------- | ----------- |
| Primer      | version | Description |
| Adapter      | version | Description | 
| Configurations      | version | Description |
| References      | version | Description |

#### Shiver init directory 2

**Name**: 
**Created**: 

| Content     | Version     | Description | 
| ----------- | ----------- | ----------- |
| Primer      | version | Description |
| Adapter      | version | Description | 
| Configurations      | version | Description |
| References      | version | Description |
