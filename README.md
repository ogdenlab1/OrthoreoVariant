# OrthoreoVariant

Code related to viral single nucleotide variant detection in Illumina datasets.

## Description

Code utilized in Ogden lab to analyze viral SNVs found in Illumina datasets.

## Getting Started

### Dependencies

* Analysis was performed on an Ubuntu 20 OS. Please note your available computational resourcing and adjust accordingly.
* Pre-processing: Trimmomatic-0.39
* Alignment: bowtie2, samtools, bbmap
* Variant calling: LoFreq
* Python 3.11+
* Python packages: pandas, os, fnmatch, argparse, pathlib

### Installing

* Code can be downloaded directly from this repo. Bash commands can be adjusted as necessary

### Preprocessing, alignment, and variant calling

* Bash script assumes Trimmomatic and BBMap are installed in the home directory
* Samtools, bowtie2, and LoFreq must be installed via an environment variable
* Script must be run in directory with raw data and reference files and reference Bowtie2 index.
* FASTQ extensions are specific to VUMC Vantage core. You may have to change them to work with your system.
* samples.txt is a file with each sample name on a line and present in the raw data file where the script is run. NOTE: If you are wanting to use an alias, this script does not rename the files.
```
sh ./VirVariant_alignment.sh SA11 SA11.fa
```

### Variant processing

* Bash script assumes Trimmomatic and BBMap are installed in the home directory
* Samtools, bowtie2, and LoFreq must be installed via an environment variable
* Script must be run in directory with raw data and reference files and reference Bowtie2 index.
* FASTQ extensions are specific to VUMC Vantage core. You may have to change them to work with your system.
```
python ./VirVariant.py samples.txt ~/Projects/SA11_variants/ SA11_analysis
```

## Help

```
python ./VirVariant.py -h
usage: VirVariant.py [-h] [--freq FREQ] [--file_tag FILE_TAG] Sample_List Working_Directory Experiment

positional arguments:
  Sample_List          A tab delineated file with each line containing sample names.
  Working_Directory    Path to directory containing data to align.
  Experiment           Experiment name.

options:
  -h, --help           show this help message and exit
  --freq FREQ          Variant frequency cutoff for filtering. Decimal between 0 and 1.
  --file_tag FILE_TAG  File naming tag can denote filters used for variant isolation.
```

## Authors

Contributors names and contact info

Jennifer Gribble (Bowser)  
jgribble.bowser@gmail.com
