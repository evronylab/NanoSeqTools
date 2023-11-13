# **NanoSeqTools**

# _R package to load all relevant files for NanoSeq data analysis_


## Description
This function loads all the relevant files one would need for analysis on NanoSeq data. For more information regarding the files refer to the "output" section in the NanoSeq GitHub here: [Link to NanoSeq GitHub](https://github.com/cancerit/NanoSeq)

The first input for this function would be the list of directories containing the NanoSeq results. The function gives a data structure as the output which has all the data frames for Mutation Burden analysis, Single base Substitutions and Indel Signature analysis in the results object as a list.

The second input of this function would be the suffix at the end of the directory name that the user wants to remove, so as to retain only the Sample IDs. For example, if the name of the directory is "sample_id_results", running the function with suffix_to_remove="_results", would remove this suffix from all file names and retain only "sample_id". If no suffix is to be removed, the input is given as "NA".

## Installation

```

## Install NanoSeqTools directly from github using devtools

install.packages("devtools")
devtools::install_github('https://github.com/evronylab/NanoSeqTools/')

```

## Usage
load_nanoseq_data(subdirs, suffix_to_remove)

## Arguments
### subdirs	
List of all directories in which the results of NanoSeq is saved. Generate the list of sub directories using list.dirs(normalizePath(directory_path), full.names = TRUE, recursive = FALSE) where directory_path is the parent directory with all the sub directories.

### suffix_to_remove	
A string input of the suffix which you want removed from the name of the directories, to retain only the Sample IDs without any suffixes. If there is no suffix to remove, give the input as "NA".

## Output
A list of data files with all the relevant results generated by NanoSeq:

* sample_id: A vector of all Sample_IDs(file names) from which the results are being loaded

* observed_trinuc_counts: The observed mutation counts with one row per sample and one column per trinucleotide substitution

* corrected_trinuc_counts: The corrected mutation counts with one row per sample and one column per trinucleotide substitution

* trinuc_bg: The trinucleotide backgrounds of each mutation with one row per sample and one column per trinucleotide substitution

* opportunities_ratio: The opportunities ratio - normalized bakgrounds/normalized human trinucleotide frequencies - of each mutation with one row per sample and one column per trinucleotide substitution

* mutation_burden: The total number of observed/corrcted mutations, total trinucleotide backgrounds, observed /corrected mutation burdens, observed/corrected lower and upper confidence of mutation burdens with one row per sample

* indel_vcf_fix: The fixed information (fix) from the indel vcf files of all the samples as a subset within this object

* indel_vcf_gt: The genotype information (gt) from the indel vcf files of all the samples as a subset within this object

* snp_vcf_fix: The fixed information (fix) from the SNP vcf files of all the samples as a subset within this object

* snp_vcf_gt: The genotype information (gt) from the SNP vcf files of all the samples as a subset within this object

* purine_trinuc_mismatches: The probability of having independent errors affecting both strands and resulting in double-strand consensus, based on the independent error rates in the purine channels with one row per sample and one column per purine mismatch

* pyrimidine_trinuc_mismatches: The probability of having independent errors affecting both strands and resulting in double-strand consensus, based on the independent error rates in the pyrimidine channels with one row per sample and one column per pyrimidine mismatch
