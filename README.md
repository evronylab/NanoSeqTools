# **NanoSeqTools**

# _R package to load all relevant files for NanoSeq data analysis_


## Description
This function loads all the files for analysis of NanoSeq data. For more information regarding NanoSeq output files, refer to the "output" section in the [NanoSeq GitHub](https://github.com/cancerit/NanoSeq).

## Prerequisites
* [Tidyverse](https://www.tidyverse.org/packages/)

* [vcfR](https://github.com/knausb/vcfR)

* BSgenome package corresponding to the reference genome used in the NanoSeq analysis

## Installation

```
## Install NanoSeqTools directly from github using devtools

install.packages("devtools")
devtools::install_github('https://github.com/evronylab/NanoSeqTools/')
```

## Usage
library(vcfR)
library(tidyverse)

load_nanoseq_data(dirs, suffix_to_remove, BSgenomepackagename, BSgenomecontigs)

## Arguments
### dirs	
List of all the directories with NanoSeq results to load.

The list of directories can be generated with:
```
list.dirs(normalizePath(directory_path), full.names = TRUE, recursive = FALSE)
```
where directory_path is a parent directory containing all the directories.

### suffix_to_remove	
A string to remove from the NanoSeq results directory names, to retain only the Sample IDs without any suffixes. If there is no suffix to remove, give the input as "NA".

### BSgenomepackagename	
A string of a BSgenome package corresponding to the reference genome used in the Nanoseq analysis (e.g., "BSgenome.Hsapiens.UCSC.hg38" or "BSgenome.Mmusculus.UCSC.mm10"). This is used to calculate the genome trinucleotide background and corrected substitution counts and burdens.

### BSgenomecontigs	
A vector of numeric indices of the contigs of the BSgenome package from which to calculate the genome trinucleotide background (e.g., 1:24 for BSgenome.Hsapiens.UCSC.hg38, or 1:21 for BSgenome.Mmusculus.UCSC.mm10).

## Outputs
A list containing the following data objects:

* sample_id: A vector of all sample IDs that were loaded

* vcf_snp.fix: List (one object per sample) containing the fixed information (fix) from the SNP vcf (only FILTER = PASS mutations)

* vcf_indel.fix: List (one object per sample) containing the fixed information (fix) from the indel vcf (only FILTER = PASS mutations)

* vcf_indel.gt: List (one object per sample) containing the genotype information (gt) from the indel vcf (only FILTER = PASS mutations)

* trinuc_bg_ratio: Data frame of the sample trinucleotide background (i.e. number of interrogated bases for each trinucleotide context), the genome trinucleotide background (i.e. number of each trinucleotide context), and the normalized ratio of these. Columns: sample, tri (trinucleotide context), sample_tri_bg, genome_tri_bg, ratio2genome.

* trinuc_bg.sigfit: Data frame in sigfit format of the trinucleotide background (i.e. number of interrogated bases for each trinucleotide context), with one row per sample and one column per trinucleotide context.

* observed_corrected_trinuc_counts: Data frame of observed and corrected mutation counts. Columns: sample, tri (trinucleotide context), trint_subst_observed, trint_subst_corrected.

* observed_trinuc_counts.sigfit: Data frame in sigfit format of observed mutation counts, with one row per sample and one column per trinucleotide substitution context.

* mutation_burden: The total number of observed and corrected mutations, total number of observed and corrected interrogated bases (note: observed and corrected are the same), observed and corrected mutation burdens, observed and corrected lower and upper confidence intervals of mutation counts, and observed and corrected lower and upper confidence intervals of mutation burdens, with one row per sample.

* purine_trinuc_mismatches: Data frame of the number of single-strand consensus purine mismatches. Columns: sample, tri (trinucleotide context), value.

* pyrimidine_trinuc_mismatches: Data frame of the number of single-strand consensus pyrimidine mismatches. Columns: sample, tri (trinucleotide context), value.

* estimated_error_rates: Data frame of the probability of having independent errors affecting both strands and resulting in a false-positive double-strand mutation and the number of estimated false positive double-strand mutations, based on the independent error rates in the purine channels. Columns: sample, total_error_rate, total_errors.