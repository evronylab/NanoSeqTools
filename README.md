# **NanoSeqTools**

This R package contains functions for NanoSeq data analysis. For more information regarding NanoSeq output files, refer to the "output" section in the [NanoSeq GitHub](https://github.com/cancerit/NanoSeq).

### Outline
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [load_nanoseq_data](#load_nanoseq_data)
  - [load_nanoseq_bedcov](#load_nanoseq_bedcov)
- [Citation](#citation)

## Prerequisites
* [Tidyverse](https://www.tidyverse.org/packages/)

* [vcfR](https://github.com/knausb/vcfR)

* [rtracklayer](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html)

* [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)

* BSgenome package corresponding to the reference genome used in the NanoSeq analysis

* [tabix](http://www.htslib.org/download/) >=v1.10; part of htslib

## Installation

```
## Install NanoSeqTools directly from github using devtools

install.packages("devtools")
devtools::install_github('https://github.com/evronylab/NanoSeqTools/')
```

## Usage
### load_nanoseq_data
Load all files for analysis of NanoSeq data
```
load_nanoseq_data(dirs, suffix_to_remove, BSgenomepackagename, BSgenomecontigs)
```

##### Arguments
* dirs: List of all the directories with NanoSeq results to load.

The list of directories can be generated with:
```
list.dirs(normalizePath(directory_path), full.names = TRUE, recursive = FALSE)
```
where directory_path is a parent directory containing all the directories.

* sample_names: A character vector of the sample names to assign to the results, in the same order as the directories in 'dirs'.

* BSgenomepackagename: A string of a BSgenome package corresponding to the reference genome used in the Nanoseq analysis (e.g., "BSgenome.Hsapiens.UCSC.hg38" or "BSgenome.Mmusculus.UCSC.mm10"). This is used to calculate the genome trinucleotide background and corrected substitution counts and burdens.

* BSgenomecontigs: A vector of numeric indices of the contigs of the BSgenome package from which to calculate the genome trinucleotide background (e.g., 1:24 for BSgenome.Hsapiens.UCSC.hg38, or 1:21 for BSgenome.Mmusculus.UCSC.mm10).

#### Outputs
A list containing the following data objects:

* sample_names: A vector of all sample IDs that were loaded

* dirs: A vector of the directories containing the NanoSeq results that were loaded

* vcf_snp.fix: List (one object per sample) containing the fixed information (fix) from the SNP vcf (only FILTER = PASS mutations)

* vcf_indel.fix: List (one object per sample) containing the fixed information (fix) from the indel vcf (only FILTER = PASS mutations)

* vcf_indel.gt: List (one object per sample) containing the genotype information (gt) from the indel vcf (only FILTER = PASS mutations)

* trinuc_bg_counts_ratio: Data frame of the sample trinucleotide background counts (i.e. number of interrogated bases for each trinucleotide context), the genome trinucleotide background counts (i.e. number of each trinucleotide context), and the normalized ratio of these. Columns: sample, tri (trinucleotide context), sample_tri_bg, genome_tri_bg, ratio2genome.

* trinuc_bg_counts.sigfit: Data frame in sigfit format of the sample trinucleotide background counts, with one row per sample and one column per trinucleotide context.

* trinuc_bg_ratio.sigfit: Data frame in sigfit format of the ratio of the sample trinucleotide background counts (normalized to a sum of 1) to the genome trinucleotide background counts (normalized to a sum of 1), with one row per sample and one column per trinucleotide context.

* genome_trinuc_counts.sigfit: Vector of the genome trinucleotide background counts, in the same order as columns in sigfit format columns.

* observed_corrected_trinuc_counts: Data frame of observed and corrected mutation counts (for all mutations and for unique mutations). Columns: sample, tri (trinucleotide context), trint_subst_observed, trint_subst_unique_observed, ratio2genome, trint_subst_corrected, trint_subst_unique_corrected.

* observed_trinuc_counts.sigfit: Data frame in sigfit format of unique observed mutation counts, with one row per sample and one column per trinucleotide substitution context.

* mutation_burden: The total number of observed and corrected mutations, total number of observed and corrected interrogated bases (note: observed and corrected are the same), observed and corrected mutation burdens, observed and corrected lower and upper confidence intervals of mutation counts, and observed and corrected lower and upper confidence intervals of mutation burdens, with one row per sample. All these statistics include all mutations, not just unique mutations.

* purine_trinuc_mismatches: Data frame of the number of single-strand consensus purine mismatches. Columns: sample, tri (trinucleotide context), value.

* pyrimidine_trinuc_mismatches: Data frame of the number of single-strand consensus pyrimidine mismatches. Columns: sample, tri (trinucleotide context), value.

* estimated_error_rates: Data frame of the probability of having independent errors affecting both strands and resulting in a false-positive double-strand mutation and the number of estimated false positive double-strand mutations, based on the independent error rates in the purine channels. Columns: sample, total_error_rate, total_errors.

### load_nanoseq_regions
Load all files for analysis of NanoSeq data

##### Arguments
* 

#### Outputs
* 

## Citation
If you use NanoSeqTools, please cite Srinivasa A, Evrony GD, [TBD].
