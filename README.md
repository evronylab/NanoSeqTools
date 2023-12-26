# **NanoSeqTools**

This R package contains functions for analysis of NanoSeq data that was processed with the standard [NanoSeq pipeline](https://github.com/cancerit/NanoSeq).

### Outline
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [load_nanoseq_data](#load_nanoseq_data)
  - [load_nanoseq_regions](#load_nanoseq_regions)
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
Load NanoSeq data for genome-wide analysis.

##### Arguments
* dirs: A character vector of the directories containing NanoSeq results to load (one directory per sample).

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

* mutation_burden, with one row per sample. All substitution mutation statistics include all mutations, not just unique mutations.
  - The number of observed and corrected substitution mutations (muts_observed and muts_corrected)
  - Number of observed indels (indels_observed)
  - Total number of observed and corrected interrogated bases (total_observed and total_corrected; note: observed and corrected are the same)
  - Observed and corrected substitution mutation burdens (burden_observed and burden_corrected)
  - Observed indel mutation burden (burden_indels_observed)
  - Observed and corrected lower and upper confidence intervals of substitution mutation counts and observed lower and upper confidence intervals of indel counts (muts_lci_observed, muts_lci_corrected, indels_lci_observed, muts_uci_observed, muts_uci_corrected, indels_uci_observed)
  - Observed and corrected lower and upper confidence intervals of substitution mutation burdens and observed lower and upper confidence intervals of indel mutation burdens (burden_lci_observed, burden_lci_corrected, burden_indels_lci_observed, burden_uci_observed, burden_uci_corrected, burden_indels_uci_observed)

* purine_trinuc_mismatches: Data frame of the number of single-strand consensus purine mismatches. Columns: sample, tri (trinucleotide context), value.

* pyrimidine_trinuc_mismatches: Data frame of the number of single-strand consensus pyrimidine mismatches. Columns: sample, tri (trinucleotide context), value.

* estimated_error_rates: Data frame of the probability of having independent errors affecting both strands and resulting in a false-positive double-strand mutation and the number of estimated false positive double-strand mutations, based on the independent error rates in the purine channels. Columns: sample, total_error_rate, total_errors.

### load_nanoseq_regions
Load NanoSeq data for region-specific analysis.

##### Arguments
* nanoseq_data: Dataset resulting from load_nanoseq_data function

* regions.list: GRangesList object, comprised of GRanges that each contains a 'region set' to jointly analyze. The regions within each 'region set' can have overlaps (the functions handle this). The strand of each region in the region set specifies for each region, which mutations to include: '+ and '-' strand include mutations where central pyrimidine is on the '+' and '-' strands of the reference genome, respectively, and '*' includes all mutations. When there are overlapping regions with opposite strands, the mutations are counted only once regardless of the strand. Best practice is to name the elements of regions.list, since these names are carried forward to the output.

* excluded_samples: Names of samples excluded from the results because they do not have NanoSeq read coverage in any regions

* tabix_bin: Full path of tabix binary

#### Outputs
* sample_names: A vector of all sample IDs that were loaded

* dir: A vector of the directories containing the NanoSeq results that were loaded

* regions.list: Copy of input regions.list

* trinuc_bg_counts_ratio: Data frame of the sample trinucleotide background counts (i.e. number of interrogated bases for each trinucleotide context), the genome trinucleotide background counts (i.e. number of each trinucleotide context), and the normalized ratio of these for each sample/region combination. Columns: sample, region, tri (trinucleotide context), sample_tri_bg, genome_tri_bg, ratio2genome.

* trinuc_bg_counts.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of the sample trinucleotide background counts, with one row per sample and one column per trinucleotide context.

* trinuc_bg_ratio.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of the ratio of the sample trinucleotide background counts (normalized to a sum of 1) to the genome trinucleotide background counts (normalized to a sum of 1), with one row per sample and one column per trinucleotide context.

* genome_trinuc_counts.sigfit: Vector of the genome trinucleotide background counts, in the same order as columns in sigfit format columns.

* observed_corrected_trinuc_counts: Data frame of observed and corrected mutation counts (for all mutations and for unique mutations) for each sample/region combination. Columns: sample, region, tri (trinucleotide context), trint_subst_observed, trint_subst_unique_observed, ratio2genome, trint_subst_corrected, trint_subst_unique_corrected.

* observed_trinuc_counts.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of unique observed mutation counts, with one row per sample and one column per trinucleotide substitution context.

* mutation_burden: Data frame of the total number of observed and corrected mutations, total number of observed and corrected interrogated bases (note: observed and corrected are the same), observed and corrected mutation burdens, observed and corrected lower and upper confidence intervals of mutation counts, and observed and corrected lower and upper confidence intervals of mutation burdens, for each sample/region combination. All these statistics include all mutations, not just unique mutations.

## Citation
If you use NanoSeqTools, please cite Srinivasa A and Evrony GD. (2024). NanoSeqTools [Computer software]. https://github.com/evronylab/NanoSeqTools

NanoSeqTools also incorporates code from [indelwald](https://github.com/MaximilianStammnitz/Indelwald). Citation: The evolution of two transmissible cancers in Tasmanian devils (Stammnitz et al. 2023, Science 380:6642)