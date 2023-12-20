#' Shared functions for NanoSeqTools
#' @NoRD

#Load packages
suppressPackageStartupMessages(library(tidyverse))

# Order of trinucleotide context labels
trint_subs_labels <- c("ACA>AAA","ACC>AAC","ACG>AAG","ACT>AAT","CCA>CAA","CCC>CAC","CCG>CAG","CCT>CAT","GCA>GAA","GCC>GAC","GCG>GAG","GCT>GAT","TCA>TAA","TCC>TAC","TCG>TAG","TCT>TAT","ACA>AGA","ACC>AGC","ACG>AGG","ACT>AGT","CCA>CGA","CCC>CGC","CCG>CGG","CCT>CGT","GCA>GGA","GCC>GGC","GCG>GGG","GCT>GGT","TCA>TGA","TCC>TGC","TCG>TGG","TCT>TGT","ACA>ATA","ACC>ATC","ACG>ATG","ACT>ATT","CCA>CTA","CCC>CTC","CCG>CTG","CCT>CTT","GCA>GTA","GCC>GTC","GCG>GTG","GCT>GTT","TCA>TTA","TCC>TTC","TCG>TTG","TCT>TTT","ATA>AAA","ATC>AAC","ATG>AAG","ATT>AAT","CTA>CAA","CTC>CAC","CTG>CAG","CTT>CAT","GTA>GAA","GTC>GAC","GTG>GAG","GTT>GAT","TTA>TAA","TTC>TAC","TTG>TAG","TTT>TAT","ATA>ACA","ATC>ACC","ATG>ACG","ATT>ACT","CTA>CCA","CTC>CCC","CTG>CCG","CTT>CCT","GTA>GCA","GTC>GCC","GTG>GCG","GTT>GCT","TTA>TCA","TTC>TCC","TTG>TCG","TTT>TCT","ATA>AGA","ATC>AGC","ATG>AGG","ATT>AGT","CTA>CGA","CTC>CGC","CTG>CGG","CTT>CGT","GTA>GGA","GTC>GGC","GTG>GGG","GTT>GGT","TTA>TGA","TTC>TGC","TTG>TGG","TTT>TGT")

genome_freqs_labels <- str_sub(trint_subs_labels,1,3)

#All possible trinucleotides
trinucleotides_64 <- apply(expand.grid(c("A","C","G","T"),c("A","C","G","T"),c("A","C","G","T")),1,paste,collapse="")
trinucleotides_32_pyr <- apply(expand.grid(c("A","C","G","T"),c("C","T"),c("A","C","G","T")),1,paste,collapse="")
trinucleotides_32_pur <- setdiff(trinucleotides_64,trinucleotides_32_pyr)

#Function to reduce 64 to 32 trinucleotide frequency with central pyrimidine. Input is integer array with named elements that results from the trinucleotideFrequency function of Biostrings.
trinucleotide64to32 <- function(x){
	trinucleotides_64 <- apply(expand.grid(c("A","C","G","T"),c("A","C","G","T"),c("A","C","G","T")),1,paste,collapse="")
	trinucleotides_32_pyr <- apply(expand.grid(c("A","C","G","T"),c("C","T"),c("A","C","G","T")),1,paste,collapse="")
	trinucleotides_32_pur <- setdiff(trinucleotides_64,trinucleotides_32_pyr)
	
	y <- x[trinucleotides_32_pur]
	x <- x[trinucleotides_32_pyr]
	names(y) <- reverseComplement(DNAStringSet(names(y)))
	result <- merge(x,y,by=0)
	row.names(result) <- result[,1]
	return(apply(result[,-1],1,sum))
}