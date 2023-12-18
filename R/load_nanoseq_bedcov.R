#' Load NanoSeq coverage and trinculeotide spectrum for burden analysis of specific genomic regions
#' 
#' Package requirements: tidyverse, rtracklayer, GenomicRanges, tabix binary
#'
#' For more information regarding NanoSeq output files, refer to the "output" section in the [NanoSeq GitHub](https://github.com/cancerit/NanoSeq).
#'
#' @param nanoseq_data Dataset resulting from load_nanoseq_data function
#' @param regions.list GRangesList object, comprised of GRanges that each contains a 'region set' to jointly analyze. The regions within each 'region set' can have overlaps (the functions handle this). The strand of each region in the region set specifies for each region, which mutations to include: '+ and '-' strand include mutations where central pyrimidine is on the '+' and '-' strands of the reference genome, respectively, and '*' includes all mutations. When there are overlapping regions with opposite strands, the mutations are counted only once regardless of the strand. Best practice is to name the elements of regions.list, since these names are carried forward to the output.
#' @param tabix_bin Full path to tabix binary
#' @return Returns NanoSeq results for each 'region set'.
#' * sample_names: A vector of all sample IDs that were loaded
#' * dir: A vector of the directories containing the NanoSeq results that were loaded
#' * regions.list: Copy of input regions.list
#' * trinuc_bg_counts_ratio: List with one object per region set, each comprised of a data frame of the sample trinucleotide background counts (i.e. number of interrogated bases for each trinucleotide context), the genome trinucleotide background counts (i.e. number of each trinucleotide context), and the normalized ratio of these. Columns: sample, tri (trinucleotide context), sample_tri_bg, genome_tri_bg, ratio2genome.
#' * trinuc_bg_counts.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of the sample trinucleotide background counts, with one row per sample and one column per trinucleotide context.
#' * trinuc_bg_ratio.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of the ratio of the sample trinucleotide background counts (normalized to a sum of 1) to the genome trinucleotide background counts (normalized to a sum of 1), with one row per sample and one column per trinucleotide context.
#' * genome_trinuc_counts.sigfit: Vector of the genome trinucleotide background counts, in the same order as columns in sigfit format columns.
#' * observed_corrected_trinuc_counts: List with one object per region set, each comprised of a data frame of observed and corrected mutation counts (for all mutations and for unique mutations). Columns: sample, tri (trinucleotide context), trint_subst_observed, trint_subst_unique_observed, ratio2genome, trint_subst_corrected, trint_subst_unique_corrected.
#' * observed_trinuc_counts.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of unique observed mutation counts, with one row per sample and one column per trinucleotide substitution context.
#' * mutation_burden: List with one object per region set, each comprised of a data frame with the total number of observed and corrected mutations, total number of observed and corrected interrogated bases (note: observed and corrected are the same), observed and corrected mutation burdens, observed and corrected lower and upper confidence intervals of mutation counts, and observed and corrected lower and upper confidence intervals of mutation burdens, with one row per sample. All these statistics include all mutations, not just unique mutations.
#' @export

load_nanoseq_regions <- function(nanoseq_data,regions.list,tabix_bin){

	#Load packages
	suppressPackageStartupMessages(library(tidyverse))
	suppressPackageStartupMessages(library(rtracklayer))
	suppressPackageStartupMessages(library(GenomicRanges))
	
	#Initialize lists for results
	bedcov <- list()
	trinuc_bg_counts_ratio <- list()
	observed_corrected_trinuc_counts <- list()
	mutation_burden <- list()
	
	dirs <- nanoseq_data$dirs
	
	
	for (i in 1:length(dirs)) {
		
		dir <- dirs[i]
		sample_name <- sample_names[i]
		
		message(paste0("  ",sample_name))
		
		#Load NanoSeq bed coverage data for each region after first reducing overlapping ranges.
		# Note, using tabix binary since it is much faster than R tabix tools.
		# Note, import function transforms bed coordinates to 1-based coordinates.
		tmp.regions.all <- tempfile()
		regions.list %>% unlist %>% reduce %>% export(con=tmp.regions.all,format="bed")
		
		tmp.bedcov.all <- tempfile()
		system(paste(tabix_bin,paste0(dir,"results.cov.bed.gz"),"-R",tmp.regions,"| sed 's/;/\t/g' | awk 'BEGIN{OFS=\"\t\"}{print $1,$2,$3,$6,$4,$5}'>",tmp.bedcov.all))
		
		bedcov.all <- import(tmp.bedcov.all,format="bedgraph")
		colnames(mcols(bedcov.all)) <- c("coverage","tri","ref")
		
		file.remove(tmp.regions.all,tmp.bedcov.all)
		
		bedcov[[sample_name]] <- findOverlaps(**)
		
	}
	
	
	results <- list(
		sample_names = sample_names,
		dirs = dirs,
		regions.list = regions.list,
		trinuc_bg_counts_ratio = ,
		trinuc_bg_counts.sigfit = ,
		trinuc_bg_ratio.sigfit = ,
		genome_trinuc_counts.sigfit = ,
		observed_corrected_trinuc_counts = ,
		observed_trinuc_counts.sigfit = ,
		mutation_burden = 
	)
	
	message("DONE")
	
	return(results)
}