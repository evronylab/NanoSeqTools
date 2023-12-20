#' Load NanoSeq data for region-specific analysis
#' 
#' Package requirements: tidyverse, rtracklayer, GenomicRanges, tabix binary
#'
#' For more information regarding NanoSeq output files, refer to the "output" section in the [NanoSeq GitHub](https://github.com/cancerit/NanoSeq).
#'
#' @param nanoseq_data Dataset resulting from load_nanoseq_data function
#' @param regions.list GRangesList object, comprised of GRanges that each contains a 'region set' to jointly analyze. The regions within each 'region set' can have overlaps (the functions handle this). The strand of each region in the region set specifies which mutations to include: '+ and '-' strand include mutations where central pyrimidine is on the '+' and '-' strands of the reference genome, respectively, and '*' includes all mutations. When there are overlapping regions with opposite strands within the same 'region set', the mutations in those overlapping regions are counted only once, because each mutation is a central pyrimidine on only one strand. Best practice is to name the elements of regions.list, since these names are carried forward to the output.
#' @param tabix_bin Full path to tabix binary
#' @return Returns NanoSeq results for each 'region set'.
#' * sample_names: A vector of all sample IDs that were loaded
#' * dir: A vector of the directories containing the NanoSeq results that were loaded
#' * regions.list: Copy of input regions.list
#' * trinuc_bg_counts_ratio: Data frame of the sample trinucleotide background counts (i.e. number of interrogated bases for each trinucleotide context), the genome trinucleotide background counts (i.e. number of each trinucleotide context), and the normalized ratio of these for each sample/region combination. Columns: sample, region, tri (trinucleotide context), sample_tri_bg, genome_tri_bg, ratio2genome.
#' * trinuc_bg_counts.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of the sample trinucleotide background counts, with one row per sample and one column per trinucleotide context.
#' * trinuc_bg_ratio.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of the ratio of the sample trinucleotide background counts (normalized to a sum of 1) to the genome trinucleotide background counts (normalized to a sum of 1), with one row per sample and one column per trinucleotide context.
#' * genome_trinuc_counts.sigfit: Vector of the genome trinucleotide background counts, in the same order as columns in sigfit format columns.
#' * observed_corrected_trinuc_counts: Data frame of observed and corrected mutation counts (for all mutations and for unique mutations) for each sample/region combination. Columns: sample, region, tri (trinucleotide context), trint_subst_observed, trint_subst_unique_observed, ratio2genome, trint_subst_corrected, trint_subst_unique_corrected.
#' * observed_trinuc_counts.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of unique observed mutation counts, with one row per sample and one column per trinucleotide substitution context.
#' * mutation_burden: Data frame of the total number of observed and corrected mutations, total number of observed and corrected interrogated bases (note: observed and corrected are the same), observed and corrected mutation burdens, observed and corrected lower and upper confidence intervals of mutation counts, and observed and corrected lower and upper confidence intervals of mutation burdens, for each sample/region combination. All these statistics include all mutations, not just unique mutations.
#' @export

load_nanoseq_regions <- function(nanoseq_data,regions.list,tabix_bin){

	#Load packages
	suppressPackageStartupMessages(library(tidyverse))
	suppressPackageStartupMessages(library(rtracklayer))
	suppressPackageStartupMessages(library(GenomicRanges))
	
	#Initialize lists for results
	trinuc_bg_counts_ratio <- list()
	vcf_snp.fix.gr <- list()
	mutation_burden <- list()
	
	dirs <- nanoseq_data$dirs
	
	message("Loading sample data...")
	for (i in 1:length(dirs)) {
		
		dir <- dirs[i]
		sample_name <- sample_names[i]
		
		message(paste0("  ",sample_name))
		
		#Load NanoSeq bed coverage data for each 'region set'. Begins by first importing coverage data for all regions together for efficiency,
		# and then separating them by region.
		# Note, NanoSeq bed coverage data does not depend on strand information, since the coverage and trinucleotide background are the same regardless of strand.
		# Note, using tabix binary since it is much faster than R tabix tools.
		# Note, import function transforms bed coordinates to 1-based coordinates.
		
		 #Load bed coverage information for all regions across all region sets.
		tmp.regions.all <- tempfile()
		regions.list %>% unlist %>% reduce %>% export(con=tmp.regions.all,format="bed")
		
		tmp.bedcov.all <- tempfile()
		system(paste(tabix_bin,paste0(dir,"/results.cov.bed.gz"),"-R",tmp.regions.all,"| sed 's/;/\t/g' | awk 'BEGIN{OFS=\"\t\"}{print $1,$2,$3,$6,$4,$5}'>",tmp.bedcov.all))
		
		bedcov.all <- import(tmp.bedcov.all,format="bedgraph")
		colnames(mcols(bedcov.all)) <- c("coverage","tri","ref")
		
		invisible(file.remove(tmp.regions.all,tmp.bedcov.all))
		
		 #Extract bed coverage information for each region set
		trinuc_bg_counts_ratio[[sample_name]] <- map(regions.list,function(x) subsetByOverlaps(bedcov.all,x,type="within"))
		
		#Calculate trinucleotide counts for each region set
		trinuc_bg_counts_ratio[[sample_name]] <- map(trinuc_bg_counts_ratio[[sample_name]],function(x){
			rep(x$tri,x$coverage) %>%
				as.data.frame %>%
				set_names("tri") %>%
				count(tri) %>%
				left_join(data.frame(tri=trinucleotides_64),.,by="tri") %>%
				replace(is.na(.),0) %>%
				deframe %>%
				trinucleotide64to32 %>%
				as_tibble(rownames="tri") %>%
				dplyr::rename(sample_tri_bg=value)
		})
		
		#Annotate with genome trinucleotide background counts and calculate normalized ratio of sample to genome trinucleotide counts
		genome_trinuc_counts <- nanoseq_data$genome_trinuc_counts.sigfit %>%
		  enframe %>%
		  distinct %>%
		  set_names(c("tri","genome_tri_bg"))
		
		trinuc_bg_counts_ratio[[sample_name]] <- map(trinuc_bg_counts_ratio[[sample_name]],function(x){
		  left_join(x,genome_trinuc_counts,by="tri") %>%
		    mutate(ratio2genome=(sample_tri_bg/sum(sample_tri_bg))/(genome_tri_bg/sum(genome_tri_bg)) )
		})
		
		#Extract mutation information for each region set
		 #Load mutations in each region set, and annotate strand of reference genome that is a central pyrimidine
		vcf_snp.fix.gr[[sample_name]] <- nanoseq_data$vcf_snp.fix[[sample_name]] %>%
		    select(-c(ID,QUAL,FILTER)) %>%
		    mutate(tri=str_replace(INFO,".*TRI=(...>.).*","\\1"),
		           tri=str_c(str_sub(tri,1,4),str_sub(tri,1,1),str_sub(tri,5,5),str_sub(tri,3,3)),
		           strand=if_else(REF %in% c("C","T"),"+","-")) %>%
		    select(-c(REF,ALT,INFO)) %>%
		    makeGRangesFromDataFrame(seqnames.field="CHROM",start.field="POS",end.field="POS",keep.extra.columns=TRUE)
		
		#Extract mutations in each region set, taking into account strand information.
		# Note, trint_subst_corrected and trint_subst_unique_corrected can have values of NaN when both the numerator and denominator
		# values used to calculate them are both 0, and they can have values of Inf when only the denominator is 0.
		vcf_snp.fix.gr[[sample_name]] <- map(regions.list,function(x){
		  subsetByOverlaps(vcf_snp.fix.gr[[sample_name]],x,type="within",ignore.strand=FALSE) %>%
		    as_tibble %>%
		    add_count(tri,name="trint_subst_observed") %>%
		    group_by(seqnames,start,tri) %>%
		    mutate(trint_subst_unique_observed=n_distinct(tri)) %>%
		    ungroup %>%
		    distinct %>%
		    dplyr::select(-c(seqnames,start,end,width,strand)) %>%
		    left_join(data.frame(tri=trint_subs_labels),.,by="tri") %>%
		    replace(is.na(.),0)
		})
		
		vcf_snp.fix.gr[[sample_name]] <- map2(vcf_snp.fix.gr[[sample_name]],trinuc_bg_counts_ratio[[sample_name]],function(x,y){
		  left_join(x %>% mutate(tri_short=str_sub(tri,1,3)),
		            y %>% dplyr::rename(tri_short=tri),by="tri_short") %>%
		    mutate(trint_subst_corrected = trint_subst_observed / ratio2genome,
		           trint_subst_unique_corrected = trint_subst_unique_observed / ratio2genome
		           )
		})
		
		#Calculate mutation burdens for each 'region set'.
		# Note, corrected burdens may be NaN when there are 0 observed/corrected mutations in a region.
		mutation_burden[[sample_name]] <- map(vcf_snp.fix.gr[[sample_name]],function(x){
		  data.frame(muts_observed = sum(x$trint_subst_observed,na.rm=TRUE),
		             muts_corrected = sum(x$trint_subst_corrected,na.rm=TRUE),
		             total_observed = sum(x$sample_tri_bg,na.rm=TRUE),
		             total_corrected = sum(x$sample_tri_bg,na.rm=TRUE)
      ) %>%
			mutate(
			  burden_observed = muts_observed / total_observed,
			  burden_corrected = muts_corrected / total_corrected,
			  muts_lci_observed = poisson.test(muts_observed)$conf.int[1],
			  muts_lci_corrected = muts_lci_observed / muts_observed * muts_corrected,
			  muts_uci_observed = poisson.test(muts_observed)$conf.int[2],
			  muts_uci_corrected = muts_uci_observed / muts_observed * muts_corrected,
			  burden_lci_observed = muts_lci_observed / total_observed,
			  burden_lci_corrected = muts_lci_corrected / total_corrected,
			  burden_uci_observed = muts_uci_observed / total_observed,
			  burden_uci_corrected = muts_uci_corrected / total_corrected
      )
		})
		
	}

	
	# Collapse lists to data frames
	message("Combining sample data into data frames...")
	trinuc_bg_counts_ratio <- trinuc_bg_counts_ratio %>%
	  map(function(x) bind_rows(x,.id="region")) %>%
	  bind_rows(.id="sample")
	
	observed_corrected_trinuc_counts <- vcf_snp.fix.gr %>%
	  map(function(x) bind_rows(x,.id="region")) %>%
	  bind_rows(.id="sample") %>%
		dplyr::select(-c(tri_short,sample_tri_bg,genome_tri_bg))
	
	mutation_burden <- mutation_burden %>%
	  map(function(x) bind_rows(x,.id="region")) %>%
	  bind_rows(.id="sample")
	
	#For each region set, create sigfit format data frames, with samples in rows and trinucleotide contexts in columns, for:
	# a) sample trinucleotide background counts
	# b) ratios of the sample trinucleotide background counts to the genome trinucleotide background counts
	# c) observed unique mutation counts. Note: using unique mutation counts, since that is a more faithful representation of the mutational process.
	trinuc_bg_counts.sigfit <- trinuc_bg_counts_ratio %>%
	  split(.$region) %>%
	  map(function(x){
	    x <- x %>% dplyr::select(-region) %>%
	      dplyr::select(sample,tri,sample_tri_bg) %>%
	      pivot_wider(names_from=tri,values_from=sample_tri_bg) %>%
	      column_to_rownames("sample")
	    x <- x[,genome_freqs_labels]
	    colnames(x) <- genome_freqs_labels
	    return(x)
	  })
	
	trinuc_bg_ratio.sigfit <- trinuc_bg_counts_ratio %>%
	  split(.$region) %>%
	  map(function(x){
	    x <- x %>% dplyr::select(sample,tri,ratio2genome) %>%
	      pivot_wider(names_from=tri,values_from=ratio2genome) %>%
	      column_to_rownames("sample")
	    x <- x[,genome_freqs_labels]
	    colnames(x) <- genome_freqs_labels
	    return(x)
	  })
	
	observed_trinuc_counts.sigfit <- observed_corrected_trinuc_counts %>%
	  split(.$region) %>%
	  map(function(x){
	    x <- x %>% dplyr::select(sample,tri,trint_subst_unique_observed) %>%
	      pivot_wider(names_from=tri,values_from=trint_subst_unique_observed) %>%
	      column_to_rownames("sample")
	    return(x)
	  })
	
	results <- list(
		sample_names = sample_names,
		dirs = dirs,
		regions.list = regions.list,
		trinuc_bg_counts_ratio = trinuc_bg_counts_ratio,
		trinuc_bg_counts.sigfit = trinuc_bg_counts.sigfit,
		trinuc_bg_ratio.sigfit = trinuc_bg_ratio.sigfit,
		genome_trinuc_counts.sigfit = nanoseq_data$genome_trinuc_counts.sigfit,
		observed_corrected_trinuc_counts = observed_corrected_trinuc_counts,
		observed_trinuc_counts.sigfit = observed_trinuc_counts.sigfit,
		mutation_burden = mutation_burden
	)
	
	message("DONE")
	
	return(results)
}