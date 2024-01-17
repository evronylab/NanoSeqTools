#' Load NanoSeq data for region-specific analysis
#' 
#' Package requirements: tidyverse, rtracklayer, GenomicRanges, tabix binary
#'
#' This function loads NanoSeq data required for region-specific analysis.
#'
#' @param nanoseq_data Dataset resulting from load_nanoseq_data function
#' @param regions.list GRangesList object, comprised of GRanges that each contains a 'region set' to jointly analyze. The regions within each 'region set' can have overlaps (the functions handle this). Regions that were excluded when running load_nanoseq_data ('exclude_regions') are excluded from all region sets. If ignore.strand = FALSE, the strand of each region in the region set specifies which mutations to include: '+ and '-' strand include mutations where central pyrimidine is on the '+' and '-' strands of the reference genome, respectively, and '*' includes all mutations. If ignore.strand = TRUE, all mutations in the region are included regardless of strand. When there are overlapping regions with opposite strands within the same 'region set', the mutations in those overlapping regions are counted only once, because each mutation is a central pyrimidine on only one strand. Regions that are not in the contigs analyzed by the NanoSeq pipeline can be included in a 'region set', but they do not contribute any aspect of the data analysis. Best practice is to name the elements of regions.list, since these names are carried forward to the output.
#' @param ignore.strand TRUE or FALSE (default). Whether to ignore strand information in regions.list.
#' @param tabix_bin Full path of tabix binary
#' @return Returns NanoSeq results for each 'region set'. Samples without coverage of any of the regions are skipped and absent from the output.
#' * sample_names: A vector of all sample IDs that were loaded
#' * dir: A vector of the directories containing the NanoSeq results that were loaded
#' * regions.list: Copy of input regions.list
#' * excluded_samples: Names of samples excluded from the results because they do not have NanoSeq read coverage in any regions
#' * indel_counts.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of unique observed indel counts (created with indelwald indel.spectrum function), with one row per sample and one column per indel context. Indel counts do not take into account strand information regardless of region strand and the ignore.strand setting.
#' * trinuc_bg_counts_ratio: Data frame of the sample trinucleotide background counts (i.e. for each sample and region set, the number of interrogated bases for each trinucleotide context), the genome trinucleotide background counts (i.e. for each sample/region set, the number of each trinucleotide context in the genomic region set), and the normalized ratio of these for each sample/region set combination. Columns: sample, region, tri (trinucleotide context), sample_tri_bg, genome_tri_bg, ratio2genome.
#' * trinuc_bg_counts.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of the sample trinucleotide background counts, with one row per sample and one column per trinucleotide context.
#' * trinuc_bg_ratio.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of the ratio of the sample trinucleotide background counts (normalized to a sum of 1) to the genome trinucleotide background counts (normalized to a sum of 1), with one row per sample and one column per trinucleotide context.
#' * genome_trinuc_counts.sigfit: Vector of the genome trinucleotide background counts, in the same order as columns in sigfit format columns.
#' * observed_corrected_trinuc_counts: Data frame of observed and corrected mutation counts (for all mutations and for unique mutations) for each sample/region set combination. Columns: sample, region, tri (trinucleotide context), trint_subst_observed, trint_subst_unique_observed, ratio2genome, trint_subst_corrected, trint_subst_unique_corrected.
#' * observed_trinuc_counts.sigfit: List with one object per region set, each comprised of a data frame in sigfit format of unique observed mutation counts, with one row per sample and one column per trinucleotide substitution context.
#' * mutation_burden, with one row per sample/region combination. All substitution mutation statistics include all mutations, not just unique mutations. Indel statistics are calculated both for all and for unique mutations and do not take into account strand regardless of region strand and the ignore.strand setting.
#'  - The number of observed and corrected substitution mutations (muts_observed and muts_corrected)
#'  - Number of all and unique observed indels (indels_observed, indels_unique_observed)
#'  - Total number of observed and corrected interrogated bases (total_observed and total_corrected; note: observed and corrected are the same)
#'  - Observed and corrected substitution mutation burdens (burden_observed and burden_corrected)
#'  - All and unique observed indel mutation burden (burden_indels_observed and burden_indels_unique_observed)
#'  - Observed and corrected lower and upper confidence intervals of substitution mutation counts and all and unique observed lower and upper confidence intervals of indel counts (muts_lci_observed, muts_lci_corrected, indels_lci_observed, indels_unique_lci_observed, muts_uci_observed, muts_uci_corrected, indels_uci_observed, indels_unique_uci_observed)
#'  - Observed and corrected lower and upper confidence intervals of substitution mutation burdens and all and unique observed lower and upper confidence intervals of indel mutation burdens (burden_lci_observed, burden_lci_corrected, burden_indels_lci_observed, burden_indels_unique_lci_observed, burden_uci_observed, burden_uci_corrected, burden_indels_uci_observed, burden_indels_unique_uci_observed)
#' @export

load_nanoseq_regions <- function(nanoseq_data,regions.list,ignore.strand = FALSE, tabix_bin){

	#Load packages required only by this function
	suppressPackageStartupMessages(library(rtracklayer))
	suppressPackageStartupMessages(library(GenomicRanges))
	suppressPackageStartupMessages(library(nanoseq_data$BSgenomepackagename,character.only=TRUE))
	
	#Initialize lists for results
	excluded_samples <- c()
	trinuc_bg_counts_ratio <- list()
	vcf_snp.fix <- list()
	indel_counts <- list()
	mutation_burden <- list()
	
	dirs <- nanoseq_data$dirs
	sample_names <- nanoseq_data$sample_names
	
	#Filter out exclude_regions from region sets, and remove strand information if ignore.strand = TRUE
	regions.list <- map(regions.list,function(x) subtract(x,nanoseq_data$exclude_regions) %>% unlist) %>% GRangesList
	
	if(ignore.strand){
		regions.list <- map(regions.list,function(x) {strand(x) <- "*"; return(x)}) %>% GRangesList
	}
	
	#Convert BSgenome to StringSet object for indel spectrum loading
	BSgenome.StringSet <- as(sapply(seqnames(eval(parse(text=nanoseq_data$BSgenomepackagename))),function(x){eval(parse(text=nanoseq_data$BSgenomepackagename))[[x]]}),"DNAStringSet")
	
	message("Loading sample data...")
	pb <- txtProgressBar(min=0,max=100,style=3)
	
	for(i in 1:length(dirs)){
		
		dir <- dirs[i]
		sample_name <- sample_names[i]
		
		setTxtProgressBar(pb,i/length(dirs)*100)
		
		#Load NanoSeq bed coverage data for each 'region set'. First import coverage data for all regions combined for efficiency,
		# and then separate the results by region set.
		# Note, NanoSeq bed coverage data does not depend on strand information, since the coverage and trinucleotide background are the same regardless of strand.
		# Note, using tabix binary since it is much faster than R tabix functions.
		# Note, import function transforms bed coordinates to 1-based coordinates, which matched the mutation VCF 1-based coordinates.
		
		 #Load bed coverage information for all regions across all region sets.
		tmp.regions.all <- tempfile(fileext=".bed")
		regions.list %>% unlist %>% reduce(ignore.strand=TRUE) %>% sort %>% export(con=tmp.regions.all,format="bed")
		
		tmp.bedcov.all <- tempfile()
		system(paste(tabix_bin,paste0(dir,"/results.cov.bed.gz"),"-R",tmp.regions.all,"| sed 's/;/\t/g' | awk 'BEGIN{OFS=\"\t\"}{print $1,$2,$3,$6,$4,$5}' >",tmp.bedcov.all))
		
		bedcov.all <- import(tmp.bedcov.all,format="bedgraph")
		invisible(file.remove(tmp.regions.all,tmp.bedcov.all))
		
		#Skip samples without coverage in any regions
		if(length(bedcov.all)==0){
			message(paste("    ... Skipping sample",sample_name,"- no coverage in any regions."))
			excluded_samples <- c(excluded_samples,i)
			next
		}
		colnames(mcols(bedcov.all)) <- c("coverage","tri","ref")
		
		#Extract bed coverage information for each region set
		 # Note, suppressing warnings for situations when there is a chromosome in the region.list that is not in the coverage data and
		 # a chromosome in the coverage data that is not in the region.list.
		trinuc_bg_counts_ratio[[sample_name]] <- map(regions.list,function(x) suppressWarnings(subsetByOverlaps(bedcov.all,x)))
		
		#Calculate trinucleotide counts for each region set
		trinuc_bg_counts_ratio[[sample_name]] <- map(trinuc_bg_counts_ratio[[sample_name]],function(x){
			mcols(x)[,c("tri","coverage")] %>%
				as.data.frame %>%
				group_by(tri) %>%
				summarize(n=sum(coverage)) %>%
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
		vcf_snp.fix.gr <- nanoseq_data$vcf_snp.fix[[sample_name]] %>%
		    select(-c(ID,QUAL,FILTER)) %>%
		    mutate(tri=str_replace(INFO,".*TRI=(...>.).*","\\1"),
		           tri=str_c(str_sub(tri,1,4),str_sub(tri,1,1),str_sub(tri,5,5),str_sub(tri,3,3)),
		           strand=if_else(REF %in% c("C","T"),"+","-")) %>%
		    select(-c(REF,ALT,INFO)) %>%
		    makeGRangesFromDataFrame(seqnames.field="CHROM",start.field="POS",end.field="POS",keep.extra.columns=TRUE)
		vcf_indel.fix.gr <- nanoseq_data$vcf_indel.fix[[sample_name]] %>% vcf_indel_toGRanges
		
		 #Extract mutations in each region set, taking into account strand information for substitutions but not for indels.
		 # Note, trint_subst_corrected and trint_subst_unique_corrected can have values of NaN when both the numerator and denominator
		 #  values used to calculate them are both 0, and they can have values of Inf when only the denominator is 0.
		 # Note, suppressing warnings for situations when the mutations granges and the region set both have at least one
		 #  chromosome that is not in the other.
		vcf_snp.fix[[sample_name]] <- map(regions.list,function(x){
			left_join(
				data.frame(tri=trint_subs_labels),
				suppressWarnings(subsetByOverlaps(vcf_snp.fix.gr,x,ignore.strand=FALSE)) %>%
					as_tibble %>%
					dplyr::count(tri,name="trint_subst_observed"),
				by="tri"
			) %>%
			left_join(
	    	suppressWarnings(subsetByOverlaps(vcf_snp.fix.gr,x,ignore.strand=FALSE)) %>%
	    		as_tibble %>%
	    		distinct %>%
	    		dplyr::count(tri,name="trint_subst_unique_observed"),
	    	by="tri"
	    ) %>%
	    replace(is.na(.),0)
		})
		
		vcf_snp.fix[[sample_name]] <- map2(vcf_snp.fix[[sample_name]],trinuc_bg_counts_ratio[[sample_name]],function(x,y){
		  left_join(x %>% mutate(tri_short=str_sub(tri,1,3)),
		            y %>% dplyr::rename(tri_short=tri),
		  					by="tri_short") %>%
		    mutate(trint_subst_corrected = trint_subst_observed / ratio2genome,
		           trint_subst_unique_corrected = trint_subst_unique_observed / ratio2genome
		           )
		})
		
		vcf_indel.fix <- map(regions.list,function(x){
			nanoseq_data$vcf_indel.fix[[sample_name]] %>% filter(suppressWarnings(countOverlaps(vcf_indel.fix.gr,x,type="within",ignore.strand=TRUE)) > 0)
		})
		
		# Calculate number of unique indel counts for each indel context for each 'region set'
		indel_counts[[sample_name]] <- map(vcf_indel.fix,function(x){
			x %>%
				distinct(CHROM,POS,REF,ALT) %>%
				indel.spectrum(BSgenome.StringSet) %>%
				indelwald.to.sigfit
		})
		
		#Calculate mutation burdens for each 'region set'.
		# Note, corrected burdens may be NaN when there are 0 observed/corrected mutations in a region.
		mutation_burden[[sample_name]] <- pmap(
			list(x=vcf_snp.fix[[sample_name]],y=vcf_indel.fix,z=trinuc_bg_counts_ratio[[sample_name]]),
			function(x,y,z){
		  	data.frame(
		  		muts_observed = sum(x$trint_subst_observed),
		  		muts_corrected = sum(x$trint_subst_corrected,na.rm=TRUE),
		  		indels_observed = y %>% nrow,
		  		indels_unique_observed = y %>% distinct(CHROM,POS,REF,ALT) %>% nrow,
		  		total_observed = sum(z$sample_tri_bg),
		  		total_corrected = sum(z$sample_tri_bg)
	      ) %>%
					mutate(
		      	burden_observed = muts_observed / total_observed,
		      	burden_corrected = muts_corrected / total_corrected,
		      	burden_indels_observed = indels_observed / total_observed,
		      	burden_indels_unique_observed = indels_unique_observed / total_observed,
		      	muts_lci_observed = poisson.test(muts_observed)$conf.int[1],
		      	muts_lci_corrected = muts_lci_observed / muts_observed * muts_corrected,
		      	indels_lci_observed = poisson.test(indels_observed)$conf.int[1],
		      	indels_unique_lci_observed = poisson.test(indels_unique_observed)$conf.int[1],
		      	muts_uci_observed = poisson.test(muts_observed)$conf.int[2],
		      	muts_uci_corrected = muts_uci_observed / muts_observed * muts_corrected,
		      	indels_uci_observed = poisson.test(indels_observed)$conf.int[2],
		      	indels_unique_uci_observed = poisson.test(indels_unique_observed)$conf.int[2],
		      	burden_lci_observed = muts_lci_observed / total_observed,
		      	burden_lci_corrected = muts_lci_corrected / total_corrected,
		      	burden_indels_lci_observed = indels_lci_observed / total_observed,
		      	burden_indels_unique_lci_observed = indels_unique_lci_observed / total_observed,
		      	burden_uci_observed = muts_uci_observed / total_observed,
		      	burden_uci_corrected = muts_uci_corrected / total_corrected,
		      	burden_indels_uci_observed = indels_uci_observed / total_observed,
		      	burden_indels_unique_uci_observed = indels_unique_uci_observed / total_observed
	      	)
			})
		
	}
	close(pb)

	# Collapse lists to data frames
	message("Combining sample data into data frames...")
	
	indel_counts <- indel_counts %>%
		map(function(x) bind_rows(x,.id="region")) %>%
		bind_rows(.id="sample")
	
	trinuc_bg_counts_ratio <- trinuc_bg_counts_ratio %>%
	  map(function(x) bind_rows(x,.id="region")) %>%
	  bind_rows(.id="sample")
	
	observed_corrected_trinuc_counts <- vcf_snp.fix %>%
	  map(function(x) bind_rows(x,.id="region")) %>%
	  bind_rows(.id="sample") %>%
		dplyr::select(-c(tri_short,sample_tri_bg,genome_tri_bg)) %>%
		as_tibble
	
	mutation_burden <- mutation_burden %>%
	  map(function(x) bind_rows(x,.id="region")) %>%
	  bind_rows(.id="sample") %>%
		as_tibble
	
	#For each region set, create sigfit format data frames, with samples in rows and trinucleotide contexts in columns, for:
	# a) indel counts spectra
	# b) sample trinucleotide background counts
	# c) ratios of the sample trinucleotide background counts to the genome trinucleotide background counts
	# d) observed unique mutation counts. Note: using unique mutation counts, since that is a more faithful representation of the mutational process.
	indel_counts.sigfit <- indel_counts %>%
		split(.$region) %>%
		map(function(x){
			x %>%
				select(-region) %>%
				remove_rownames %>%
				column_to_rownames("sample")
			})
	
	trinuc_bg_counts.sigfit <- trinuc_bg_counts_ratio %>%
	  split(.$region) %>%
	  map(function(x){
	    x <- x %>%
	    	dplyr::select(sample,tri,sample_tri_bg) %>%
	      pivot_wider(names_from=tri,values_from=sample_tri_bg) %>%
	      column_to_rownames("sample")
	    x <- x[,genome_freqs_labels] %>% set_names(genome_freqs_labels)
	    return(x)
	  })
	
	trinuc_bg_ratio.sigfit <- trinuc_bg_counts_ratio %>%
	  split(.$region) %>%
	  map(function(x){
	    x <- x %>% dplyr::select(sample,tri,ratio2genome) %>%
	      pivot_wider(names_from=tri,values_from=ratio2genome) %>%
	      column_to_rownames("sample")
	    x <- x[,genome_freqs_labels] %>% set_names(genome_freqs_labels)
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
	
	#Remove excluded samples from output
	if(!is.null(excluded_samples)){
	  sample_names <- sample_names[-excluded_samples]
	  dirs <- dirs[-excluded_samples]
	}
	
	results <- list(
		sample_names = sample_names,
		dirs = dirs,
		regions.list = regions.list,
		excluded_samples = sample_names[excluded_samples],
		indel_counts.sigfit = indel_counts.sigfit,
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