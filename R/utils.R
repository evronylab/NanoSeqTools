#' Shared functions for NanoSeqTools
#' @NoRD

#Load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))

# Order of trinucleotide context labels
trint_subs_labels <- c("ACA>AAA","ACC>AAC","ACG>AAG","ACT>AAT","CCA>CAA","CCC>CAC","CCG>CAG","CCT>CAT","GCA>GAA","GCC>GAC","GCG>GAG","GCT>GAT","TCA>TAA","TCC>TAC","TCG>TAG","TCT>TAT","ACA>AGA","ACC>AGC","ACG>AGG","ACT>AGT","CCA>CGA","CCC>CGC","CCG>CGG","CCT>CGT","GCA>GGA","GCC>GGC","GCG>GGG","GCT>GGT","TCA>TGA","TCC>TGC","TCG>TGG","TCT>TGT","ACA>ATA","ACC>ATC","ACG>ATG","ACT>ATT","CCA>CTA","CCC>CTC","CCG>CTG","CCT>CTT","GCA>GTA","GCC>GTC","GCG>GTG","GCT>GTT","TCA>TTA","TCC>TTC","TCG>TTG","TCT>TTT","ATA>AAA","ATC>AAC","ATG>AAG","ATT>AAT","CTA>CAA","CTC>CAC","CTG>CAG","CTT>CAT","GTA>GAA","GTC>GAC","GTG>GAG","GTT>GAT","TTA>TAA","TTC>TAC","TTG>TAG","TTT>TAT","ATA>ACA","ATC>ACC","ATG>ACG","ATT>ACT","CTA>CCA","CTC>CCC","CTG>CCG","CTT>CCT","GTA>GCA","GTC>GCC","GTG>GCG","GTT>GCT","TTA>TCA","TTC>TCC","TTG>TCG","TTT>TCT","ATA>AGA","ATC>AGC","ATG>AGG","ATT>AGT","CTA>CGA","CTC>CGC","CTG>CGG","CTT>CGT","GTA>GGA","GTC>GGC","GTG>GGG","GTT>GGT","TTA>TGA","TTC>TGC","TTG>TGG","TTT>TGT")

genome_freqs_labels <- str_sub(trint_subs_labels,1,3)

#Order of indel context labels
indel_labels <- c("1:Del:C:0","1:Del:C:1","1:Del:C:2","1:Del:C:3","1:Del:C:4","1:Del:C:5","1:Del:T:0","1:Del:T:1","1:Del:T:2","1:Del:T:3","1:Del:T:4","1:Del:T:5","1:Ins:C:0","1:Ins:C:1","1:Ins:C:2","1:Ins:C:3","1:Ins:C:4","1:Ins:C:5","1:Ins:T:0","1:Ins:T:1","1:Ins:T:2","1:Ins:T:3","1:Ins:T:4","1:Ins:T:5","2:Del:R:0","2:Del:R:1","2:Del:R:2","2:Del:R:3","2:Del:R:4","2:Del:R:5","3:Del:R:0","3:Del:R:1","3:Del:R:2","3:Del:R:3","3:Del:R:4","3:Del:R:5","4:Del:R:0","4:Del:R:1","4:Del:R:2","4:Del:R:3","4:Del:R:4","4:Del:R:5","5:Del:R:0","5:Del:R:1","5:Del:R:2","5:Del:R:3","5:Del:R:4","5:Del:R:5","2:Ins:R:0","2:Ins:R:1","2:Ins:R:2","2:Ins:R:3","2:Ins:R:4","2:Ins:R:5","3:Ins:R:0","3:Ins:R:1","3:Ins:R:2","3:Ins:R:3","3:Ins:R:4","3:Ins:R:5","4:Ins:R:0","4:Ins:R:1","4:Ins:R:2","4:Ins:R:3","4:Ins:R:4","4:Ins:R:5","5:Ins:R:0","5:Ins:R:1","5:Ins:R:2","5:Ins:R:3","5:Ins:R:4","5:Ins:R:5","2:Del:M:1","3:Del:M:1","3:Del:M:2","4:Del:M:1","4:Del:M:2","4:Del:M:3","5:Del:M:1","5:Del:M:2","5:Del:M:3","5:Del:M:4","5:Del:M:5")

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

#Function to convert vcf indel data frame to GRanges
# Note: in NanoSeq indel vcfs, always either REF or ALT has length = 1 and the other has length > 1.
vcf_indel_toGRanges <- function(x){
	x <- x %>%
		mutate(
			START=case_when(
				nchar(REF) > nchar(ALT) ~ POS + 1,
				.default = POS
			),
			END=case_when(
				nchar(REF) > nchar(ALT) ~ POS + nchar(REF) - 1,
				.default = POS
			)
		) %>%
		makeGRangesFromDataFrame(seqnames.field="CHROM",start.field="START",end.field="END")
	return(x)
}

# Internal readers and calculations for the DupCaller 1.2.11 output schema.
dupcaller_reference_trinuc <- function(chrom,pos,genome){
	if(any(!chrom %in% seqnames(genome))){
		stop("DupCaller contigs do not match the supplied BSgenome.")
	}
	result <- rep(NA_character_,length(pos))
	for(chr in unique(chrom)){
		idx <- which(chrom == chr & pos > 1 & pos < seqlengths(genome)[chr])
		if(length(idx) > 0){
			result[idx] <- as.character(Biostrings::Views(genome[[chr]],start=pos[idx] - 1,end=pos[idx] + 1))
		}
	}
	return(result)
}

dupcaller_version <- function(path){
	x <- grep("^Version:",readLines(path,warn=FALSE),value=TRUE)
	x <- trimws(sub("^Version:","",x))
	if(length(x) != 1 || !grepl("^[0-9]+[.][0-9]+[.][0-9]+$",x) || utils::compareVersion(x,"1.2.11") < 0){
		stop("DupCaller call and estimate outputs must be version 1.2.11 or newer: ",path)
	}
	return(x)
}

dupcaller_scalar <- function(path,key){
	x <- read_tsv(path,col_names=c("key","value"),col_types="cc",show_col_types=FALSE)
	value <- suppressWarnings(as.numeric(x$value[x$key == key]))
	if(length(value) != 1 || !is.finite(value) || value < 0){
		stop("Missing or invalid '",key,"' in ",path,". Run DupCaller estimate first.")
	}
	return(value)
}

dupcaller_indel_multiplier <- function(path,burden_path){
	x <- read_tsv(path,show_col_types=FALSE,name_repair="minimal")
	value <- x$indel83_number_genome
	if(nrow(x) != 83 || !is.numeric(value) || any(!is.finite(value) | value < 0)){
		stop("Expected DupCaller ID83 reference opportunities in ",path)
	}
	multiplier <- sum(value) / dupcaller_scalar(burden_path,"Reference base number")
	if(!is.finite(multiplier) || multiplier <= 0){
		stop("Invalid DupCaller reference indel opportunity multiplier.")
	}
	return(multiplier)
}

dupcaller_fold_sbs <- function(tri,alt,value,count_col){
	keep <- tri %in% trinucleotides_64 & alt != substr(tri,2,2)
	tri <- tri[keep]
	alt <- alt[keep]
	value <- value[keep]
	purine <- substr(tri,2,2) %in% c("A","G")
	tri[purine] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(tri[purine])))
	alt[purine] <- as.character(Biostrings::complement(Biostrings::DNAStringSet(alt[purine])))
	x <- tibble(tri=paste0(tri,">",substr(tri,1,1),alt,substr(tri,3,3)),value=value) %>%
		group_by(tri) %>% summarize(value=sum(value),.groups="drop")
	left_join(tibble(tri=trint_subs_labels),x,by="tri") %>%
		mutate(value=replace_na(value,0)) %>% dplyr::rename(!!count_col := value)
}

dupcaller_read_sbs_opportunities <- function(path){
	x <- read_tsv(path,show_col_types=FALSE,name_repair="minimal")
	names(x)[1] <- "label"
	expected <- expand_grid(tri=trinucleotides_64,alt=c("A","T","C","G")) %>%
		filter(alt != substr(tri,2,2)) %>% mutate(label=paste0(tri,">",alt))
	if(nrow(x) != 192 || anyDuplicated(x$label) || !setequal(x$label,expected$label) ||
		ncol(x) < 2 || !all(vapply(x[-1],is.numeric,logical(1)))){
		stop("Expected 192 alternate-specific DupCaller SBS opportunity rows in ",path)
	}
	values <- as.matrix(x[-1])
	if(any(!is.finite(values) | values < 0)) stop("Invalid DupCaller SBS opportunities in ",path)
	dupcaller_fold_sbs(substr(x$label,1,3),substr(x$label,5,5),rowSums(values),"sample_tri_bg")
}

dupcaller_read_coverage <- function(path,str_h5,indel_multiplier){
	if(length(indel_multiplier) != 1 || !is.finite(indel_multiplier) || indel_multiplier <= 0){
		stop("Missing DupCaller indel normalization; reload with load_nanoseq_data.")
	}
	columns <- c("CHROM","START","END",paste0("cov_",c("A","T","C","G")),
		"del_unit","ins_unit","del_2","del_3","del_4","del_5plus",
		"ins_A","ins_T","ins_C","ins_G","ins_len2","ins_len3","ins_len4","ins_len5plus","del_hp","ins_hp")
	x <- read_tsv(path,col_names=columns,col_types=paste0("c",strrep("d",22)),show_col_types=FALSE)
	if(ncol(x) != 23 || nrow(problems(x)) > 0 || anyNA(x) ||
		any(!is.finite(as.matrix(x[-1]))) || any(as.matrix(x[-1]) < 0) ||
		any(x$END != x$START + 1 | x$START != floor(x$START))){
		stop("Expected the 23-column, one-base DupCaller coverage BED: ",path)
	}
	# The raw BED channels have different multiplicities after expansion to ID83.
	# Deletion length bins contribute to STR plus 1/2/3/5 MH channels;
	# insertion length bins contribute to both repeat-count 0 and 1.
	weights <- c(1,1,2,3,4,6,1,1,1,1,2,2,2,2,1,1)
	x$indel_coverage <- as.vector(as.matrix(x[,8:23]) %*% weights)
	idx <- which(x$ins_unit > 0)
	if(length(idx) > 0){
		if(is.null(str_h5) || length(str_h5) != 1 || !file.exists(str_h5)){
			stop("Supply dupcaller_str_h5 (the matching reference .str.h5 file) for regional/excluded indel coverage.")
		}
		if(!requireNamespace("rhdf5",quietly=TRUE)) stop("Install Bioconductor rhdf5 to read dupcaller_str_h5.")
		for(chr in unique(x$CHROM[idx])){
			rows <- idx[x$CHROM[idx] == chr]
			positions <- x$START[rows] + 1
			# rhdf5 reverses h5py dimensions: Python row 1 becomes R column 2.
			repeats <- as.vector(rhdf5::h5read(str_h5,chr,index=list(positions,2)))
			if(length(repeats) != length(rows) || any(!is.finite(repeats) | repeats < 0)){
				stop("Invalid repeat annotations in dupcaller_str_h5 for ",chr)
			}
			x$indel_coverage[rows] <- x$indel_coverage[rows] + x$ins_unit[rows] * (repeats == 1)
		}
	}
	x$indel_coverage <- x$indel_coverage / indel_multiplier
	return(x)
}

dupcaller_coverage_sbs <- function(x,count_col){
	x <- x %>% dplyr::select(tri,all_of(paste0("cov_",c("A","T","C","G")))) %>%
		pivot_longer(-tri,names_to="alt",names_prefix="cov_",values_to="value")
	dupcaller_fold_sbs(x$tri,x$alt,x$value,count_col)
}

dupcaller_weights <- function(ratio){
	ifelse(is.finite(ratio) & ratio > 0,1 / ratio,0)
}

dupcaller_fay_feuer <- function(counts,weights){
	y <- sum(counts * weights)
	v <- sum(counts * weights^2)
	w_max <- max(weights)
	lower <- if(y > 0) qgamma(0.025,shape=y^2/v,scale=v/y) else 0
	y_star <- y + w_max
	v_star <- v + w_max^2
	upper <- if(y_star > 0) qgamma(0.975,shape=y_star^2/v_star,scale=v_star/y_star) else 0
	return(c(lower,upper))
}

dupcaller_corrected_ci <- function(burden,counts){
	ci <- dupcaller_fay_feuer(counts$trint_subst_observed,dupcaller_weights(counts$ratio2genome))
	burden$muts_lci_corrected <- ci[1]
	burden$muts_uci_corrected <- ci[2]
	burden$burden_lci_corrected <- ci[1] / burden$total_corrected
	burden$burden_uci_corrected <- ci[2] / burden$total_corrected
	return(burden)
}

dupcaller_dbs_ranges <- function(x){
	if(any(nchar(x$REF) != 2 | nchar(x$ALT) != 2)) stop("Expected two-base DupCaller DBS records.")
	x %>% mutate(END=POS + 1) %>%
		makeGRangesFromDataFrame(seqnames.field="CHROM",start.field="POS",end.field="END")
}

dupcaller_check_remaining <- function(x){
	if(any(!is.finite(x)) || any(x < -1e-6)){
		stop("Excluded coverage exceeds total coverage; check that DupCaller outputs and reference match.")
	}
	return(pmax(x,0))
}

#Function to annotate PASS SBS records from DupCaller 1.2.11
annotate_dupcaller_sbs_records <- function(x,BSgenome_object){
	if(any(!grepl("^[ACGT]$",x$REF) | !grepl("^[ACGT]$",x$ALT) | x$REF == x$ALT)){
		stop("Expected single-base PASS records in DupCaller SBS VCF; DBS must be loaded separately.")
	}
	x$TN <- dupcaller_reference_trinuc(x$CHROM,x$POS,BSgenome_object)
	return(x)
}

#Function to keep DupCaller SNV records that can be represented in SBS96 space
dupcaller_sbs_variants <- function(x){
	x <- x %>%
		mutate(tri=NA_character_)
	ref <- toupper(as.character(x$REF))
	alt <- toupper(as.character(x$ALT))
	tri <- toupper(as.character(x$TN))
	tri_sbs <- rep(NA_character_,nrow(x))
	valid <- !is.na(ref) & !is.na(alt) & !is.na(tri) &
		nchar(ref) == 1 & nchar(alt) == 1 & nchar(tri) == 3 &
		ref %in% c("A","C","G","T") & alt %in% c("A","C","G","T") &
		tri %in% trinucleotides_64
	if(any(valid)){
		tri_norm <- tri
		alt_norm <- alt
		purine_ref <- valid & ref %in% c("A","G")
		if(any(purine_ref)){
			tri_norm[purine_ref] <- as.character(reverseComplement(DNAStringSet(tri[purine_ref])))
			alt_norm[purine_ref] <- as.character(complement(DNAStringSet(alt[purine_ref])))
		}
		tri_sbs[valid] <- str_c(tri_norm[valid],">",str_sub(tri_norm[valid],1,1),alt_norm[valid],str_sub(tri_norm[valid],3,3))
	}
	x <- x %>%
		mutate(tri=tri_sbs) %>%
		mutate(strand=if_else(REF %in% c("C","T"),"+","-")) %>%
		filter(!is.na(tri))
	return(x)
}

#Function to annotate DupCaller per-position coverage with reference trinucleotide context
annotate_dupcaller_coverage_trinuc <- function(x,BSgenome_object){
	x$tri <- dupcaller_reference_trinuc(x$CHROM,x$START + 1,BSgenome_object)
	return(x)
}

#Function to import indels for spectrum analysis, modified from INDELWALD package:
## Max Stammnitz; maxrupsta@gmail.com; University of Cambridge  ##
## Citation: The evolution of two transmissible cancers in Tasmanian devils (Stammnitz et al. 2023, Science 380:6642)
indel.spectrum <- function(x,reference){
  
  ## 1. split VCF indels into:
  # (i) 1-bp deletions
  dels.1bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-1,,drop=F]
  
  # (ii) 1-bp insertions
  ins.1bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-1,,drop=F]
  
  # (iii) >1 bp deletions
  
  ## 2
  dels.2bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-2,,drop=F]
  
  ## 3
  dels.3bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-3,,drop=F]
  
  ## 4
  dels.4bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-4,,drop=F]
  
  ## 5+
  dels.5bp <- x[nchar(as.character(x[,'ALT'])) <= nchar(as.character(x[,'REF']))-5,,drop=F]
  
  # (iv) >1 bp insertions
  
  ## 2
  ins.2bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-2,,drop=F]
  
  ## 3
  ins.3bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-3,,drop=F]
  
  ## 4
  ins.4bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-4,,drop=F]
  
  ## 5+
  ins.5bp <- x[nchar(as.character(x[,'REF'])) <= nchar(as.character(x[,'ALT']))-5,,drop=F]
  
  ## 2. classify 1 bp events into:
  
  # (i) 1 bp deletions at homopolymers (length 1 == "no neighbouring homopolymers")
  if(nrow(dels.1bp) > 0){
    
    ## extract 10 bp upstream/downstream sequence context from reference
    dels.1bp.context <- as.character(subseq(x = reference[as.character(dels.1bp[,'CHROM'])], 
                                            start = as.numeric(dels.1bp[,'POS']) - 9, 
                                            end = as.numeric(dels.1bp[,'POS']) + 11))
    dels.1bp.context.middle <- paste0('[', str_split_fixed(dels.1bp.context, '', 21)[,11,drop=F], ']')
    dels.1bp.context.start <- str_split_fixed(dels.1bp.context, '', 11)[,1:10,drop=F]
    dels.1bp.context.start <- paste(dels.1bp.context.start[,1], dels.1bp.context.start[,2], dels.1bp.context.start[,3],
                                    dels.1bp.context.start[,4], dels.1bp.context.start[,5], dels.1bp.context.start[,6],
                                    dels.1bp.context.start[,7], dels.1bp.context.start[,8], dels.1bp.context.start[,9],
                                    dels.1bp.context.start[,10], sep = '')
    dels.1bp.context.end <- str_split_fixed(dels.1bp.context, '', 12)[,12,drop=F]
    dels.1bp[,'TRIPLET'] <- paste(dels.1bp.context.start, dels.1bp.context.middle, dels.1bp.context.end, sep = '')
    colnames(dels.1bp)[5] <- 'CONTEXT FW' 
    
    ## need the central base to be pyrimidine-centred, i.e. C or T, hence also run reverseComplements on full contexts
    dels.1bp.context.rc <- as.character(reverseComplement(subseq(x = reference[as.character(dels.1bp[,'CHROM'])], 
                                                                 start = as.numeric(dels.1bp[,'POS']) - 9, 
                                                                 end = as.numeric(dels.1bp[,'POS']) + 11)))
    dels.1bp.context.rc.middle <- paste0('[', str_split_fixed(dels.1bp.context.rc, '', 21)[,11,drop=F], ']')
    dels.1bp.context.rc.start <- str_split_fixed(dels.1bp.context.rc, '', 11)[,1:10,drop=F]
    dels.1bp.context.rc.start <- paste(dels.1bp.context.rc.start[,1], dels.1bp.context.rc.start[,2], dels.1bp.context.rc.start[,3],
                                       dels.1bp.context.rc.start[,4], dels.1bp.context.rc.start[,5], dels.1bp.context.rc.start[,6],
                                       dels.1bp.context.rc.start[,7], dels.1bp.context.rc.start[,8], dels.1bp.context.rc.start[,9],
                                       dels.1bp.context.rc.start[,10], sep = '')
    dels.1bp.context.rc.end <- str_split_fixed(dels.1bp.context.rc, '', 12)[,12,drop=F]
    dels.1bp <- cbind(dels.1bp[,1:5,drop=F], paste(dels.1bp.context.rc.start, dels.1bp.context.rc.middle, dels.1bp.context.rc.end, sep = ''))
    colnames(dels.1bp)[6] <- 'CONTEXT RC'
    
    ## summarise 1 bp deletions in matrix format
    dels.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(dels.1bp.summary) <- c('C', 'T') ## pyrimidine-centred deleted base
    rownames(dels.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    fw <- str_split_fixed(str_split_fixed(dels.1bp[,'CONTEXT FW'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    rc <- str_split_fixed(str_split_fixed(dels.1bp[,'CONTEXT RC'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    fw.if <- fw == 'C' | fw == 'T'
    for (i in 1:nrow(dels.1bp)){
      
      if(fw.if[i] == T){
        
        ## look at forward context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT FW']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == fw[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT FW']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == fw[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer (upstream vs. downstream context) is longer
        dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] <- dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] + 1
        
      } else {
        
        ## look at reverse complement context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT RC']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == rc[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT RC']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == rc[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer (upstream vs. downstream context) is longer
        dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] <- dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] + 1
      }
      
    }
    
  }else{
    
    ## summarise 1 bp deletions in matrix format
    dels.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(dels.1bp.summary) <- c('C', 'T') ## pyrimidine-centred deleted base
    rownames(dels.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    
  }
  
  # (ii) 1 bp insertions at homopolymers (length 0 == "no neighbouring homopolymers")
  if(nrow(ins.1bp) > 0){
    
    ## extract 10 bp upstream/downstream sequence context from reference
    ins.1bp.context <- as.character(subseq(x = reference[as.character(ins.1bp[,'CHROM'])], 
                                           start = as.numeric(ins.1bp[,'POS']) - 9, 
                                           end = as.numeric(ins.1bp[,'POS']) + 10))
    ## insertion after base pos. 10
    ins.1bp.context.middle <- paste0('[', str_split_fixed(as.character(ins.1bp[,'ALT']), '', 2)[,2,drop=F], ']')
    ins.1bp.context.start <- str_split_fixed(ins.1bp.context, '', 11)[,1:10,drop=F]
    ins.1bp.context.start <- paste(ins.1bp.context.start[,1], ins.1bp.context.start[,2], ins.1bp.context.start[,3],
                                   ins.1bp.context.start[,4], ins.1bp.context.start[,5], ins.1bp.context.start[,6],
                                   ins.1bp.context.start[,7], ins.1bp.context.start[,8], ins.1bp.context.start[,9],
                                   ins.1bp.context.start[,10], sep = '')
    ins.1bp.context.end <- str_split_fixed(ins.1bp.context, '', 11)[,11,drop=F]
    ins.1bp[,'TRIPLET'] <- paste(ins.1bp.context.start, ins.1bp.context.middle, ins.1bp.context.end, sep = '')
    colnames(ins.1bp)[5] <- 'CONTEXT FW'
    
    ## need the central base to be pyrimidine-centred, i.e. C or T, hence also run reverseComplements on full contexts
    ins.1bp.context.rc <- as.character(reverseComplement(subseq(x = reference[as.character(ins.1bp[,'CHROM'])], 
                                                                start = as.numeric(ins.1bp[,'POS']) - 9, 
                                                                end = as.numeric(ins.1bp[,'POS']) + 10)))
    ins.1bp.context.rc.middle <- paste0('[', as.character(reverseComplement(DNAStringSet(str_split_fixed(as.character(ins.1bp[,'ALT']), '', 2)[,2,drop=F]))), ']')
    ins.1bp.context.rc.start <- str_split_fixed(ins.1bp.context.rc, '', 11)[,1:10,drop=F]
    ins.1bp.context.rc.start <- paste(ins.1bp.context.rc.start[,1], ins.1bp.context.rc.start[,2], ins.1bp.context.rc.start[,3],
                                      ins.1bp.context.rc.start[,4], ins.1bp.context.rc.start[,5], ins.1bp.context.rc.start[,6],
                                      ins.1bp.context.rc.start[,7], ins.1bp.context.rc.start[,8], ins.1bp.context.rc.start[,9],
                                      ins.1bp.context.rc.start[,10], sep = '')
    ins.1bp.context.rc.end <- str_split_fixed(ins.1bp.context.rc, '', 11)[,11,drop=F]
    ins.1bp <- cbind(ins.1bp[,1:5,drop=F], paste(ins.1bp.context.rc.start, ins.1bp.context.rc.middle, ins.1bp.context.rc.end, sep = ''))
    colnames(ins.1bp)[6] <- 'CONTEXT RC'
    
    ## summarise 1 bp insertions in matrix format
    ins.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(ins.1bp.summary) <- c('C', 'T') ## pyrimidine-centred inserted base
    rownames(ins.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    fw <- str_split_fixed(str_split_fixed(ins.1bp[,'CONTEXT FW'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    rc <- str_split_fixed(str_split_fixed(ins.1bp[,'CONTEXT RC'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    fw.if <- fw == 'C' | fw == 'T'
    for (i in 1:nrow(ins.1bp)){
      
      if(fw.if[i] == T){
        
        ## look at forward context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT FW']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == fw[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT FW']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == fw[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer is longer
        ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] <- ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] + 1
        
      } else {
        
        ## look at reverse complement context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT RC']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == rc[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT RC']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == rc[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer is longer
        ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] <- ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] + 1
      }
      
    } 
    
  }else{
    
    ## summarise 1 bp insertions in matrix format
    ins.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(ins.1bp.summary) <- c('C', 'T') ## pyrimidine-centred inserted base
    rownames(ins.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    
  }
  
  ## 3. classify >=2 bp deletions into:
  
  # (i) 2 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.2bp) > 0){
    
    ## extract 5 x 2 bp upstream/downstream sequence context from reference
    dels.2bp.context <- as.character(subseq(x = reference[as.character(dels.2bp[,'CHROM'])], 
                                            start = as.numeric(dels.2bp[,'POS']) - 9, 
                                            end = as.numeric(dels.2bp[,'POS']) + 2 + 10))
    dels.2bp.context.middle <- str_split_fixed(dels.2bp.context, '', 22)[,11:12,drop=F]
    dels.2bp.context.middle <- paste0('[', dels.2bp.context.middle[,1], dels.2bp.context.middle[,2], ']')
    dels.2bp.context.start <- str_split_fixed(dels.2bp.context, '', 11)[,1:10,drop=F]
    dels.2bp.context.start <- paste(dels.2bp.context.start[,1], dels.2bp.context.start[,2], dels.2bp.context.start[,3],
                                    dels.2bp.context.start[,4], dels.2bp.context.start[,5], dels.2bp.context.start[,6],
                                    dels.2bp.context.start[,7], dels.2bp.context.start[,8], dels.2bp.context.start[,9],
                                    dels.2bp.context.start[,10], sep = '')
    dels.2bp.context.end <- str_split_fixed(dels.2bp.context, '', 13)[,13,drop=F]
    dels.2bp[,'TRIPLET'] <- paste(dels.2bp.context.start, dels.2bp.context.middle, dels.2bp.context.end, sep = '')
    colnames(dels.2bp)[5] <- 'CONTEXT' 
    
  }
  
  # (ii) 3 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.3bp) > 0){
    
    ## extract 5 x 3 bp upstream/downstream sequence context from reference
    dels.3bp.context <- as.character(subseq(x = reference[as.character(dels.3bp[,'CHROM'])], 
                                            start = as.numeric(dels.3bp[,'POS']) - 14, 
                                            end = as.numeric(dels.3bp[,'POS']) + 3 + 15))
    dels.3bp.context.middle <- str_split_fixed(dels.3bp.context, '', 33)[,16:18,drop=F]
    dels.3bp.context.middle <- paste0('[', dels.3bp.context.middle[,1], dels.3bp.context.middle[,2], dels.3bp.context.middle[,3], ']')
    dels.3bp.context.start <- str_split_fixed(dels.3bp.context, '', 33)[,1:15,drop=F]
    dels.3bp.context.start <- paste(dels.3bp.context.start[,1], dels.3bp.context.start[,2], dels.3bp.context.start[,3],
                                    dels.3bp.context.start[,4], dels.3bp.context.start[,5], dels.3bp.context.start[,6],
                                    dels.3bp.context.start[,7], dels.3bp.context.start[,8], dels.3bp.context.start[,9],
                                    dels.3bp.context.start[,10], dels.3bp.context.start[,11], dels.3bp.context.start[,12], 
                                    dels.3bp.context.start[,13], dels.3bp.context.start[,14], dels.3bp.context.start[,15], sep = '')
    dels.3bp.context.end <- str_split_fixed(dels.3bp.context, '', 19)[,19,drop=F]
    dels.3bp[,'TRIPLET'] <- paste(dels.3bp.context.start, dels.3bp.context.middle, dels.3bp.context.end, sep = '')
    colnames(dels.3bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iii) 4 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.4bp) > 0){
    
    ## extract 5 x 4 bp upstream/downstream sequence context from reference
    dels.4bp.context <- as.character(subseq(x = reference[as.character(dels.4bp[,'CHROM'])], 
                                            start = as.numeric(dels.4bp[,'POS']) - 19, 
                                            end = as.numeric(dels.4bp[,'POS']) + 4 + 20))
    dels.4bp.context.middle <- str_split_fixed(dels.4bp.context, '', 41)[,21:24,drop=F]
    dels.4bp.context.middle <- paste0('[', dels.4bp.context.middle[,1], 
                                      dels.4bp.context.middle[,2], 
                                      dels.4bp.context.middle[,3], 
                                      dels.4bp.context.middle[,4], ']')
    dels.4bp.context.start <- str_split_fixed(dels.4bp.context, '', 21)[,1:20,drop=F]
    dels.4bp.context.start <- paste(dels.4bp.context.start[,1], dels.4bp.context.start[,2], dels.4bp.context.start[,3],
                                    dels.4bp.context.start[,4], dels.4bp.context.start[,5], dels.4bp.context.start[,6],
                                    dels.4bp.context.start[,7], dels.4bp.context.start[,8], dels.4bp.context.start[,9],
                                    dels.4bp.context.start[,10], dels.4bp.context.start[,11], dels.4bp.context.start[,12], 
                                    dels.4bp.context.start[,13], dels.4bp.context.start[,14], dels.4bp.context.start[,15], 
                                    dels.4bp.context.start[,16], dels.4bp.context.start[,17], dels.4bp.context.start[,18], 
                                    dels.4bp.context.start[,19], dels.4bp.context.start[,20], sep = '')
    dels.4bp.context.end <- str_split_fixed(dels.4bp.context, '', 25)[,25,drop=F]
    dels.4bp[,'TRIPLET'] <- paste(dels.4bp.context.start, dels.4bp.context.middle, dels.4bp.context.end, sep = '')
    colnames(dels.4bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iv) 5+ bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.5bp) > 0){
    
    ## extract 5 x 100 bp upstream/downstream sequence context from reference
    ## i.e. max. 100 bp repeat motif
    dels.5bp.context <- as.character(subseq(x = reference[as.character(dels.5bp[,'CHROM'])], 
                                            start = as.numeric(dels.5bp[,'POS']) - 499, 
                                            end = as.numeric(dels.5bp[,'POS']) + 5 + 500))
    
    ### account for the different repeat lengths: iterate
    dels.5bp.context.middle <- as.character(dels.5bp[,'REF'])
    dels.5bp.context.middle <- str_split_fixed(dels.5bp.context.middle, '', 2)[,2,drop=F]
    dels.5bp.context.middle <- paste0('[', dels.5bp.context.middle, ']')
    dels.5bp.context.middle.lengths <- nchar(as.character(dels.5bp[,'REF'])) - 1
    dels.5bp.context.start <- rep(NA, nrow(dels.5bp))
    for (i in 1:length(dels.5bp.context.start)){
      tmp.dels.5bp.context.start <- str_split_fixed(dels.5bp.context[i], '', 501)[,c(1 + 500-c(5*dels.5bp.context.middle.lengths[i])):500,drop=F]
      dels.5bp.context.start[i] <- paste(tmp.dels.5bp.context.start, collapse = '')
    }
    dels.5bp.context.end <- rep(NA, nrow(dels.5bp))
    for (i in 1:length(dels.5bp.context.end)){
      tmp.dels.5bp.context.end <- str_split_fixed(dels.5bp.context[i], '', 1006)[,c(501+dels.5bp.context.middle.lengths[i]):c(501 + dels.5bp.context.middle.lengths[i]*6 - 1),drop=F]
      dels.5bp.context.end[i] <- paste(tmp.dels.5bp.context.end, collapse = '')
    }
    dels.5bp[,'TRIPLET'] <- paste(dels.5bp.context.start, dels.5bp.context.middle, dels.5bp.context.end, sep = '')
    colnames(dels.5bp)[5] <- 'CONTEXT'
    
  }
  
  ## summarise >=2 bp deletions at simple repeats in matrix format
  ## also create sub-tables for repeat length = 1 hits, to later assess these for microhomologies
  dels.greater.2bp.summary <- matrix(0, ncol = 4, nrow = 6)
  colnames(dels.greater.2bp.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## deletion size
  rownames(dels.greater.2bp.summary) <- c('1 or MH', '2', '3', '4', '5', '6+') ## number of repeats
  
  ## 2 bp
  if(nrow(dels.2bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.2bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.2bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '2 bp'] <- dels.greater.2bp.summary['6+', '2 bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '2 bp'] <- dels.greater.2bp.summary['5', '2 bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '2 bp'] <- dels.greater.2bp.summary['4', '2 bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '2 bp'] <- dels.greater.2bp.summary['3', '2 bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '2 bp'] <- dels.greater.2bp.summary['2', '2 bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '2 bp'] <- dels.greater.2bp.summary['1 or MH', '2 bp'] + 1
        dels.2bp.pot.MH.ind <- c(dels.2bp.pot.MH.ind, i)
        
      }
      
    }
    dels.2bp.pot.MH <- dels.2bp[dels.2bp.pot.MH.ind,,drop=F] 
    
  }
  
  ## 3 bp
  if(nrow(dels.3bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.3bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.3bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '3 bp'] <- dels.greater.2bp.summary['6+', '3 bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '3 bp'] <- dels.greater.2bp.summary['5', '3 bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '3 bp'] <- dels.greater.2bp.summary['4', '3 bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '3 bp'] <- dels.greater.2bp.summary['3', '3 bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '3 bp'] <- dels.greater.2bp.summary['2', '3 bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] + 1
        dels.3bp.pot.MH.ind <- c(dels.3bp.pot.MH.ind, i)
        
      }
      
    }
    dels.3bp.pot.MH <- dels.3bp[dels.3bp.pot.MH.ind,,drop=F] 
    
  }
  
  ## 4 bp
  if(nrow(dels.4bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.4bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.4bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '4 bp'] <- dels.greater.2bp.summary['6+', '4 bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '4 bp'] <- dels.greater.2bp.summary['5', '4 bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '4 bp'] <- dels.greater.2bp.summary['4', '4 bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '4 bp'] <- dels.greater.2bp.summary['3', '4 bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '4 bp'] <- dels.greater.2bp.summary['2', '4 bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] + 1
        dels.4bp.pot.MH.ind <- c(dels.4bp.pot.MH.ind, i)
        
      }
      
    }
    dels.4bp.pot.MH <- dels.4bp[dels.4bp.pot.MH.ind,,drop=F]
    
  }
  
  ## 5+ bp
  if(nrow(dels.5bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.5bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.5bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '5+ bp'] <- dels.greater.2bp.summary['6+', '5+ bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '5+ bp'] <- dels.greater.2bp.summary['5', '5+ bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '5+ bp'] <- dels.greater.2bp.summary['4', '5+ bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '5+ bp'] <- dels.greater.2bp.summary['3', '5+ bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '5+ bp'] <- dels.greater.2bp.summary['2', '5+ bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] + 1
        dels.5bp.pot.MH.ind <- c(dels.5bp.pot.MH.ind, i)
        
      }
      
    }
    dels.5bp.pot.MH <- dels.5bp[dels.5bp.pot.MH.ind,,drop=F] 
    
  }
  
  # (v) >=2 bp deletions at microhomologies
  ## go directly into in matrix format, thereby also re-classify repeat length = 1 deletions
  dels.greater.2bp.MH.summary <- matrix(0, ncol = 4, nrow = 5)
  colnames(dels.greater.2bp.MH.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## deletion size
  rownames(dels.greater.2bp.MH.summary) <- c('1 bp MH', '2 bp MH', '3 bp MH', '4 bp MH', '5+ bp MH') ## microhomology length
  dels.greater.2bp.MH.summary[2:5, '2 bp'] <- NA
  dels.greater.2bp.MH.summary[3:5, '3 bp'] <- NA
  dels.greater.2bp.MH.summary[4:5, '4 bp'] <- NA
  
  ## 2 bp
  if(nrow(dels.2bp) > 0){
    
    if(nrow(dels.2bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.2bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        
        ## isolate first nucleotide in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first nucleotide in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1]
        
        if(strsplit(tmp.repeat.nts, '')[[1]][2] == tmp.upstream.context | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context){
          
          dels.greater.2bp.MH.summary['1 bp MH', '2 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '2 bp'] + 1
          dels.greater.2bp.summary['1 or MH', '2 bp'] <- dels.greater.2bp.summary['1 or MH', '2 bp'] - 1
          
        }
      } 
      
    }
    
  }
  
  ## 3 bp
  if(nrow(dels.3bp) > 0){
    
    if(nrow(dels.3bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.3bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        
        ## isolate first two nucleotides in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 1):length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first two nucleotides in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:2]
        
        ### MH length = 2 (important to start with the highest possible MH length)
        
        if(all(c(strsplit(tmp.repeat.nts, '')[[1]][2:3] == tmp.upstream.context) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context) == T)){
          
          dels.greater.2bp.MH.summary['2 bp MH', '3 bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '3 bp'] + 1
          dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] - 1
          
        }else{
          
          ### MH length = 1
          
          if(strsplit(tmp.repeat.nts, '')[[1]][3] == tmp.upstream.context[2] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
            
            dels.greater.2bp.MH.summary['1 bp MH', '3 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '3 bp'] + 1
            dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] - 1
            
          }
          
        }
        
      } 
      
    }
    
  }
  
  ## 4 bp
  if(nrow(dels.4bp) > 0){
    
    if(nrow(dels.4bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.4bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        
        ## isolate first three nucleotides in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 2):length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first three nucleotides in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:3]
        
        ### MH length = 3 (important to start with the highest possible MH length)
        
        if(all(c(strsplit(tmp.repeat.nts, '')[[1]][2:4] == tmp.upstream.context) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:3] == tmp.downstream.context) == T)){
          
          dels.greater.2bp.MH.summary['3 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['3 bp MH', '4 bp'] + 1
          dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
          
        }else{
          
          ### MH length = 2
          
          if(all(c(strsplit(tmp.repeat.nts, '')[[1]][3:4] == tmp.upstream.context[2:3]) ==T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context[1:2]) == T)){
            
            dels.greater.2bp.MH.summary['2 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '4 bp'] + 1
            dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
            
          }else{
            
            ### MH length = 1
            
            if(strsplit(tmp.repeat.nts, '')[[1]][4] == tmp.upstream.context[3] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
              
              dels.greater.2bp.MH.summary['1 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '4 bp'] + 1
              dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
              
            }
            
          }
          
        }
        
      } 
      
    }
    
  }
  
  ## 5+ bp
  if(nrow(dels.5bp) > 0){
    
    if(nrow(dels.5bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.5bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        tmp.repeat.length <- nchar(tmp.repeat.nts)
        
        ## isolate first X (= repeat length - 1) nucleotides in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 4):length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first X (= repeat length - 1) nucleotides in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:5]
        
        ### MH length = 5+ (important to start with the highest possible MH length); only need to look at first 5 up/downstream NTs
        
        if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-4):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 4):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:5] == tmp.downstream.context[1:5]) == T)){
          
          dels.greater.2bp.MH.summary['5+ bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['5+ bp MH', '5+ bp'] + 1
          dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
          
        }else{
          
          ### MH length = 4
          
          if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-3):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 3):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:4] == tmp.downstream.context[1:4]) == T)){
            
            dels.greater.2bp.MH.summary['4 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['4 bp MH', '5+ bp'] + 1
            dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
            
          }else{
            
            ### MH length = 3
            
            if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-2):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 2):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:3] == tmp.downstream.context[1:3]) == T)){
              
              dels.greater.2bp.MH.summary['3 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['3 bp MH', '5+ bp'] + 1
              dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
              
            }else{
              
              ### MH length = 2
              
              if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-1):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 1):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context[1:2]) == T)){
                
                dels.greater.2bp.MH.summary['2 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '5+ bp'] + 1
                dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
                
              }else{
                
                if(strsplit(tmp.repeat.nts, '')[[1]][tmp.repeat.length] == tmp.upstream.context[length(tmp.upstream.context)] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
                  
                  dels.greater.2bp.MH.summary['1 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '5+ bp'] + 1
                  dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
                  
                }
                
              }
              
            }
            
          }
          
        }
        
      } 
      
    }
    
  }
  
  ## after MH cases have been taken out, rename row of deletions at 1-unit repeats in original table
  rownames(dels.greater.2bp.summary)[1] <- '1'
  
  ## 4. classify >=2 bp insertions into:
  
  # (i) 2 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.2bp) > 0){
    
    ins.2bp.context <- as.character(subseq(x = reference[as.character(ins.2bp[,'CHROM'])], 
                                           start = as.numeric(ins.2bp[,'POS']) - 9, 
                                           end = as.numeric(ins.2bp[,'POS']) + 10))
    ins.2bp.context.middle <- as.character(ins.2bp[,'ALT'])
    ins.2bp.context.middle <- paste0('[', str_split_fixed(ins.2bp.context.middle, '', 2)[,2,drop=F], ']')
    ins.2bp.context.start <- str_split_fixed(ins.2bp.context, '', 11)[,1:10,drop=F]
    ins.2bp.context.start <- paste(ins.2bp.context.start[,1], ins.2bp.context.start[,2], ins.2bp.context.start[,3],
                                   ins.2bp.context.start[,4], ins.2bp.context.start[,5], ins.2bp.context.start[,6],
                                   ins.2bp.context.start[,7], ins.2bp.context.start[,8], ins.2bp.context.start[,9],
                                   ins.2bp.context.start[,10], sep = '')
    ins.2bp.context.end <- str_split_fixed(ins.2bp.context, '', 11)[,11,drop=F]
    ins.2bp[,'TRIPLET'] <- paste(ins.2bp.context.start, ins.2bp.context.middle, ins.2bp.context.end, sep = '')
    colnames(ins.2bp)[5] <- 'CONTEXT' 
    
  }
  
  # (ii) 3 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.3bp) > 0){
    
    ins.3bp.context <- as.character(subseq(x = reference[as.character(ins.3bp[,'CHROM'])], 
                                           start = as.numeric(ins.3bp[,'POS']) - 14, 
                                           end = as.numeric(ins.3bp[,'POS']) + 15))
    ins.3bp.context.middle <- as.character(ins.3bp[,'ALT'])
    ins.3bp.context.middle <- paste0('[', str_split_fixed(ins.3bp.context.middle, '', 2)[,2,drop=F], ']')
    ins.3bp.context.start <- str_split_fixed(ins.3bp.context, '', 16)[,1:15,drop=F]
    ins.3bp.context.start <- paste(ins.3bp.context.start[,1], ins.3bp.context.start[,2], ins.3bp.context.start[,3],
                                   ins.3bp.context.start[,4], ins.3bp.context.start[,5], ins.3bp.context.start[,6],
                                   ins.3bp.context.start[,7], ins.3bp.context.start[,8], ins.3bp.context.start[,9],
                                   ins.3bp.context.start[,10], ins.3bp.context.start[,11], ins.3bp.context.start[,12], 
                                   ins.3bp.context.start[,13], ins.3bp.context.start[,14], ins.3bp.context.start[,15], sep = '')
    ins.3bp.context.end <- str_split_fixed(ins.3bp.context, '', 16)[,16,drop=F]
    ins.3bp[,'TRIPLET'] <- paste(ins.3bp.context.start, ins.3bp.context.middle, ins.3bp.context.end, sep = '')
    colnames(ins.3bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iii) 4 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.4bp) > 0){
    
    ins.4bp.context <- as.character(subseq(x = reference[as.character(ins.4bp[,'CHROM'])], 
                                           start = as.numeric(ins.4bp[,'POS']) - 19, 
                                           end = as.numeric(ins.4bp[,'POS']) + 20))
    ins.4bp.context.middle <- as.character(ins.4bp[,'ALT'])
    ins.4bp.context.middle <- paste0('[', str_split_fixed(ins.4bp.context.middle, '', 2)[,2,drop=F], ']')
    ins.4bp.context.start <- str_split_fixed(ins.4bp.context, '', 21)[,1:20,drop=F]
    ins.4bp.context.start <- paste(ins.4bp.context.start[,1], ins.4bp.context.start[,2], ins.4bp.context.start[,3],
                                   ins.4bp.context.start[,4], ins.4bp.context.start[,5], ins.4bp.context.start[,6],
                                   ins.4bp.context.start[,7], ins.4bp.context.start[,8], ins.4bp.context.start[,9],
                                   ins.4bp.context.start[,10], ins.4bp.context.start[,11], ins.4bp.context.start[,12], 
                                   ins.4bp.context.start[,13], ins.4bp.context.start[,14], ins.4bp.context.start[,15], 
                                   ins.4bp.context.start[,16], ins.4bp.context.start[,17], ins.4bp.context.start[,18], 
                                   ins.4bp.context.start[,19], ins.4bp.context.start[,20], sep = '')
    ins.4bp.context.end <- str_split_fixed(ins.4bp.context, '', 21)[,21,drop=F]
    ins.4bp[,'TRIPLET'] <- paste(ins.4bp.context.start, ins.4bp.context.middle, ins.4bp.context.end, sep = '')
    colnames(ins.4bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iv) 5+ bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.5bp) > 0){
    
    ins.5bp.context <- as.character(subseq(x = reference[as.character(ins.5bp[,'CHROM'])], 
                                           start = as.numeric(ins.5bp[,'POS']) - 499, 
                                           end = as.numeric(ins.5bp[,'POS']) + 500))
    
    ### account for the different repeat lengths: iterate
    ins.5bp.context.middle <- as.character(ins.5bp[,'ALT'])
    ins.5bp.context.middle <- str_split_fixed(ins.5bp.context.middle, '', 2)[,2,drop=F]
    ins.5bp.context.middle <- paste0('[', ins.5bp.context.middle, ']')
    ins.5bp.context.middle.lengths <- nchar(as.character(ins.5bp[,'ALT'])) - 1
    ins.5bp.context.start <- rep(NA, nrow(ins.5bp))
    for (i in 1:length(ins.5bp.context.start)){
      tmp.ins.5bp.context.start <- str_split_fixed(ins.5bp.context[i], '', 501)[,c(1 + 500-c(5*ins.5bp.context.middle.lengths[i])):500,drop=F]
      ins.5bp.context.start[i] <- paste(tmp.ins.5bp.context.start, collapse = '')
    }
    ins.5bp.context.end <- rep(NA, nrow(ins.5bp))
    for (i in 1:length(ins.5bp.context.end)){
      tmp.ins.5bp.context.end <- str_split_fixed(ins.5bp.context[i], '', 1000)[,501:c(501 + ins.5bp.context.middle.lengths[i]*6 - 1),drop=F]
      ins.5bp.context.end[i] <- paste(tmp.ins.5bp.context.end, collapse = '')
    }
    ins.5bp[,'TRIPLET'] <- paste(ins.5bp.context.start, ins.5bp.context.middle, ins.5bp.context.end, sep = '')
    colnames(ins.5bp)[5] <- 'CONTEXT'
    
  }
  
  ## summarise >=2 bp insertions at simple repeats in matrix format
  ins.greater.2bp.summary <- matrix(0, ncol = 4, nrow = 6)
  colnames(ins.greater.2bp.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## insertion size
  rownames(ins.greater.2bp.summary) <- c('0', '1', '2', '3', '4', '5+') ## number of repeats
  
  ## 2 bp
  if(nrow(ins.2bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.2bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '2 bp'] <- ins.greater.2bp.summary['5+', '2 bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '2 bp'] <- ins.greater.2bp.summary['4', '2 bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '2 bp'] <- ins.greater.2bp.summary['3', '2 bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '2 bp'] <- ins.greater.2bp.summary['2', '2 bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '2 bp'] <- ins.greater.2bp.summary['1', '2 bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '2 bp'] <- ins.greater.2bp.summary['0', '2 bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 3 bp
  if(nrow(ins.3bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.3bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '3 bp'] <- ins.greater.2bp.summary['5+', '3 bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '3 bp'] <- ins.greater.2bp.summary['4', '3 bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '3 bp'] <- ins.greater.2bp.summary['3', '3 bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '3 bp'] <- ins.greater.2bp.summary['2', '3 bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '3 bp'] <- ins.greater.2bp.summary['1', '3 bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '3 bp'] <- ins.greater.2bp.summary['0', '3 bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 4 bp
  if(nrow(ins.4bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.4bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '4 bp'] <- ins.greater.2bp.summary['5+', '4 bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '4 bp'] <- ins.greater.2bp.summary['4', '4 bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '4 bp'] <- ins.greater.2bp.summary['3', '4 bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '4 bp'] <- ins.greater.2bp.summary['2', '4 bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '4 bp'] <- ins.greater.2bp.summary['1', '4 bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '4 bp'] <- ins.greater.2bp.summary['0', '4 bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 5+ bp
  if(nrow(ins.5bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.5bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '5+ bp'] <- ins.greater.2bp.summary['5+', '5+ bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '5+ bp'] <- ins.greater.2bp.summary['4', '5+ bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '5+ bp'] <- ins.greater.2bp.summary['3', '5+ bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '5+ bp'] <- ins.greater.2bp.summary['2', '5+ bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '5+ bp'] <- ins.greater.2bp.summary['1', '5+ bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '5+ bp'] <- ins.greater.2bp.summary['0', '5+ bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 5. summarise outputs as five lists, each featuring one table per major ID category
  out <- list('1bp del' = dels.1bp.summary,
              '1bp ins' = ins.1bp.summary,
              '>2bp del' = dels.greater.2bp.summary,
              '>2bp ins' = ins.greater.2bp.summary,
              'MH dels' = dels.greater.2bp.MH.summary)
  return(out)
  
}

#Function to convert indelwald spectrum (produced by indel.spectrum) to sigfit format
indelwald.to.sigfit <- function(indelwald.spectrum){
  result <- unlist(lapply(indelwald.spectrum,c)) %>%
    as.data.frame %>%
    na.omit %>%
    set_names("count") %>%
    bind_cols(label=indel_labels,.) %>%
    pivot_wider(names_from=label,values_from=count) %>%
    as.data.frame
  return(result)
}
