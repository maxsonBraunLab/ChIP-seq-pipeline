# Goal: find overlapping intervals among replicates, split by ChIP-Seq factor (histone or TF). #
library(dplyr)
library(plyranges)
library(stringr)
library(purrr)
library(GenomicRanges)
library(tidyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

#This is the location for our merged narrow peaks file
np_file <- args[1]
blacklist_file <- args[2]
presence_in_samples <- as.numeric(args[3]) # e.g. 2
genome_name <- args[4]

#This command reads in the merged broadnarrowpeak files and adds the information about the genome
#we used 
#peaks <- read_narrowPeaks(np_file, genome_info = genome_name)
peaks <- fread(np_file)

#lets load the peaks file, define cond and factor from name column, and convert it to a granges object
colnames(peaks) <- c("seqnames", "start", "end", "name", "score", "no_strand", "signalValue", "pValue","qvalue", "peak")
peaks <- peaks %>%
	separate(name, into = c("rep", "factor"), sep = "_", remove = FALSE, extra = "drop") %>% 
	mutate(cond = str_sub(rep, end = -2))
peaks <- as_granges(peaks) %>% set_genome_info(genome = genome_name)


#now lets read in the blacklist file so we can remove blacklisted regions
bl <- read.delim(blacklist_file, header=FALSE)
colnames(bl) <- c("seqnames","start","end","reason")
blacklist <- as_granges(bl) %>% set_genome_info(genome = genome_name)

#filter out only the chromosomes that we want and that we care about
#This filtering command keeps all chromosomes that are numbered with nothing after the
#number, as well as the X and Y chromosomes
filtered_peaks <- peaks %>% filter(., str_detect(seqnames,"^chr[\\dXY]+$"))

# define chr seqlengths. Note chrM not included for epigenetics-focused analysis.
if (genome_name == "hg38") {
	sl <- c(chr1=248956422,
			chr10=133797422,
			chr11=135086622,
			chr12=133275309,
			chr13=114364328,
			chr14=107043718,
			chr15=101991189,
			chr16=90338345,
			chr17=83257441,
			chr18=80373285,
			chr19=58617616,
			chr2=242193529,
			chr20=64444167,
			chr21=46709983,
			chr22=50818468,
			chr3=198295559,
			chr4=190214555,
			chr5=181538259,
			chr6=170805979,
			chr7=159345973,
			chr8=145138636,
			chr9=138394717,
			chrX=156040895,
			chrY=57227415)
} else if (genome_name == "mm10") {
	sl <- c(chr1=195471971,
			chr10=130694993,
			chr11=122082543,
			chr12=120129022,
			chr13=120421639,
			chr14=124902244,
			chr15=104043685,
			chr16=98207768,
			chr17=94987271,
			chr18=90702639,
			chr19=61431566,
			chr2=182113224,
			chr3=160039680,
			chr4=156508116,
			chr5=151834684,
			chr6=149736546,
			chr7=145441459,
			chr8=129401213,
			chr9=124595110,
			chrX=171031299,
			chrY=91744698)
}

# for each factor and for each condition (cond), find peaks that appear in presence_in_samples number of replicates.
Factors <- unique(filtered_peaks$factor)

stats_catalog <- data.frame()

for (f in Factors) {
	# given all peaks, subset by factor, split by condition, count no. of overlaps between peaks per reps
	outfile <- paste0('samples/macs/', f, '_peaks.bed')
	print(paste("Finding consensus peak in >=", presence_in_samples, "number of replicates from factor", f))
	temp_df <- filtered_peaks %>% filter(factor == f)
	temp_split <- temp_df %>% split(., temp_df$cond) # split factor by conditions
	seqlengths(temp_split) <- sl
	temp_split <- lapply(temp_split, function(x) {
			compute_coverage(x) %>% 
			filter(score >= presence_in_samples) %>% 
			reduce_ranges() %>% 
			filter_by_non_overlaps(blacklist)
			}) %>% GRangesList()
	# print message if no overlaps in at least one cond.
	if (any(unlist( lapply(temp_split, length) ) == 0 )) {
		message("Factor ", f, " has no consensus peaks within a condition.")
		print(unlist( lapply(temp_split, length) ))
	}
	# within a factor, combine all condition's peaks again and merge intervals. name peaks
	consensus_per_factor <- unlist(temp_split) %>% reduce_ranges()
	consensus_per_factor <- consensus_per_factor %>% as.data.frame() %>% mutate(name = paste0(f, "_peak_",row_number())) %>% as_granges()
	write_bed(consensus_per_factor, outfile)


	# write consensus statistics over all replicates.
	rep_split <- temp_df %>% split(., temp_df$rep)
	stats_list <- lapply(rep_split, function(x) {
		replicate <- toString(unique(x$rep))
		condition <- unique(x$cond)
		peaks_in_replicate <- as.integer(nrow(as.data.frame(x)))
		peaks_in_consensus <- nrow(as.data.frame(  join_overlap_inner(x, consensus_per_factor)  ))
		peaks_in_consensus_over_peaks_in_consensus <- peaks_in_consensus / peaks_in_replicate
		number_consensus_peaks <- nrow(as.data.frame(consensus_per_factor))
		temp_stats <- data.frame(
			factor=f,
			cond=condition,
			rep=replicate,
			tot_peaks=peaks_in_replicate,
			peaks_in_consensus=peaks_in_consensus,
			prop_consensus_peaks_in_tot_peaks=peaks_in_consensus_over_peaks_in_consensus,
			tot_consensus_peaks=number_consensus_peaks)
	})
	stats_per_factor <- do.call(rbind, stats_list)
	stats_catalog <- rbind(stats_catalog, stats_per_factor)
}

print("Consensus statistics")
print(stats_catalog)
fwrite(stats_catalog, "samples/macs/consensus_stats.txt", sep = "\t", quote = FALSE, )