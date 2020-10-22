library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tibble)
library(stringr)

parallel <- FALSE
if (snakemake@threads > 1) {
	library("BiocParallel")
	# setup parallelization
	register(MulticoreParam(snakemake@threads))
	parallel <- TRUE
}

# import data ----------------------------------------------------------
dat <- read.delim(snakemake@input[["catalog"]], header = T, stringsAsFactors=F)
md <- read.delim(snakemake@input[["contrasts"]], header = T, stringsAsFactors=F)
md <- md %>% filter(Factor == snakemake@params[["factor"]])
dat_bed <- dat[,1:6]

counts <- dat[,md$SampleID]
rownames(counts) <- dat$V4
rownames(md) <- md$SampleID

# deseq2 ----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
	countData = as.matrix(counts),
	colData = md,
	design = ~ Condition)

dds <- DESeq(dds, parallel = parallel)
res <- results(dds)

# export deseq2-normalized counts ---------------------------------------
normCounts <- counts(dds, normalized=TRUE)
normCounts <- 0.00001+(normCounts) # ensure nonzero 
lnormCounts <- log2(normCounts)
lnormCounts <- as.data.frame(lnormCounts)
lnormCounts$ID <- counts$peak
write.table(normCounts, snakemake@output[['normCounts']], quote = FALSE, sep = "\t")
write.table(lnormCounts, snakemake@output[['lnormCounts']], quote = FALSE, sep = "\t")

# create + export all unique contrasts ----------------------------------
combs <- noquote(t(combn(unique(md$Condition), 2)))
colnames(combs) <- c("c1", "c2")
print("Contrast Combinations for Conditions")
print(combs)

factor <- snakemake@params[["factor"]]

if (!dir.exists(paste0("results/de/", factor))) {
	dir.create(paste0("results/de/", factor))
}

for (i in 1:nrow(combs)) {
	# define contrasts
	c1 = combs[i, "c1"]
	c2 = combs[i, "c2"]
	# output differential interval tables
	outfile <- paste0("results/de/", factor, "/", factor, "-", c1, "-vs-", c2, ".txt")
	temp_res <- results(dds, c("Condition", c1, c2)) %>% 
			as.data.frame() %>% 
			rownames_to_column("V4") %>% 
			arrange(padj) %>% 
			inner_join(., dat_bed, by = "V4") %>% 
			select(Chr, start, stop, V4, V5, V6, everything())
	print(paste(c1, c2))
	print(head(temp_res))
	write.table(temp_res, outfile, quote = FALSE, sep = "\t", row.names = F)
	# output MA plots
	# GO analysis
}