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

counts <- dat[,md$SampleID]
rownames(counts) <- dat$V4
rownames(md) <- md$SampleID

# deseq2 ----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
	countData = as.matrix(counts),
	colData = md,
	design = as.formula(snakemake@params[["linear_model"]]) )

dds <- DESeq(dds, parallel = parallel)
res <- results(dds)

# pca -------------------------------------------------------------------
# pca = (1 pca with provided labels + all marks + all reps) + (1 pca per mark with all reps)
vsd <- vst(dds, blind = FALSE)
intgroups <- unlist(str_extract_all(snakemake@params[["linear_model"]], "[:alpha:]+")) # select only strings in the linear_model.

# all_pca. should cluster by TF
all_pca <- DESeq2::plotPCA(vsd, intgroup = intgroups) + geom_text_repel(aes(label=name))
ggsave(snakemake@output[['all_pca']], all_pca, width = 16, height = 9, dpi = 300, units = "in")

# pca per factor. should cluster by cond + reps.
# for each factor, subset the metadata table + counts table and rerun through deseq2.
for (factor in unique(md$Factor)) {
	print(paste("factor:", factor))
	temp_md <- md %>% 
		rownames_to_column('info') %>% 
		filter(Factor == factor) %>% 
		column_to_rownames('info')
	temp_counts <- dat[,temp_md$SampleID]
	rownames(temp_counts) <- dat$V4

	temp_dds <- DESeqDataSetFromMatrix(
		countData = as.matrix(temp_counts),
		colData = temp_md,
		design = ~ Condition )

	temp_dds <- DESeq(temp_dds)
	temp_res <- results(temp_dds)
	temp_vsd <- vst(temp_dds, blind = FALSE)
	temp_pca <- DESeq2::plotPCA(temp_vsd, intgroup = "Condition") +
		geom_text_repel(aes(label=name)) +
		ggtitle(factor)
	outfile <- paste0("data/de/", factor, ".png")
	ggsave(outfile, temp_pca, width = 16, height = 9, dpi = 300, units = "in")
}
# can parallelize with foreach but generally pretty fast.

# foreach(factor = unique(md$Factor)) %do% {
# 	print(paste("factor:", factor))
# 	temp_md <- md %>% 
# 		rownames_to_column('info') %>% 
# 		filter(Factor == factor) %>% 
# 		column_to_rownames('info')
# 	temp_counts <- dat[,temp_md$SampleID]
# 	rownames(temp_counts) <- dat$V4

# 	temp_dds <- DESeqDataSetFromMatrix(
# 		countData = as.matrix(temp_counts),
# 		colData = temp_md,
# 		design = ~ Condition )

# 	temp_dds <- DESeq(temp_dds)
# 	temp_res <- results(temp_dds)
# 	temp_vsd <- vst(temp_dds, blind = FALSE)
# 	temp_pca <- DESeq2::plotPCA(temp_vsd, intgroup = "Condition") + 
# 		geom_text_repel(aes(label=name)) +
# 		ggtitle(factor)
# 	outfile <- paste0("data/de/", factor, ".png")
# 	ggsave(outfile, temp_pca, width = 16, height = 9, dpi = 300, units = "in")
# }