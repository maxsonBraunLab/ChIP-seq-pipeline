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

# pca -------------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
pca <- DESeq2::plotPCA(vsd, intgroup = "Condition") + geom_text_repel(aes(label=name))
ggsave(snakemake@output[['factor_pca']], pca, width = 16, height = 9, dpi = 300, units = "in")