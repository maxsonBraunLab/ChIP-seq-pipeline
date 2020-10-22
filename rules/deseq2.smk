# for each factor input, apply DESeq2 on Conditions + PCA output.
rule deseq2_plots:
	input:
		catalog = "results/counts/{factor}_counts.txt",
		contrasts = "config/contrasts.txt"
	output:
		factor_pca = "results/de/plots/{factor}_pca.png",
		factor_ma = "results/de/plots/{factor}_ma.png"
	params:
		factor = "{factor}"
	conda:
		"../envs/deseq2.yaml"
	threads: 2
	script:
		"../scripts/deseq2_pca.R"

rule deseq2_pairwise:
	input:
		catalog = "results/counts/{factor}_counts.txt",
		contrasts = "config/contrasts.txt"
	output:
		# deseq2-normalized counts
		normCounts = "results/de/{factor}/{factor}_deseq2_norm_counts.txt",
		lnormCounts = "results/de/{factor}/{factor}_deseq2_lognorm_counts.txt",
		# differential peaks + gene ontology + ma plots
		d = directory("results/de/{factor}")
	params:
		linear_model = config['linear_model'],
		factor = "{factor}"
	conda:
		"../envs/deseq2.yaml"
	threads: 4
	log: "logs/deseq2_pairwise/{factor}.log"
	script:
		"../scripts/deseq2_pairwise.R"