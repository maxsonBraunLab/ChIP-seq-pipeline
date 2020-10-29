# for each factor input, apply DESeq2 on Conditions + PCA output.
rule deseq2_plots:
	input:
		catalog = "results/counts/{factor}_counts.txt",
		contrasts = "config/contrasts.txt"
	output:
		factor_pca = "results/de/plots/{factor}_pca.png"
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
		directory("results/de/{factor}"),
		# deseq2-normalized counts
		normCounts = "results/de/{factor}/{factor}_deseq2_norm_counts.txt",
		lnormCounts = "results/de/{factor}/{factor}_deseq2_lognorm_counts.txt",
		# differential peaks (contains DE, up, down, and regions).
		stats = temp("results/de/{factor}_stats.txt"),
		all_sig_intervals = "results/de/{factor}/{factor}_sig_intervals.bed"
	params:
		linear_model = config['linear_model'],
		factor = "{factor}",
		significance = config["significance"]
	conda:
		"../envs/deseq2.yaml"
	log: 
		"logs/deseq2_pairwise/{factor}.log"
	threads: 4
	script:
		"../scripts/deseq2_pairwise.R"

rule de_stats:
	input:
		expand("results/de/{factor}_stats.txt", factor = FACTORS)
	output:
		"results/de/de_stats.txt"
	params:
		header = "\t".join(['factor', 'cond1', 'cond2', 'tot_peaks', 'sig_peaks', 'sig_up', 'sig_down'])
	shell:
		"cat {input} | sort > {output}; sed -i '1i {params.header}' {output}"
