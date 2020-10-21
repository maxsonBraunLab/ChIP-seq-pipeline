# pca of all samples * conds * factor and pca per factor
rule deseq2_pca:
	input:
		catalog = "data/counts_table.txt",
		contrasts = "config/contrasts.txt"
	output:
		all_pca = "data/de/all_pca.png",
		mark_pca = expand("data/de/{factor}.png", factor = set(contrasts.Factor))
	params:
		linear_model = config['linear_model']
	conda:
		"../envs/deseq2.yaml"
	threads: 8
	script:
		"../scripts/deseq2_pca.R"

rule deseq2:
	input:
		catalog = "data/counts_table.txt",
		contrasts = "config/contrasts.txt"
	output:
		# pca and ma plots
		"data/de/{mark}/{case}_{control}_upregulate.txt",
		"data/de/{mark}/{case}_{control}_downregulate.txt",
	params:
		linear_model = config['linear_model']
	conda:
		"../envs/deseq2.yaml"
	script:
		"../scripts/deseq2_pairwise.R"