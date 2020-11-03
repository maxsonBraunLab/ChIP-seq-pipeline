rule deseq2:
	input:
		catalog = "results/counts/{factor}_counts.txt",
		contrasts = "config/contrasts.txt"
	output:
		directory("results/de/{factor}"),
		factor_pca = "results/de/plots/{factor}_pca.png",
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
		"../scripts/deseq2.R"

rule de_stats:
	input:
		expand("results/de/{factor}_stats.txt", factor = FACTORS)
	output:
		"results/de/de_stats.txt"
	params:
		header = "\t".join(['factor', 'cond1', 'cond2', 'tot_peaks', 'sig_peaks', 'sig_up', 'sig_down'])
	shell:
		"cat {input} | sort > {output}; sed -i '1i {params.header}' {output}"

rule essential_report:
	input:
		# pipeline QC metrics
		align_stats = "results/qc/align_stats.txt",
		consensus_stats = "samples/macs/consensus_stats.txt",
		de_stats = "results/de/de_stats.txt",
		frip_folder = directory("samples/qc/frip")
	output:
		"results/essential_report.html"
	conda:
		"../envs/report.yaml"
	params:
		# user info
		title = config["title"],
		authors = config["authors"],
		# written content
		intro = config["intro"],
		analysis = config["analysis"],
		takeaways = config["takeaways"],
		notes = config["notes"],
		# config files
		peak_md = config["samples"],
		contrasts = config["contrasts"],
	script:
		"../scripts/essential_report.Rmd"