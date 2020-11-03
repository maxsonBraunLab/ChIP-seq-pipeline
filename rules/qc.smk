# pre-alignment QC ------------------------------------------------------

# fastq-screen
rule fastq_screen_cases:
	input:
		r1 = "samples/cases/{sample}_R1.fastq.gz",
		r2 = "samples/cases/{sample}_R2.fastq.gz"
	output:
		"results/qc/fastq_screen/{sample}/{sample}_R1_screen.txt",
		"results/qc/fastq_screen/{sample}/{sample}_R1_screen.png",
		"results/qc/fastq_screen/{sample}/{sample}_R1_screen.html",
		"results/qc/fastq_screen/{sample}/{sample}_R2_screen.txt",
		"results/qc/fastq_screen/{sample}/{sample}_R2_screen.png",
		"results/qc/fastq_screen/{sample}/{sample}_R2_screen.html"
	params:
		conf = config["fastq_screen_conf"],
		outdir = "results/qc/fastq_screen/{sample}"
	conda:
		"../envs/fastq_screen.yaml"
	message: " -- Screening {wildcards.sample} --"
	log: "logs/fastq_screen/{sample}.log"
	threads: 16
	shell:
		"fastq_screen --aligner bowtie2 --threads {threads} --conf {params.conf} --outdir {params.outdir} {input.r1} {input.r2}"

rule fastq_screen_controls:
	input:
		r1 = "samples/controls/{sample}_R1.fastq.gz",
		r2 = "samples/controls/{sample}_R2.fastq.gz"
	output:
		"results/qc/fastq_screen/{sample}/{sample}_R1_screen.txt",
		"results/qc/fastq_screen/{sample}/{sample}_R1_screen.png",
		"results/qc/fastq_screen/{sample}/{sample}_R1_screen.html",
		"results/qc/fastq_screen/{sample}/{sample}_R2_screen.txt",
		"results/qc/fastq_screen/{sample}/{sample}_R2_screen.png",
		"results/qc/fastq_screen/{sample}/{sample}_R2_screen.html"
	params:
		conf = config["fastq_screen_conf"],
		outdir = "results/qc/fastq_screen/{sample}"
	conda:
		"../envs/fastq_screen.yaml"
	message: " -- Screening {wildcards.sample} --"
	log: "logs/fastq_screen/{sample}.log"
	threads: 16
	shell:
		"fastq_screen --aligner bowtie2 --threads {threads} --conf {params.conf} --outdir {params.outdir} {input.r1} {input.r2}"

# post-alignment QC -----------------------------------------------------
rule frip_count:
	input:
		bam = "samples/bams/{sample}.mapped.dedup.sorted.bam",
		peaks = "samples/macs/{sample}/{sample}_peaks.narrowPeak"
	output:
		"samples/qc/frip/{sample}_stats.txt"
	conda:
		"../envs/chip.yaml"
	shell:
		"""
		rip=$(bedtools sort -i {input.peaks} | bedtools merge -i stdin | bedtools intersect -u -a {input.bam} -b stdin -ubam | samtools view -c); \
		counts=$(samtools view -c {input.bam}); \
		echo -e "{wildcards.sample}\n$counts\n$rip" > {output}
		"""

rule frip_plot:
	input:
		expand("samples/qc/frip/{sample}_stats.txt", sample = CASES)
	output:
		"results/qc/frip.html"
	run:
		pd.options.plotting.backend = "plotly"
		# import data. row = outside, inside, and ratio. cols = samples. 
		df = pd.concat([ pd.read_csv(i) for i in sorted(input) ], axis = 1)
		df = df.rename(index={0: 'outside', 1: 'inside'})
		df.loc['ratio'] = df.loc['inside'] / df.loc['outside']
		# plot graph. plot ratio as bottom as percent, and plot to max value of 1.
		fig = go.Figure(data=[
			go.Bar(name='inside_peaks', x=df.columns, y=df.loc['ratio'], marker_color='rgb(255,201,57)'),
			go.Bar(name='outside_peaks', x=df.columns, y= ([1] * df.shape[1]) - df.loc['ratio'], marker_color='rgb(0,39,118)')])
		# Change the bar mode
		fig.update_layout(barmode='stack', title='Fraction of Reads in Peaks by Sample', xaxis_tickfont_size=14,
			yaxis=dict(title='Fraction of reads in peaks', titlefont_size=16, tickfont_size=14),
			xaxis=dict(title='Samples'))
		fig.write_html(str(output))

rule chip_screen:
	input:
		factor_peaks = "samples/macs/{factor}_peaks.bed",
		chip_db = "/home/groups/MaxsonLab/kongg/chip_seq/data/beds/{db}/{cell_line}/tf_chip.bed.gz"
	output:
		"results/chip_screen/{db}/{cell_line}/{factor}_screen.bed.gz"
	conda:
		"../envs/chip.yaml"
	shell:
		"bedtools intersect -wa -wb -f 0.50 -a {input.factor_peaks} -b {input.chip_db} | awk -F '[\t_]' '$4 ~ $12 || $12 ~ $4 {{print}}' | gzip > {output}"
# use tab and "_" as delimiters, then query if either our antibody matches db antibody, or other way around.

rule chip_fisher:
	input:
		factor_peaks = "samples/macs/{factor}_peaks.bed",
		chip_db = "/home/groups/MaxsonLab/kongg/chip_seq/data/beds/{db}/{cell_line}/tf_chip.bed.gz"
	output:
		"results/chip_screen/{db}/{cell_line}/{factor}_screen.out"
	params:
		chrom_sizes = config["chrom_sizes"]
	conda:
		"../envs/chip.yaml"
	shell:
		"bedtools merge -i {input.chip_db} | bedtools fisher -f 0.50 -a {input.factor_peaks} -b stdin -g {params.chrom_sizes} > {output}"
# merge all db intervals into one interval to reduce multiple intersects + account for different size db's.

def get_control_bam(wildcards):
    return md.loc[(wildcards.sample),["Control_bam"]]

rule plotFingerprint_factor:
	input:
		bam = "samples/bams/{sample}.mapped.dedup.sorted.bam",
		control = get_control_bam,
		bl = config["blacklist"]
	output:
		"samples/qc/{factor}_fingerprint.png"
	conda:
		"../envs/deeptools.yaml"
	threads: 8
	shell:
		"plotFingerprint -b {input.bam} {input.control} -o {output} -bl {input.bl} --smartLabels -p {threads}"

rule computeMatrix_factor:
	input:
		bw = expand("samples/bigwig/{sample}_{{factor}}.bw", sample = COND_REPS),
		r = "samples/macs/{factor}_peaks.bed",
		bl = config["blacklist"]
	output:
		gz = "samples/qc/computeMatrix_factor/{factor}_matrix.gz",
		matx = "samples/qc/computeMatrix_factor/{factor}_matrix.txt"
	params:
		before = config["before_region"],
		after = config["after_region"],
	conda:
		"../envs/deeptools.yaml"
	threads: 16
	shell:
		"computeMatrix reference-point -S {input.bw} -R {input.r} \
		-a {params.after} -b {params.before} \
		-bl {input.bl} -p {threads} --smartLabels -o {output.gz} --outFileNameMatrix {output.matx}"

rule plotHeatmap_factor:
	input:
		rules.computeMatrix_factor.output
	output:
		"results/qc/plotHeatmap/{factor}_all_samples.png"
	conda:
		"../envs/deeptools.yaml"
	threads: 4
	shell:
		"plotHeatmap -m {input} -o {output} --colorMap RdBu"
