rule frip_count:
	input:
		bam = "samples/bams/{sample}.sorted.bam",
		peaks = "samples/macs/{sample}_peaks.narrowPeak"
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
		"data/qc/frip.html"
	run:
		pd.options.plotting.backend = "plotly"
		# import data. row = outside, inside, and ratio. cols = samples. 
		df = pd.concat([ pd.read_csv(i) for i in sorted(input) ], axis = 1)
		df = df.rename(index={0: 'outside', 1: 'inside'})
		df.loc['ratio'] = df.loc['inside'] / df.loc['outside']
		# plot graph. plot ratio as bottom as percent, and plot to max value of 1.
		fig = go.Figure(data=[
			go.Bar(name='FRiP', x=df.columns, y=df.loc['ratio'], marker_color='rgb(255,201,57)'),
			go.Bar(name='outside', x=df.columns, y= ([1] * df.shape[1]) - df.loc['ratio'], marker_color='rgb(0,39,118)')])
		# Change the bar mode
		fig.update_layout(barmode='stack', title='Fraction of Reads in Peaks by Sample', xaxis_tickfont_size=14,
			yaxis=dict(title='Fraction of reads in peaks', titlefont_size=16, tickfont_size=14),
			xaxis=dict(title='Samples'))
		fig.write_html(str(output))

# fastqc-screen

# library complexity

# deduplication rate