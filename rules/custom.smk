# custom analyses: replicate-merged heatmaps ----------------------------------------
rule merge_bams:
    input:
        rep1 = "samples/bams/{cond}1_{factor}.mapped.dedup.sorted.bam",
        rep2 = "samples/bams/{cond}2_{factor}.mapped.dedup.sorted.bam"
    output:
        "custom/bams/{cond}_{factor}.mapped.dedup.sorted.bam"
    conda:
        "envs/deeptools.yaml"
    log:
        "logs/merge_bams/{cond}_{factor}.log"
    threads: 4
    shell:
        """
        samtools merge -@ {threads} -f {output} {input.rep1} {input.rep2}
        sleep 10
        samtools index {output}
        echo $(date) 'samtools merge -@ {threads} -f {output} {input}' > {log}
        """

rule merged_bw:
    input:
        "custom/bams/{sample}.mapped.dedup.sorted.bam"
    output:
        "custom/bigwigs/{sample}.bw"
    conda:
        "envs/deeptools.yaml"
    threads: 16
    shell:
        "bamCoverage -b {input} -o {output} -p {threads}"

rule computeMatrix_merged:
    input:
        bw = "custom/bigwigs/{sample}_{factor}.bw",
        r = "samples/macs/{factor}_peaks.bed",
        bl = config["blacklist"]
    output:
        gz = "custom/qc/computeMatrix_merged/{sample}_{factor}_matrix.gz",
        matx = "custom/qc/computeMatrix_merged/{sample}_{factor}_matrix.txt"
    params:
        before = config["before_region"],
        after = config["after_region"],
    conda:
        "envs/deeptools.yaml"
    threads: 16
    shell:
        "computeMatrix reference-point -S {input.bw} -R {input.r} -a {params.after} -b {params.before} -bl {input.bl} -p {threads}  --referencePoint center --smartLabels -o {output.gz} --outFileNameMatrix {output.matx}"

rule plotHeatmap_merged:
    input:
        rules.computeMatrix_merged.output.gz
    output:
        "custom/results/plotHeatmap/{sample}_{factor}_heatmap.png"
    conda:
        "envs/deeptools.yaml"
    threads: 4
    shell:
        "plotHeatmap -m {input} -o {output} --colorMap RdBu --boxAroundHeatmaps no"

rule plotHeatmap_merged_plotly:
    input:
        expand("custom/qc/computeMatrix_merged/{sample}_{{factor}}_matrix.txt", sample = set([i.split("_")[0] for i in MERGED]))
    output:
        "custom/results/{factor}_heatmap.html"
    run:
        pd.options.plotting.backend = "plotly"
        signal_df = []
        for file in input:
            file = str(file)
            basename = os.path.basename(file).split(".")[0] # no extension
            # grab the bin size
            with open(file, "r") as fi:
                for i, line in enumerate(fi):
                    if i == 1:
                        bin_size = line.split("\t")[3].split(":")[1]
            # read in data as table, take mean of all signal at consensus intervals, collate col to signal_df.
            temp_df = pd.read_csv(file, sep = "\t", skiprows = 3, header = None)
            avg_signal = temp_df.mean(0)
            avg_signal = avg_signal.rename("_".join(basename.split("_")[0:2]))
            signal_df.append(avg_signal)
        signal_df = pd.concat(signal_df, axis = 1)
        # re-zero x-axis. 0 is center.
        half_window = int(len(signal_df.index) / 2)
        whole_window = list(range(-half_window, half_window))
        whole_window = [i * int(bin_size) for i in whole_window]
        new_index = dict()
        for i, j in zip(signal_df.index, whole_window):
            new_index[i] = j
        signal_df = signal_df.rename(new_index)
        print("signal dataframe")
        print(signal_df)
        # Create title and plot
        title = "Effect of drug treatment with " + basename.split("_")[1] + " factor"
        fig = signal_df.plot()
        fig.update_layout( 
            title=title, 
            xaxis_title='Distance away from consensus intervals (bp)', 
            yaxis_title="Signal", 
            legend_title_text='Treatment_Factor')
        fig.write_html(str(output))