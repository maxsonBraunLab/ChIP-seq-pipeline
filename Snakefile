__author__ = "Breeshey Roskams-Hieter"
__email__ = "roskamsh@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""

import datetime
import sys
import os
import pandas as pd
import plotly as plt
import plotly.express as px
import plotly.graph_objects as go

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"config.yaml"
project_id = config["project_id"]
seq_type = config["seq_type"]
md = pd.read_table(config["samples"], index_col=["SampleID"], dtype=str)
contrasts = pd.read_csv(config['contrasts'], sep = "\t")

if config["seq_type"]=="SE":
    CASES, = glob_wildcards("samples/cases/{sample}.fastq.gz")
    COND_REPS, FACTORS, = glob_wildcards("samples/cases/{cond_rep}_{factor}.fastq.gz")
else:
    CASES, = glob_wildcards("samples/cases/{sample}_R1.fastq.gz")
    COND_REPS, FACTORS, = glob_wildcards("samples/cases/{cond_rep}_{factor}_R1.fastq.gz")
# CASES glob for "{cond}{rep}_{factor}"
# COND_REPS glob for "{cond}{rep}" only
# FACTORS glob for "{factor}" only

if config["seq_type"]=="SE":
    CONTROLS, = glob_wildcards("samples/controls/{sample}.fastq.gz")
else:
    CONTROLS, = glob_wildcards("samples/controls/{sample}_R1.fastq.gz")
CASES = sorted(CASES)
CONTROLS = sorted(CONTROLS)
COND_REPS = sorted(set(COND_REPS))

## multiple samples may use the same control input/IgG files
CONTROLS_UNIQUE = list(set(CONTROLS))

SAMPLES = CASES + CONTROLS_UNIQUE
FACTORS = set(contrasts.Factor)

rule_dirs = ['mapReads','makeTracks','bb2bed','call_peaks']
for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

def get_peak(wildcards):
    if md.loc[wildcards.sample,["Peak_call"]].values == "broad":
        return "broad"
    elif md.loc[wildcards.sample,["Peak_call"]].values == "narrow":
        return "narrow"
    else:
        return "ERROR"

def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

for case in CASES:
    message("case " + case + " will be expanded")

DB, CELL_LINE, FILE, = glob_wildcards("/home/groups/MaxsonLab/kongg/chip_seq/data/beds/{db}/{cell_line}/{file}.bed.gz")
MERGED, = glob_wildcards("custom/bams/{sample}.mapped.dedup.sorted.bam")

# snakemake -j 64 --use-conda --rerun-incomplete --latency-wait 60 --keep-going --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N} -o {cluster.o} -e {cluster.e} -t {cluster.t} -J {cluster.J} -c {cluster.c} --mem={cluster.mem}" -s Snakefile

# snakemake -j 99 --use-conda --rerun-incomplete --latency-wait 60 --keep-going --profile slurm --cluster-config cluster2.json

localrules: frip_plot, consensus_peaks, merge_counts, plotHeatmap_merged_plotly

rule all:
    input:
        # pre-process = align + make tracks + macs2 + homer
        expand("samples/bams/{sample}.mapped.dedup.sorted.bam", sample = SAMPLES),
        expand("samples/bigBed/{sample}.all.bb", sample = SAMPLES),
        expand("samples/bigwig/{sample}.bw", sample = SAMPLES),
        expand("samples/macs/{sample}/{sample}_peaks.narrowPeak", sample  = CASES),
        expand("results/motifs/{sample}/homerResults.html", sample = CASES),
        # pre-DESeq2 = consensus peaks + counts table
        "samples/macs/all_peaks.bed",
        "samples/macs/consensus_stats.txt",
        expand("samples/macs/{factor}_peaks.bed", factor = FACTORS),
        expand("samples/macs/counts/{sample}_{factor}.txt", sample = COND_REPS, factor = FACTORS),
        # DESeq2 = pca + normalized counts + analyze contrasts
        expand("results/de/plots/{factor}_pca.png", factor = FACTORS),
        expand(["results/de/{factor}/{factor}_deseq2_norm_counts.txt",
                "results/de/{factor}/{factor}_deseq2_lognorm_counts.txt"], factor = FACTORS),
        directory(expand("results/de/{factor}", factor = FACTORS)),
        # quality control = frip + fastqc-screen + dispersion of coverage + library complexity + duplicate rate
        expand(["samples/qc/fastq_screen/{cases}/{cases}_{dir}_screen.{ext}",
                "samples/qc/fastq_screen/{controls}/{controls}_{dir}_screen.{ext}"],
                dir = ["R1", "R2"],
                cases = CASES,
                controls = CONTROLS_UNIQUE,
                ext = ["png", "txt", "html"]),
        "results/qc/frip.html",
        # custom analysis -------------------------------------------------
        # merge bam replicates + convert to bigwig + compute consensus peaks heatmaps + plot in one stacked plot.
        expand("custom/bams/{cond}_{factor}.mapped.dedup.sorted.bam", 
            cond = set( [i.split("_")[0][:-1] for i in CASES] ), 
            factor = set(contrasts.Factor) ),
        expand("custom/bigwigs/{sample}.bw", sample = MERGED),
        expand("custom/qc/computeMatrix_merged/{sample}_{factor}_matrix.{ext}", 
            sample = [i.split("_")[0] for i in MERGED], 
            factor = FACTORS, 
            ext = ["gz", "txt"]),
        expand("custom/results/plotHeatmap/{sample}_{factor}_heatmap.png", 
            sample = [i.split("_")[0] for i in MERGED],
            factor = FACTORS),
        expand("custom/results/{factor}_heatmap.html",
            factor = FACTORS)

include: "rules/align.smk"
include: "rules/peaks.smk"
include: "rules/qc.smk"
include: "rules/deseq2.smk"

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