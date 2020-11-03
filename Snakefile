__author__ = "Breeshey Roskams-Hieter"
__email__ = "roskamsh@ohsu.edu"
__license__ = "MIT"
# improved by Garth Kong

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

essential_report = []
if config["gen_report"] == True: essential_report.append("results/essential_report.html"); message("Generating essential report")

# custom analysis
MERGED, = glob_wildcards("custom/bams/{sample}.mapped.dedup.sorted.bam")

# snakemake -j 99 --use-conda --rerun-incomplete --latency-wait 60 --keep-going --profile slurm --cluster-config cluster2.json

localrules: frip_plot, de_stats, align_stats, merge_counts, essential_report

rule all:
    input:
        # pre-process = align + make tracks + macs2 + homer
        expand("samples/bams/{sample}.mapped.dedup.sorted.bam", sample = SAMPLES),
        expand("samples/bigBed/{sample}.all.bb", sample = SAMPLES),
        expand("samples/bigwig/{sample}.bw", sample = SAMPLES),
        expand("samples/macs/{sample}/{sample}_peaks.narrowPeak", sample  = CASES),
        expand("results/motifs/{sample}/homerResults.html", sample = CASES),
        "results/qc/align_stats.txt",
        # pre-DESeq2 = consensus peaks + counts table
        "samples/macs/all_peaks.bed",
        "results/qc/consensus_stats.txt",
        expand("samples/macs/{factor}_peaks.bed", factor = FACTORS),
        expand("samples/macs/counts/{sample}_{factor}.txt", sample = COND_REPS, factor = FACTORS),
        # DESeq2 = pca + normalized counts + analyze contrasts
        "results/de/de_stats.txt",
        expand("results/de/plots/{factor}_pca.png", factor = FACTORS),
        expand(["results/de/{factor}/{factor}_deseq2_norm_counts.txt",
                "results/de/{factor}/{factor}_deseq2_lognorm_counts.txt"], factor = FACTORS),
        directory(expand("results/de/{factor}", factor = FACTORS)),
        # quality control = frip + fastqc-screen + duplicate rate
        "results/qc/frip.html",
        expand(["results/qc/fastq_screen/{cases}/{cases}_{dir}_screen.{ext}",
                "results/qc/fastq_screen/{controls}/{controls}_{dir}_screen.{ext}"],
                dir = ["R1", "R2"],
                cases = CASES,
                controls = CONTROLS_UNIQUE,
                ext = ["png", "txt", "html"]),
        # essential report
        essential_report,
        # custom analysis -------------------------------------------------
        # merge bam replicates + convert to bigwig + compute consensus peaks heatmaps + plot in one stacked plot.
        # expand("custom/bams/{cond}_{factor}.mapped.dedup.sorted.bam", 
        #     cond = set( [i.split("_")[0][:-1] for i in CASES] ), 
        #     factor = set(contrasts.Factor) ),
        # expand("custom/bigwigs/{sample}.bw", sample = MERGED),
        # expand("custom/qc/computeMatrix_merged/{sample}_{factor}_matrix.{ext}", 
        #     sample = [i.split("_")[0] for i in MERGED], 
        #     factor = FACTORS, 
        #     ext = ["gz", "txt"]),
        # expand("custom/results/plotHeatmap/{sample}_{factor}_heatmap.png", 
        #     sample = [i.split("_")[0] for i in MERGED],
        #     factor = FACTORS),
        # expand("custom/results/{factor}_heatmap.html",
        #     factor = FACTORS),
        # expand("custom/qc/computeMatrix_LSD1_DE/{sample}_LSD1_matrix.{ext}",
        #     sample = [i.split("_")[0] for i in MERGED],
        #     ext = ["gz", "txt"])

# workflow in this order
include: "rules/align.smk"
include: "rules/peaks.smk"
include: "rules/qc.smk"
include: "rules/deseq2_report.smk"
include: "rules/custom.smk"