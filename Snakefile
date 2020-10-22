__author__ = "Breeshey Roskams-Hieter"
__email__ = "roskamsh@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""

import datetime
import sys
import os
import pandas as pd
import plotly as plt
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

# snakemake -j 64 --use-conda --rerun-incomplete --latency-wait 60 --keep-going --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N} -o {cluster.o} -e {cluster.e} -t {cluster.t} -J {cluster.J} -c {cluster.c} --mem={cluster.mem}" -s Snakefile

# snakemake -j 99 --use-conda --rerun-incomplete --latency-wait 60 --keep-going --profile slurm --cluster-config cluster2.json

localrules: frip_plot, consensus_peaks, merge_counts

rule all:
    input:
        # pre-process = align + make tracks + macs2 + homer
        expand("samples/bams/{sample}.sorted.bam", sample = SAMPLES),
        expand("samples/bigBed/{sample}.all.bb", sample = SAMPLES),
        expand("samples/bigwig/{sample}.bw", sample = SAMPLES),
        expand("samples/macs/{sample}/{sample}_peaks.narrowPeak", sample  = CASES),
        expand("results/motifs/{sample}/homerResults.html", sample = CASES),
        # pre-DESeq2 = consensus peaks + counts table
        "samples/macs/all_peaks.bed",
        expand("samples/macs/{factor}_peaks.bed", factor = FACTORS),
        expand("results/counts/{factor}_counts.txt", factor = FACTORS),
        # DESeq2 = pca + normalized counts + analyze contrasts
        expand("results/de/pca/{factor}.png", factor = FACTORS),
        expand(["results/de/{factor}/{factor}_deseq2_norm_counts.txt",
                "results/de/{factor}/{factor}_deseq2_lognorm_counts.txt"], factor = FACTORS),
        directory(expand("results/de/{factor}", factor = FACTORS)),
        # quality control = frip + fastqc-screen + dispersion of coverage + library complexity + duplicate rate
        expand("results/qc/plotHeatmap/{factor}_all_samples.png", factor = FACTORS),
        expand(["samples/qc/fastq_screen/{cases}/{cases}_{dir}_screen.{ext}",
                "samples/qc/fastq_screen/{controls}/{controls}_{dir}_screen.{ext}"],
                dir = ["R1", "R2"],
                cases = CASES,
                controls = CONTROLS_UNIQUE,
                ext = ["png", "txt", "html"]),
        expand("samples/qc/{factor}_fingerprint.png", factor = FACTORS),
        "results/qc/frip.html",

# downstream = chip_screen + chip_fisher
# chip screen - intersect consensus peaks with public chip data
# expand(["data/chip_screen/{db}/{cell_line}/consensus_peaks_{file}.bed.gz",
#         "data/chip_screen/{db}/{cell_line}/consensus_peaks_{file}.out"], 
#         zip, db = DB, cell_line = CELL_LINE, file = FILE)

include: "rules/align.smk"
include: "rules/peaks.smk"
include: "rules/qc.smk"
include: "rules/deseq2.smk"