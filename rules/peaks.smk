def get_control(wildcards):
    return md.loc[(wildcards.sample),["Control"]]

# call peaks ------------------------------------------------------------

rule bb2bed:
    input:
        "samples/bigBed/{sample}.all.bb"
    output:
        "samples/bed/{sample}.bed"
    params:
        bigBedToBed = config["bed_tool"]
    shell:
        "{params.bigBedToBed} {input} {output}"

rule call_peaks:
    input:
        case = "samples/bed/{sample}.bed",
        control = get_control
    output:
        "results/motifs/{sample}/homerResults.html",
        "samples/macs/{sample}/{sample}_peaks.narrowPeak"
    params:
        assembly = config["assembly"],
        peak = get_peak
    conda:
        "../envs/peaks.yaml"
    shell:
        """scripts/callPeaks.sh -t {input.case} -c {input.control} -n {wildcards.sample} -p {params.peak} -a {params.assembly}"""

# consensus peaks + counts ----------------------------------------------------

# output 2 files. 1 is all peaks concatenated together, 2 is consensus peaks in n_intersects number of replicates. 2 depends on 1 for finding consensus peaks.
rule consensus_peaks:
    input:
        expand("samples/macs/{sample}/{sample}_peaks.narrowPeak", sample = CASES)
    output:
        all_peaks = "samples/macs/all_peaks.bed",
        consensus_per_factor = expand("samples/macs/{factor}_peaks.bed", factor = FACTORS),
        consensus_stats = "samples/macs/consensus_stats.txt"
    params:
        present_in_number = config["n_intersects"],
        blacklist = config["blacklist"],
        genome = config["assembly"],
        script_file = "scripts/consensus_per_factor.R"
    conda:
        "../envs/consensus_peaks.yaml"
    log:
        "logs/consensus_peaks/{}.log".format(timestamp)
    threads: 8
    shell:
        """
        cat {input} > {output.all_peaks}
        Rscript --vanilla {params.script_file} {output.all_peaks} {params.blacklist} {params.present_in_number} {params.genome} {threads}
        """

rule sample_counts:
    input:
        bams = "samples/bams/{sample}_{factor}.mapped.dedup.sorted.bam",
        peaks = "samples/macs/{factor}_peaks.bed"
    output:
        "samples/macs/counts/{sample}_{factor}.txt"
    params:
        header = '\t'.join(["Chr","start","stop","V4","V5","V6"]) + '\t' + "{sample}_{factor}"
    conda:
        "../envs/chip.yaml"
    shell:
        "bedtools multicov -bams {input.bams} -bed {input.peaks} > {output}; sed -i '1i {params.header}' {output}"

rule merge_counts:
    input:
        expand("samples/macs/counts/{sample}_{{factor}}.txt", sample = COND_REPS)
    output:
        "results/counts/{factor}_counts.txt"
    run:
        dataframes = []
        for file in input:
            dataframes.append(pd.read_csv(file, sep='\t' ))
        merged = pd.concat(dataframes, axis=1)
        merged = merged.loc[:,~merged.columns.duplicated()]
        merged.to_csv(str(output), header=True, index=False, sep='\t')