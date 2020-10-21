# custom scaling factor for bam --> bigwig --------------------------------------------
# rule all: expand("samples/norm_bigwig/{sample}.bw", sample = SAMPLES),

if config["assembly"] == "hg38":
    genome_size = 2913022398
elif config["assembly"] == "mm10":
    genome_size = 2652783500
else:
    print("assembly must be hg38 or mm10. Your input: {}".format(config["assembly"]))

rule count_bams:
    input:
        "samples/bams/{sample}.sorted.bam"
    output:
        "samples/norm_bigwig/counts/{sample}.txt"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """ counts=$(samtools view -@ 4 -c {input}); echo -e "{wildcards.sample}\n$counts" > {output} """

rule gather_counts:
    input:
        expand("samples/norm_bigwig/counts/{sample}.txt", sample = SAMPLES)
    output:
        "samples/norm_bigwig/scaling_factors.txt"
    run:
        df = pd.concat([ pd.read_csv(i) for i in sorted(input) ], axis = 1)
        df = df.transpose()
        df = df.rename(columns={0: "total_reads"})
        df['scaling_factors'] = float(df.min()) / df.total_reads
        df.to_csv(str(output), sep = "\t")

rule normTracks:
    input:
        bam = "samples/bams/{sample}.sorted.bam",
        scaling_factors = rules.gather_counts.output
    output:
        "samples/norm_bigwig/{sample}.bw"
    params:
        g = genome_size
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        factor=$(grep {wildcards.sample} {input.scaling_factors} | cut -f3)
        bamCoverage -b {input.bam} -o {output} -p {threads} --scaleFactor $factor --binSize 10 --smoothLength 40 --effectiveGenomeSize {params.g} -r chr12
        """