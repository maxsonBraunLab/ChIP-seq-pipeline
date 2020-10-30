rule mapReads_paired_cases:
    input:
        R1 = "samples/cases/{sample}_R1.fastq.gz",
        R2 = "samples/cases/{sample}_R2.fastq.gz"
    output:
        "samples/bigBed/{sample}.all.bb",
        temp("samples/bams/{sample}.sorted.bam"),
        temp("samples/bams/{sample}.sorted.bam.bai")
    params:
        assembly = config["assembly"],
        filter = config["filter"],
        type = config["seq_type"]
    conda:
        "../envs/chip.yaml"
    threads: 12
    shell:
        """scripts/mapReads.sh -i {input.R1} -I {input.R2} -t {params.type} -n {wildcards.sample} -a {params.assembly} -f {params.filter} -u {threads} -o samples/bigBed/"""

rule mapReads_paired_controls:
    input:
        R1 = "samples/controls/{sample}_R1.fastq.gz",
        R2 = "samples/controls/{sample}_R2.fastq.gz"
    output:
        "samples/bigBed/{sample}.all.bb",
        temp("samples/bams/{sample}.sorted.bam"),
        temp("samples/bams/{sample}.sorted.bam.bai")
    params:
        assembly = config["assembly"],
        filter = config["filter"],
        type = config["seq_type"]
    conda:
        "../envs/chip.yaml"
    threads: 12
    shell:
        """scripts/mapReads.sh -i {input.R1} -I {input.R2} -t {params.type} -n {wildcards.sample} -a {params.assembly} -f {params.filter} -u {threads} -o samples/bigBed/"""

rule mapReads_single_cases:
    input:
        R1 = "samples/cases/{sample}.fastq.gz"
    output:
        "samples/bigBed/{sample}.all.bb",
        "samples/bams/{sample}.sorted.bam",
    params:
        assembly = config["assembly"],
        filter = config["filter"],
        type = config["seq_type"]
    conda:
        "../envs/chip.yaml"
    threads: 12
    shell:
        """scripts/mapReads.sh -i {input.R1} -t {params.type} -n {wildcards.sample} -a {params.assembly} -f {params.filter} -u {threads} -o samples/bigBed/"""

rule mapReads_single_controls:
    input:
        R1 = "samples/controls/{sample}.fastq.gz"
    output:
        "samples/bigBed/{sample}.all.bb",
        "samples/bams/{sample}.sorted.bam",
    params:
        assembly = config["assembly"],
        filter = config["filter"],
        type = config["seq_type"]
    conda:
        "../envs/chip.yaml"
    threads: 12
    shell:
        """scripts/mapReads.sh -i {input.R1} -t {params.type} -n {wildcards.sample} -a {params.assembly} -f {params.filter} -u {threads} -o samples/bigBed/"""

rule makeTracks:
    input:
        "samples/bigBed/{sample}.all.bb"
    output:
        "samples/bigwig/{sample}.bw"
    params:
        assembly = config["assembly"],
        ext = config["extension"]
    conda:
        "../envs/chip.yaml"
    shell:
        """scripts/makeTracks.sh -i {input} -o {output} -g {params.assembly} -e {params.ext}"""

# count total aligned reads
rule align_stats_1:
    input:
        "samples/bams/{sample}.sorted.bam"
    output:
        "samples/bams/stats/{sample}_tot_reads.txt"
    conda:
        "../envs/chip.yaml"
    shell:
        """
        tot_reads=$(samtools view -c {input})
        aligned_reads=$(samtools view -h {input} | awk '{{if ($3 != "*") {{print $0}} }}' | samtools view -c)
        echo -e "{wildcards.sample}\n$tot_reads\n$aligned_reads" > {output}
        """

rule markdup:
    input:
        "samples/bams/{sample}.sorted.bam"
    output:
        temp("samples/bams/{sample}.dedup.sorted.bam"),
        temp("samples/bams/{sample}.dedup.sorted.bam.bai")
    conda:
        "../envs/sambamba.yaml"
    log:
        "logs/markdup/{sample}.log"
    threads: 8
    shell:
        "sambamba markdup -r -t {threads} {input} {output}"

rule align_stats_2:
    input:
        rules.markdup.output[0]
    output:
        "samples/bams/stats/{sample}_uniq_reads.txt"
    conda:
        "../envs/chip.yaml"
    shell:
        """
        uniq_reads=$(samtools view -c {input})
        echo -e "{wildcards.sample}\n$uniq_reads" > {output}
        """

rule rm_unmapped:
    input:
        rules.markdup.output
    output:
        "samples/bams/{sample}.mapped.dedup.sorted.bam"
    conda:
        "../envs/chip.yaml"
    shell:
        """samtools view -h {input} | awk '{{  if ($3 != "*") {{print $0}}  }}' | samtools view -bS > {output}; samtools index {output}"""

# rule align_stats:
#     input:
#         tot = expand("samples/bams/stats/{sample}_tot_reads.txt", sample = SAMPLES),
#         uniq = expand("samples/bams/stats/{sample}_uniq_reads.txt", sample = SAMPLES)
#     output:
#         "samples/bams/align_stats.txt"
#     run:
#         