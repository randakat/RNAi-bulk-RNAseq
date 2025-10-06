SAMPLES = ["zfp3_1", "zfp3_2", "zfp3_3", "ctrl_1", "ctrl_2", "ctrl_3"]

# Paths to pre-created conda envs
FASTP_ENV = "/home/bret/.snakemake-conda/fastp"
STAR_ENV = "/home/bret/.snakemake-conda/star"
SUBREAD_ENV = "/home/bret/.snakemake-conda/subread"
FEATURECOUNTS_ENV = "/home/bret/.snakemake-conda/featurecounts"

# all
rule all:
    input:
        expand("trimmed/{sample}_R1.trimmed.fastq.gz", sample=SAMPLES),
        expand("trimmed/{sample}_R2.trimmed.fastq.gz", sample=SAMPLES),
        expand("aligned/{sample}.bam", sample=SAMPLES),
        expand("counts/{sample}.counts.txt", sample=SAMPLES),
        expand("counts_fractional/{sample}.counts.txt", sample=SAMPLES),
        "merged_counts.txt"

# FASTQ trimming
rule fastp:
    input:
        R1=lambda wildcards: f"/media/bret/Disk5/KR_RNA-seq/reads/{wildcards.sample}_R1_001.fastq.gz",
        R2=lambda wildcards: f"/media/bret/Disk5/KR_RNA-seq/reads/{wildcards.sample}_R2_001.fastq.gz"
    output:
        R1_trimmed="trimmed/{sample}_R1.trimmed.fastq.gz",
        R2_trimmed="trimmed/{sample}_R2.trimmed.fastq.gz"
    threads: 4
    conda:
        FASTP_ENV
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -o {output.R1_trimmed} -O {output.R2_trimmed} -w {threads}
        """

# genome index
rule star_index:
    input:
        genome="genome.fa",
        gtf="annotation.gtf"
    output:
        directory("star_index")
    threads: 16
    conda:
        STAR_ENV
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.genome} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang 99
        """

# STAR alignment
rule star_align:
    input:
        R1="trimmed/{sample}_R1.trimmed.fastq.gz",
        R2="trimmed/{sample}_R2.trimmed.fastq.gz",
        index="star_index"
    output:
        bam="aligned/{sample}.bam"
    threads: 16
    conda:
        STAR_ENV
    shell:
        """
        STAR --genomeDir {input.index} \
             --readFilesIn {input.R1} {input.R2} \
             --readFilesCommand zcat \
             --runThreadN {threads} \
             --outFileNamePrefix aligned/{wildcards.sample}_ \
             --outSAMtype BAM SortedByCoordinate
        mv aligned/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
        """

# featureCounts
rule featurecounts:
    input:
        bam="aligned/{sample}.bam",
        gtf="annotation.gtf"
    output:
        "counts/{sample}.counts.txt"
    threads: 4
    conda:
        FEATURECOUNTS_ENV
    shell:
        """
        featureCounts -T {threads} -p -s 0 \
                      -a {input.gtf} -o {output} {input.bam}
        """

# featureCounts fractional multimapping
rule featurecounts_fraction:
    input:
        bam="aligned/{sample}.bam",
        gtf="annotation.gtf"
    output:
        "counts_fractional/{sample}.counts.txt"
    threads: 4
    conda:
        FEATURECOUNTS_ENV
    shell:
        """
        featureCounts -T {threads} -p -s 0 -M --fraction \
                      -a {input.gtf} -o {output} {input.bam}
        """

# Merge counts
rule merge_counts:
    input:
        expand("counts/{sample}.counts.txt", sample=SAMPLES)
    output:
        "merged_counts.txt"
    conda:
        SUBREAD_ENV
    shell:
        """
        paste {input} > {output}
        """
