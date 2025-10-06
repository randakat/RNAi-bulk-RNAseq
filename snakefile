import pandas as pd

# -----------------------------
# Define absolute paths to binaries
# -----------------------------
FASTP = "/home/bret/.snakemake-conda/fastp/bin/fastp"
STAR = "/home/bret/.snakemake-conda/star/bin/STAR"
FEATURECOUNTS = "/home/bret/.snakemake-conda/featurecounts/bin/featureCounts"
SUBREAD = "/home/bret/.snakemake-conda/subread/bin/featureCounts"


# Directories
TRIMMED_DIR = "trimmed"
ALIGNED_DIR = "aligned"
COUNTS_DIR = "counts"
COUNTS_FRAC_DIR = "counts_fractional"

# -----------------------------
# Samples
# -----------------------------
# Load sample table
sample_table = pd.read_csv("samples.tsv", sep="\t", header=None, names=["sample", "R1", "R2"])

SAMPLES = sample_table["sample"].tolist()

# Map sample names to FASTQ paths
R1 = dict(zip(sample_table["sample"], sample_table["R1"]))
R2 = dict(zip(sample_table["sample"], sample_table["R2"]))

# -----------------------------
# Rules
# -----------------------------

rule all:
    input:
        expand("counts/{sample}.counts.txt", sample=SAMPLES),
        expand("counts_fractional/{sample}.counts.txt", sample=SAMPLES)

# -----------------------------
# Fastp: trimming raw reads
# -----------------------------
rule fastp:
    input:
        R1=lambda wc: R1[wc.sample],
        R2=lambda wc: R2[wc.sample]
    output:
        R1_trim=f"{TRIMMED_DIR}/{{sample}}_R1.trimmed.fastq.gz",
        R2_trim=f"{TRIMMED_DIR}/{{sample}}_R2.trimmed.fastq.gz"
    threads: 16
    shell:
          """
        fastp -i {input.R1} -I {input.R2} \
              -o {output.R1_trim} -O {output.R2_trim} \
              -w {threads} \
              --detect_adapter_for_pe
        """

# -----------------------------
# STAR index: genome preparation
# -----------------------------
rule star_index:
    input:
        genome="genome.fa",
        gtf="annotation.gtf"
    output:
        "star_index"
    threads: 16
    params:
        sjdbOverhang=99
    shell:
        "{STAR} --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {output} --genomeFastaFiles {input.genome} "
        "--sjdbGTFfile {input.gtf} --sjdbOverhang {params.sjdbOverhang}"

# -----------------------------
# STAR alignment: mapping reads
# -----------------------------
rule star_align:
    input:
        R1_trim=f"{TRIMMED_DIR}/{{sample}}_R1.trimmed.fastq.gz",
        R2_trim=f"{TRIMMED_DIR}/{{sample}}_R2.trimmed.fastq.gz",
        index="star_index"
    output:
        bam=f"{ALIGNED_DIR}/{{sample}}.bam"
    threads: 16
    shell:
        """
        STAR --genomeDir {input.index} --readFilesIn {input.R1_trim} {input.R2_trim} \
             --runThreadN {threads} --outFileNamePrefix {ALIGNED_DIR}/{{wildcards.sample}}_ \
             --outSAMtype BAM SortedByCoordinate
        """


# -----------------------------
# FeatureCounts: integer counts
# -----------------------------
rule featurecounts:
    input:
        bam=f"{ALIGNED_DIR}/{{sample}}.bam",
        gtf="annotation.gtf"
    output:
        counts=f"{COUNTS_DIR}/{{sample}}.counts.txt"
    threads: 4
    shell:
        """
        featureCounts -T {threads} -p -s 0 \
                      -a {input.gtf} -o {output.counts} {input.bam}
        """


# -----------------------------
# FeatureCounts: fractional counts
# -----------------------------
rule featurecounts_fraction:
    input:
        bam="aligned/{sample}.bam",
        gtf="annotation.gtf"
    output:
        counts="counts_fractional/{sample}.counts.txt"
    threads: 4
    shell:
        "{SUBREAD} -T {threads} -p -s 0 -M --fraction -a {input.gtf} -o {output.counts} {input.bam}"

