# RNAi-bulk-RNAseq

# Step 1: create envs for each package
mamba env create -f envs/fastp.yaml -p /home/bret/.snakemake-conda/fastp
mamba env create -f envs/star.yaml -p /home/bret/.snakemake-conda/star
mamba env create -f envs/subread.yaml -p /home/bret/.snakemake-conda/subread
mamba env create -f envs/featurecounts.yaml -p /home/bret/.snakemake-conda/featurecounts

# Step 2: run snakefile
snakemake --cores 24 --printshellcmds
