# Analysis of EMSeq data for project Boulard - Stork

* Until count of methylation on tiles, promoters and CpG Islands, analysis is performed using snakemake pipeline in snake-make/Snakefile.
* Snakefile is run on the EMBL cluster using pipeline_wrapper.sh in src/sh
* Snakefile uses conda environments in env/conda and singularity containers built with recipes in singularity/recipes
* config file for Snakefile and for SLURM are in config/
* important for reproducibility: 
  - conda version 4.8.3
  - singularity version used to build containers is 3.5.3
  - snakemake version to run Snakefile is 5.9.1
* Using the Rdata files produced in Snakefile, I then continue the analysis in the Rmd files, which can run both on the cluster and locally on a personal computer (having at least 16Gb of RAM).

## Steps outside of the Snakemake pipeline and the Rmd files

* Definition of promoters as HCG or LCG and creation of 'HCG_transcripts.txt' and 'LCG_transcripts.txt' files is done in R script 'src/R/TSSs_CpG.R'.
* The bed file containing genome tiles of 200 and 100 consecutive CpGs which are used in the 'CpG_tiles_counts.R' script (run in a step of the Snakemake pipeline) are produced using [SeqMonk](https://www.bioinformatics.babraham.ac.uk/projects/seqmonk/) Read Position Probe Generator, asking for a minimum of 1 read count per cytocine over all sample.
* The bed file containing coverage outliers used in the 'CpG_tiles_counts.R' script were also produced using SeqMonk, by quantifying methylation over 25 kb windows and picking those regions with methylation level 10x above median.
