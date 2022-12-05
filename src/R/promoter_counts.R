library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(methylKit)

# directory containing bismark meth extractor coverage output files
samples_dir = snakemake@input[["meth_cov_dir"]]

# directory to store the tabix files, used to index the compressed database that will be created from the coverage files
tbi_dir = paste0(snakemake@config[["experiment_folder"]], snakemake@config[["analysis_output_dir"]], "methylDB")

# directory to store the Rdata files I want to save
Rdata_dir = paste0(snakemake@config[["experiment_folder"]], snakemake@config[["script_dir"]])

# names of files without path, sorted alphabetically
files <- sort(grep(".*bismark.cov.gz", list.files(path = samples_dir), value = TRUE))

# names of files with path
files_w_path <- file.path(samples_dir, files)
"CpG methylation counted in promoters for following files:"
files_w_path

# loading samples table
samples_table <- fread(snakemake@config[["input_samples_table"]])

# adding samples' short names
samples_table$sample_name <- gsub(".*EMSeq_", "", samples_table$sample)
# sort the table by sample short names
setkey(samples_table, sample_name)

# read the files to a methylRawList DB object
#mrlDB <- methRead(as.list(files_w_path),
#                  sample.id=as.list(samples_table$sample_name),
#                  assembly="mm10",
#                  treatment=as.numeric(ifelse(samples_table$condition == "KO", 1, 0)),
#                  pipeline="bismarkCoverage",
#                  context="CpG",
#                  mincov=1,
#                  dbtype = "tabix",
#                  dbdir = tbi_dir)

# saving the R object
#save(mrlDB, file = paste0(Rdata_dir, "R/mrlDB.Rdata"))
load(paste0(Rdata_dir, "R/mrlDB.Rdata"))

# reading the bed file containing N bp regions centered on TSS as GRanges object i.e. promoters
TSS <- as(import(paste0(snakemake@config[["experiment_folder"]], snakemake@input[["promoters_bed"]])), "GRanges")
# removing 'chr' from chromosome names
seqlevels(TSS) <- substr(seqlevels(TSS), start = 4, stop = 100)
seqlevels(TSS)[seqlevels(TSS) == "M"] <- "MT"

# reading bed file containing coverage outliers (10x above median, found by quantifying over 25k genome windows)
outliers_25k <- as(import(paste0(snakemake@config[["experiment_folder"]], snakemake@config[["analysis_output_dir"]], "SeqMonk_project/25k_outliers_10_above_median.bed")), "GRanges")
# removing 'chr' from chromosome names
seqlevels(outliers_25k) <- substr(seqlevels(outliers_25k), start = 4, stop = 100)

# filtering out the outliers regions from the promoters
TSS_f <- TSS[-queryHits(findOverlaps(TSS, outliers_25k))]

### Counting methylation over promoters
#mrlDB_TSS <- regionCounts(object = mrlDB, regions = TSS_f, chunk.size = 1e6, cov.bases = 10, save.db = TRUE, suffix = "TSS", dbdir = tbi_dir)

# saving the R objects
#save(mrlDB_TSS, file = paste0(Rdata_dir, "R/mrlDB_TSS.Rdata"))

# transforming the DB object in in-memory ones for convenience
#mrl_TSS <- as(mrlDB_TSS,"methylRawList")

# saving the R objects
#save(mrl_TSS, file = paste0(Rdata_dir, "R/mrl_TSS.Rdata"))
load(paste0(Rdata_dir, "R/mrl_TSS.Rdata"))

### Finding Differentially Methylated Promoters (DMPs)

## Reorganizing the methylation object
# based on genotype
## function to reorganize by genotype
reorganize_by_geno <- function (my_geno, my_tps, mrl_obj) {
  mrl_obj_s <- reorganize(mrl_obj, sample.ids = samples_table[condition == my_geno & group_or_time_point
%in% my_tps, sample_name], treatment = ifelse(samples_table[condition == my_geno & group_or_time_point
%in% my_tps, group_or_time_point] == min(samples_table[condition == my_geno & group_or_time_point
%in% my_tps, group_or_time_point]), 0, 1))
  return(methylKit::unite(mrl_obj_s))
}
meth_by_geno_0vs7d <- lapply(X = unique(samples_table$condition), FUN = reorganize_by_geno, mrl_obj = mrl_TSS, my_tps = c(1,2))
meth_by_geno_7vs14d <- lapply(X = unique(samples_table$condition), FUN = reorganize_by_geno, mrl_obj = mrl_TSS, my_tps = c(2,3))
meth_by_geno_0vs14d <- lapply(X = unique(samples_table$condition), FUN = reorganize_by_geno, mrl_obj = mrl_TSS, my_tps = c(1,3))

# based on time point
## function to reorganize by time point
reorganize_by_tp <- function (my_tp, mrl_obj) {
  mrl_obj_s <- reorganize(mrl_obj, sample.ids = samples_table[group_or_time_point == my_tp, sample_name], treatment = ifelse(samples_table[group_or_time_point == my_tp, condition] == "WT", 0, 1))
  return(methylKit::unite(mrl_obj_s))
}
meth_by_tp <- lapply(X = unique(samples_table$group_or_time_point), FUN = reorganize_by_tp, mrl_obj = mrl_TSS)

## computing differential methylation test
## logistic regression based model with over-dispersion correction and Chi-square test
# between genotypes for each time point
dm_geno_for_each_tp <- lapply(meth_by_tp, FUN = function (tp) {calculateDiffMeth(tp, overdispersion = "MN", test ="Chisq", mc.cores = 5)})
# between time points for each genotype
dm_tp_for_each_geno_0vs7d <- lapply(meth_by_geno_0vs7d, FUN = function (tp) {calculateDiffMeth(tp, overdispersion = "MN", test ="Chisq", mc.cores = 5)})
dm_tp_for_each_geno_7vs14d <- lapply(meth_by_geno_7vs14d, FUN = function (tp) {calculateDiffMeth(tp, overdispersion = "MN", test ="Chisq", mc.cores = 5)})
dm_tp_for_each_geno_0vs14d <- lapply(meth_by_geno_0vs14d, FUN = function (tp) {calculateDiffMeth(tp, overdispersion = "MN", test ="Chisq", mc.cores = 5)})

# save the objects to Rdata files
save(dm_geno_for_each_tp, file = paste0(Rdata_dir, "R/dm_prom_geno_for_each_tp.Rdata"))
save(dm_tp_for_each_geno_0vs7d, file = paste0(Rdata_dir, "R/dm_prom_tp_for_each_geno_0vs7d.Rdata"))
save(dm_tp_for_each_geno_7vs14d, file = paste0(Rdata_dir, "R/dm_prom_tp_for_each_geno_7vs14d.Rdata"))
save(dm_tp_for_each_geno_0vs14d, file = paste0(Rdata_dir, "R/dm_prom_tp_for_each_geno_0vs14d.Rdata"))

