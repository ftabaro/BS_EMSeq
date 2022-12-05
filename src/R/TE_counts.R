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
"CpG methylation counted on  nsposable Elements for following files:"
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

TE <- as(fread(paste0(snakemake@config[["experiment_folder"]], snakemake@input[["TE_annotation"]]), col.names = c("seqnames", "start", "end", "repName", "length", "strand", "repClass/repFamily")), "GRanges")
# adding +1 to start is necessary having used a simple fread function to read the bed file
start(TE) <- start(TE) + 1
# removing 'chr' from chromosome names
seqlevels(TE) <- substr(seqlevels(TE), start = 4, stop = 100)

# reading bed file containing coverage outliers (10x above median, found by quantifying over 25k genome windows)
outliers_25k <- as(import(paste0(snakemake@config[["experiment_folder"]], snakemake@config[["analysis_output_dir"]], "SeqMonk_project/25k_outliers_10_above_median.bed")), "GRanges")
# removing 'chr' from chromosome names
seqlevels(outliers_25k) <- substr(seqlevels(outliers_25k), start = 4, stop = 100)

# filtering out the outliers regions from TE
TE_f <- TE[-queryHits(findOverlaps(TE, outliers_25k))]

### Counting methylation over TE
mrl_TE <- regionCounts(object = mrlDB, regions = TE_f, chunk.size = 1e6, cov.bases = 10, save.db = FALSE)

# saving the R objects
save(mrl_TE, file = paste0(Rdata_dir, "R/mrl_TE.Rdata"))
