#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(Biostrings)

#TSS_long_fa <- paste0(snakemake@config[["experiment_folder"]], "data/annotations/", "TSS_2kb.fa")
TSS_long_fa <- "../../data/annotations/TSS_2kb.fa"
TSS_short_fa <- "../../data/annotations/TSS_900bp.fa"
TSS_short_bed <- "../../data/annotations/TSS_900bp.bed"

############ Using the 2kb promoters for distribution plot around TSS as in Saxonov et al., PNAS (2006)

## Reading in a DNAStringSet the fasta file containing the 2kb sequences centered on transcripts TSS 
TSS_2kb_fa <- readDNAStringSet(TSS_long_fa, 
              format="fasta",
              nrec=-1L, skip=0L, seek.first.rec=FALSE,
              use.names=TRUE, with.qualities=FALSE)
## Removing those sequences which are shorter than 2kb because at the end of chromosomes (N=3)
TSS_2kb_fa <- TSS_2kb_fa[width(TSS_2kb_fa) == 2000]

## Function to create a list of IRanges of overlapping 50 nt ranges spanning the 2kb
overlapping_ranges <- function (i) {
                        if (i != 1961) {
                          my_start = i
                          my_width = 50
                        } else {
                          my_start = i
                          my_width = 2000 - 1961 + 1
                        }
                        return(IRanges(start = my_start, width = my_width))
}
overlapping_ranges_2kb_by50 <- lapply(seq(from = 1, to = 1961, by = 49), FUN = overlapping_ranges)

## Function to extract the overlapping ranges of sequence from the DNAStringSet. The Biostring function 'extractAt' would create a DNAStringSetList, which becomes a DNAStringSet (with the same informational content) with the function 'unlist' creating a list
extract_overlapping_ranges <- function (one_range, fasta_DNASS) {
  unlist(extractAt(fasta_DNASS, one_range))
}
## Applying the above function to the DNAStringSet creates a list of DNAStringSets
split_2kb <- lapply(overlapping_ranges_2kb_by50, FUN = extract_overlapping_ranges, fasta_DNASS = TSS_2kb_fa)

## Function to compute the normalized CpG fraction as described in Saxonov et al., PNAS (2006) and Weber et al., Nature Genetics (2006): numCpG divided by expected number of CpG (i.e. (number of Cs * number of Gs)/number of nucleotides).
compute_normCpG <- function (fasta_DNASS) {
  CpG_df <- data.frame(C = alphabetFrequency(fasta_DNASS)[, "C"], G = alphabetFrequency(fasta_DNASS)[, "G"], N = alphabetFrequency(fasta_DNASS)[, "N"])
  CpG_df$expCpG <- apply(CpG_df[, 1:2], 1, prod)/unique(width(fasta_DNASS))
  CpG_df$obsCpG <- dinucleotideFrequency(fasta_DNASS, step = 1)[,"CG"]
  CpG_df$obsExp <- CpG_df$obsCpG/CpG_df$expCpG
  # Writing 0 for those sequences that contain more than 10% of Ns
  CpG_df$obsExp[CpG_df[, "N"] > (unique(width(fasta_DNASS))/10)] <- 0
  return(CpG_df$obsExp)
}
## Applying the above function to the list of DNAStringSets creates a list of vectors, each vector being long as the number of transcripts in the input file (minus the non-2kb ones) and containing the normalized CpG fraction for one 50 nt subsequence
normCpG_list <- lapply(split_2kb, FUN = compute_normCpG)
## Transforming the list in a dataframe and naming the column from 1 to the last possible range
normCpG_df <- as.data.frame(normCpG_list, row.names = names(TSS_2kb_fa))
names(normCpG_df) <- c(1:41)

#### Plotting

### Distribution of CpG around TSS for all TSSs
#pdf(paste0(snakemake@config[["experiment_folder"]], snakemake@config[["analysis_output_dir"]], "CpG_analysis/CpG_distrib_around_TSS.pdf"), width = 10, height = 10)
pdf("../../analysis/CpG_distrib_around_TSS.pdf", width = 10, height = 10)
ggplot(data = data.frame(x = factor(1:41), y = colMeans(normCpG_df, na.rm = TRUE)), aes(x = x, y = y)) +
  geom_point() +
  scale_x_discrete(name = "\ndistance from TSS", breaks = c("1", "21", "41"), labels = c("-1000", "0", "1000")) +
  ylab("normalized CpG\n") +
  ggtitle("Distribution of CpGs around TSS\n") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()

############ Using the 900bp promoters for saving actual value of normalized CpG fraction and CpG percentage, as in Weber et al., Nature Genetics (2006): maximum value of normalized CpG fraction of out of all values computed for 500 bp intervals drwan in the 900bp regions with 5bps of offset between them 

## Reading in a DNAStringSet the fasta file containing the 2kb sequences centered on transcripts TSS 
TSS_900_fa <- readDNAStringSet(TSS_short_fa, 
                               format="fasta",
                               nrec=-1L, skip=0L, seek.first.rec=FALSE,
                               use.names=TRUE, with.qualities=FALSE)
## Removing those sequences which are shorter than 900bp (N=1)
TSS_900_fa <- TSS_900_fa[width(TSS_900_fa) == 900]

## Function to create a list of IRanges of overlapping 500 nt ranges with 5bp offset spanning the 900bp
overlapping_ranges_2 <- function (i) {
  if (i != 401) {
    my_start = i
    my_width = 500
  } else {
    my_start = i
    my_width = 900 - 401 + 1
  }
  return(IRanges(start = my_start, width = my_width))
}
overlapping_500_ranges_900_by5 <- lapply(seq(from = 1, to = 401, by = 5), FUN = overlapping_ranges_2)

## Applying the extract_overlapping_ranges function defined in first part of this script to the DNAStringSet creating a list of DNAStringSets
split_900bp <- lapply(overlapping_500_ranges_900_by5, FUN = extract_overlapping_ranges, fasta_DNASS = TSS_900_fa)

## Applying the compute_normCpG function defined in first part of this script to the list of DNAStringSets creating a list of vectors, each vector being long as the number of transcripts in the input file and containing the normalized CpG fraction for one 500 nt sequence
normCpG_list_900 <- lapply(split_900bp, FUN = compute_normCpG)
## Transforming the list in a dataframe and naming the column from 1 to the last possible range
normCpG_df_900 <- as.data.frame(normCpG_list_900, row.names = names(TSS_900_fa))
names(normCpG_df_900) <- c(1:81)

## Function to compute the CpG percentage
compute_perCpG <- function (fasta_DNASS) {
  CpG_df <- data.frame(C = alphabetFrequency(fasta_DNASS)[, "C"], G = alphabetFrequency(fasta_DNASS)[, "G"], N = alphabetFrequency(fasta_DNASS)[, "N"])
  CpG_df$obs_CpG <- dinucleotideFrequency(fasta_DNASS, step = 1)[,"CG"]
  # Writing NA for those sequences that contain more than 10% of Ns
  CpG_df$obs_CpG[CpG_df[, "N"] > (unique(width(fasta_DNASS))/10)] <- 0
  return(CpG_df$obs_CpG/unique(width(fasta_DNASS))*100)
}
## Applying the function to the list of DNAStringSets creating a dataframe, each row being a transcript in the input file and each column being the CpG percentage for one 500 nt subsequence
obsCpG_df_900 <- as.data.frame(sapply(split_900bp, FUN = compute_perCpG))
rownames(obsCpG_df_900) <-rownames(normCpG_df_900)
colnames(obsCpG_df_900) <- c(1:81)

## For each transcript, computing the max of norm CpG fraction and the mean of CpG percentage among all 500bp subsequences 
## The names given to the two computed parameters are the same used in the UCSC table browser CGIs track's table
obsCpG_df_900$perCpg <- apply(obsCpG_df_900, 1, mean, na.rm = TRUE)
normCpG_df_900$obsExp <- apply(normCpG_df_900, 1, max)
obsCpG_df_900$transcript_id <- gsub("\\..*", "", rownames(obsCpG_df_900))
normCpG_df_900$transcript_id <- gsub("\\..*", "", rownames(normCpG_df_900))

## Merging the 2 dataframes obtained with 900 bp promoter sequences
promoters_CpG_df <- merge(obsCpG_df_900[, 82:83], normCpG_df_900[, 82:83], by = "transcript_id")

## Writing the dataframe to a file which can be used as a database from which CpG information can be fetched for specific transcripts
#write.table(promoters_CpG_df, 
#            file = paste0(snakemake@config[["experiment_folder"]], snakemake@config[["analysis_output_dir"]], "CpG_analysis/promoters_CpG_df.txt"), 
#            sep = "\t",
#            quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(promoters_CpG_df, 
            file = "../../analysis/CpG_analysis/promoters_CpG_df.txt", 
            sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

## Relationship between percentage of CpG and normalized CpG fraction of promoters
pdf("../../analysis/CpG_analysis/perCpG_vs_normCpGvalue.pdf", width = 10, height = 8)
ggplot(data = promoters_CpG_df, aes(x = perCpg, y = obsExp)) +
  geom_point()
dev.off()

######## Separating promoters in two classes based on the normalized CpG fraction at their promoters

### Distribution of normalized CpG values among TSSs
#pdf(paste0(snakemake@config[["experiment_folder"]], snakemake@config[["analysis_output_dir"]], "CpG_analysis/CpG_distrib_among_TSSs_hist.pdf"), width = 10, height = 8)
pdf("../../analysis/CpG_analysis/CpG_distrib_among_TSSs_hist.pdf", width = 10, height = 8)
par(mar=c(5.1,5.1,4.1,2.1))
normCpG_hist_900 <- hist(apply(normCpG_df_900[, 1:(ncol(normCpG_df_900)-1)], 1, max, na.rm = TRUE), 
                         breaks = 100, 
                         xlab = "normalized CpG", 
                         ylab = "number of TSS", 
                         main = "Distribution of normalized CpG values among TSSs",
                         cex.lab=1.8, cex.axis=1.5, cex.main=1.8)

## Fitting the histogram to two gaussians and adding them to the plot
hist_to_fit <- data.frame(x = normCpG_hist_900$mids, y = normCpG_hist_900$counts)

nls_fit <- nls(y ~ (a/b)*exp(-(x-c)^2/(2*b^2)) + (d/e)*exp(-(x-f)^2/(2*e^2)),
        data=hist_to_fit,
        start=list(a=(1/sqrt(2*pi)) / 0.12, 
        b=0.12, 
        c=0.28,
        d=(1/sqrt(2*pi)) / 0.14, 
        e=0.14, 
        f=0.99),
        control=nls.control(tol=1E-5, minFactor=1/1024),
        trace=TRUE,
        algorithm="port")

lines(normCpG_hist_900$mids, predict(nls_fit), col = "red")
dev.off()

### Finding the intersection between the two gaussians
## Building two gaussians based on the mean and sd parameters found by nls fitting and computing the difference function
f <- function(x) dnorm(x, m=coef(nls_fit)["c"], sd=coef(nls_fit)["b"]) - dnorm(x, m=coef(nls_fit)["f"], sd=coef(nls_fit)["e"])
## Finding the intersection point. The ínterval'argument specifies the interval which should include the point
gaussians_intersection_point <- uniroot(f, interval=c(0.3, 0.7))$root

## Separating promoters in two classes based on the intersection between the two gaussians
HCG_transcripts <- gsub("\\..*", "", promoters_CpG_df$transcript_id[promoters_CpG_df$obsExp > gaussians_intersection_point])
LCG_transcripts <- gsub("\\..*", "", promoters_CpG_df$transcript_id[promoters_CpG_df$obsExp < gaussians_intersection_point])
## Separating promoters in two more extreme classes
HCG_transcripts_e <- gsub("\\..*", "", promoters_CpG_df$transcript_id[promoters_CpG_df$obsExp > gaussians_intersection_point + 0.1])
LCG_transcripts_e <- gsub("\\..*", "", promoters_CpG_df$transcript_id[promoters_CpG_df$obsExp < gaussians_intersection_point - 0.1])

### Writing HGC and LCG transcripts to txt files
#write.table(HCG_transcripts, 
#            file = paste0(snakemake@config[["experiment_folder"]], snakemake@config[["analysis_output_dir"]], "CpG_analysis/HCG_transcripts.txt"), 
#            quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(LCG_transcripts, 
#            file = paste0(snakemake@config[["experiment_folder"]], snakemake@config[["analysis_output_dir"]], "CpG_analysis/LCG_transcripts.txt"), 
#            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(HCG_transcripts, 
            file = "../../analysis/CpG_analysis/HCG_transcripts.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(LCG_transcripts, 
            file = "../../analysis/CpG_analysis/LCG_transcripts.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

### Distribution of CpG around TSS, separating HCG and LCG promoters

CpG_around_TSSs <- data.frame(x = factor(1:41), HCG = colMeans(normCpG_df[which(gsub("\\..*", "", rownames(normCpG_df)) %in% HCG_transcripts), ], na.rm = TRUE), LCG = colMeans(normCpG_df[which(gsub("\\..*", "", rownames(normCpG_df)) %in% LCG_transcripts), ], na.rm = TRUE))

#pdf(paste0(snakemake@config[["experiment_folder"]], snakemake@config[["analysis_output_dir"]], "CpG_analysis/CpG_distrib_around_TSS_HCGvsLCG.pdf"), width = 10, height = 8)
pdf("../../analysis/CpG_analysis/CpG_distrib_around_TSS_HCGvsLCG.pdf", width = 10, height = 8)
ggplot(data = melt(CpG_around_TSSs), aes(x = x, y = value)) +
  geom_point(aes(colour = variable)) +
  scale_x_discrete(name = "\ndistance from TSS", breaks = c("1", "21", "41"), labels = c("-1000", "0", "1000")) +
  ylab("normalized CpG\n") +
  ggtitle("Distribution of CpGs around TSS\n") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.margin=unit(c(1,1,1,1),"cm")) +
  guides(colour = guide_legend(title="Type of TSS"))
dev.off()

### Gene Ontology over-representation test for HCG and LCG promoters

library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)

## Function which: 1) perform the GO Over-representation Test for Biological Process, 2) filter for a level of GO terms, 3) simplify the result by grouping those terms which have correlation > 0.7 (keeping only the one with the minimum q-value among the GO terms of one group) 
perform_BP_go_ora <- function (my_GOIs, f_level) {
  BP <- enrichGO(gene = my_GOIs,
                 universe = keys(org.Mm.eg.db, keytype="ENSEMBL"),
                 OrgDb = org.Mm.eg.db, 
                 keyType = "ENSEMBL",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.1,
                 readable = TRUE)
  BP_f <- gofilter(BP, level = , f_level)
  return(BP_f)
  return(simplify(BP_f, cutoff=0.7, by="qvalue", select_fun=min))
}

# Creating the tx2gene IDs df; the tx2gene tsv file has been created with custom shell script tx2gene_from_encode_gtf.sh
gene_name_map <- fread("../../data/annotations/gencode.vM25.annotation.tx2gene.tsv", col.names = c("ensembl_transcript_id", "ensembl_gene_id"))

# tx2gene IDs for the two extreme classes of promoters
HCG_genes <- unique(sort(gene_name_map[gsub("\\..*", "", ensembl_transcript_id) %in% HCG_transcripts_e, gsub("\\..*", "", ensembl_gene_id)]))
LCG_genes <- unique(sort(gene_name_map[gsub("\\..*", "", ensembl_transcript_id) %in% LCG_transcripts_e, gsub("\\..*", "", ensembl_gene_id)]))

# GO test and transforming GO object to data.frame
GO_HCG <- fortify(perform_BP_go_ora(my_GOIs = HCG_genes, f_level = 4), showCategory = 10)
GO_HCG$Description <- factor(GO_HCG$Description, levels=rev(c(as.character(fortify(GO_HCG, showCategory = 10)$Description))))
GO_LCG <- fortify(perform_BP_go_ora(my_GOIs = LCG_genes, f_level = 4), showCategory = 10)
GO_LCG$Description <- factor(GO_LCG$Description, levels=rev(c(as.character(fortify(GO_LCG, showCategory = 10)$Description))))

# barplot of enriched terms
pdf("../../analysis/CpG_analysis/GO_HCG_LCG_genes.pdf", width = 10, height = 8)
ggplot(GO_HCG, aes_string(x = "Count", y = "Description", fill = "qvalue")) +
  scale_fill_continuous(low="red", high="blue", name = "qvalue", guide = guide_colorbar(reverse=TRUE)) +
  geom_col() +
  xlab("Gene count") + ylab(NULL) +
  theme_dose(13) +
  ggtitle("Genes with High CpG promoters") +
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))
ggplot(GO_LCG, aes_string(x = "Count", y = "Description", fill = "qvalue")) +
  scale_fill_continuous(low="red", high="blue", name = "qvalue", guide = guide_colorbar(reverse=TRUE)) +
  geom_col() +
  xlab("Gene count") + ylab(NULL) +
  ggtitle("Genes with Low CpG promoters") +
  theme_dose(13) +
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))
dev.off()

