#!/bin/bash

#- take TSS from mm10 gtf file, keeping the gene information
#- for both 900bp and 2kb wide sequences around TSS:
#   - get bed and fasta files
#   - compute obsExp and perCpg values as in Weber et al. Nature Genetics (2006)
#   - creates bed files of unique coordinates for methylation count
#- create bed file of non-promoter CGIs

gene_annotation=${1}
chr_sizes=${2}
genome_fasta=${3}
CGI_annotation=${4}

bedtools flank -s -i ${gene_annotation} -g ${chr_sizes} -l 1000 -r 0 | bedtools slop -s -i - -g ${chr_sizes} -l 0 -r 1000 | bedtools getfasta -s -bedOut -name -fi ${genome_fasta} -bed - > data/annotations/TSS_2kb.bed 
bedtools flank -s -i ${gene_annotation} -g ${chr_sizes} -l 1000 -r 0 | bedtools slop -s -i - -g ${chr_sizes} -l 0 -r 1000 | bedtools getfasta -s -name -fi ${genome_fasta} -bed - > data/annotations/TSS_2kb.fa 

bedtools flank -s -i ${gene_annotation} -g ${chr_sizes} -l 700 -r 0 | bedtools slop -s -i - -g ${chr_sizes} -l 0 -r 200 | bedtools getfasta -s -bedOut -name -fi ${genome_fasta} -bed - > data/annotations/TSS_900bp.bed 
bedtools flank -s -i ${gene_annotation} -g ${chr_sizes} -l 700 -r 0 | bedtools slop -s -i - -g ${chr_sizes} -l 0 -r 200 | bedtools getfasta -s -name -fi ${genome_fasta} -bed - > data/annotations/TSS_900bp.fa 
# Some times, the same 2kb sequence around TSS belongs to more than one annotated transcript, therefore this would create duplicated rows when counting methylation over these sequences. In order to avoid this, I create a bed file of unique 2kb sequences coordinates that does not contain info of transcripts' names and will be used exclusively to count methylation over promoters in script src/R/promoter_counts.R
cut -f1-3 data/annotations/TSS_900bp.bed | sort | uniq -u > data/annotations/TSS_900bp_uniq.bed

## Rscript which computes the normalized CpG fraction as described in Saxonov et al., PNAS (2006) and creates a data.frame containing this info for the 50 bp overlapping intervals around TSSs, plus informative plots 
#Rscript --vanilla src/R/TSSs_CpG.R 

## Creating bed files of CGIs which do not overlap by more than 50% of their sequence to the 900bp promoter segments
intersectBed -a ${CGI_annotation} -b data/annotations/TSS_900bp.bed -v -f 0.5 | grep -v _GL - | grep -v random - > data/annotations/UCSC_CGI_nonPromoters.bed 
