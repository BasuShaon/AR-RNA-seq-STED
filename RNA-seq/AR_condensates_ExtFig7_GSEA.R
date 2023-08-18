################################################################################
## GSEA code to produce Extended Figure 7C and 7E
## Basu et al. 2023
## Rational optimization of a transcription factor optimization domain inhibitor
## Author: Shaon Basu
################################################################################

# dependencies and environment
install.packages('readxl')
install.packages('gglplot2')
install.packages('tidyverse')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler", force = TRUE)
library(ggplot2)
library(readxl)
library(tidyverse)
library(msigdbr)
library(org.Hs.eg.db)
library(DOSE)
library(clusterProfiler)
library(enrichplot)

#set organism
organism = 'org.Hs.eg.db'

# read in differential expression lists from Supplemental Table S3
setwd("C:/Users/Shaon/Desktop/AR NSMB/GEO_submission/")
df1 <- read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 1)
df2 <- read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 2)

# extract log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$Gene

# omit any NA values, sort
gene_list<-na.omit(original_gene_list)
gene_list_sorted <- sort(gene_list, decreasing = TRUE)

# ensembl anontation from msigdbr, hallmark pathways
H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

# GSEA with 10,000 permutations, p-val cut off 0.05
result <- GSEA(gene_list_sorted, TERM2GENE = H, 
               nPerm = 10000, pvalueCutoff = 0.05)

# Plot GSEA curve (Ext Fig 7C)
gseaplot2(result, geneSetID = "HALLMARK_ANDROGEN_RESPONSE", subplots = 1:2)

# Plot Dotplots of Activated / Supressed pathways (Ext Fig 7E)
dotplot(result, showCategory = 10, font.size = 8,
        x = "GeneRatio",   # option -> c("GeneRatio", "Count")
        color = "p.adjust") 


