################################################################################
## Diff exp volcano plot code to produce Figure 6D
## Basu et al. 2023
## Rational optimization of a transcription factor optimization domain inhibitor
## Author: Shaon Basu
################################################################################

# dependencies and environment
install.packages('readxl')
install.packages('gglplot2')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")
BiocManager::install("EnhancedVolcano")
library(ggplot2)
library(readxl)
library(DESeq2)
library(EnhancedVolcano)

# read in differential expression lists from Supplemental Table S3, peform GSEA

DE_24h_1ae_5uM_vsDMSO <- read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 1)
DE_24h_EPI_25uM_vsDMSO  <- read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 2)
DE_24h_1ae_0p5uM_vsDMSO  <- read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 3)
DE_24h_EPI_2p5uM_vsDMSO  <- read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 4)
DE_6h_1ae_5uM_vsDMSO <- read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 5)
DE_6h_EPI_25uM_vsDMSO <- read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 6)
DE_6h_1ae_0p5uM_vsDMSO <- read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 7)
DE_6h_EPI_2p5uM_vsDMSO <- read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 8)

# generate volcano plots for comparisons shown in Figure 7D 
EnhancedVolcano(DE_24h_1ae_5uM_vsDMSO,
                lab = DE_24h_1ae_5uM_vsDMSO$Gene,
                x = 'log2FoldChange',
                ylim = c(0, 200),
                xlim = c(-5, 5),
                y = 'padj')

EnhancedVolcano(DE_24h_EPI_25uM_vsDMSO,
                lab = NA,
                x = 'log2FoldChange',
                ylim = c(0, 200),
                xlim = c(-5, 5),
                y = 'padj')
