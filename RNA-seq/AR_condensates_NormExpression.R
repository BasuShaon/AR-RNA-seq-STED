################################################################################
## Normalized expression code to produce Figure 6F and extended Figure 7D
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
BiocManager::install(version = "3.17")
BiocManager::install('msigdbr', force = TRUE)
library(ggplot2)
library(tidyverse)
library(msigdbr)
library(readxl)
library(reshape2)
library(dplyr)

# read in normalized readcounts
# read in count data from GSE206853
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206853
start = 
  read_excel("GEO_submission/GSE206853_SupplementTable3_ARcondensates_220531a.xlsx",
             sheet = 2)

df <- start[,3:ncol(start)]

# norm functions
z_score_scale_row <- function(row) {
  (row - mean(row)) / sd(row)
}

mean_normalize_row <- function(row) {
  row - mean(row, na.rm = TRUE)
}

scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

# format data
scaled_df <- as.data.frame(t(apply(df, 1, z_score_scale_row)))
scaled_df$GeneID <- start$GeneID
scaled_df <- na.omit(scaled_df)
input <- scaled_df
input_6 <-input[,c(grep('6',names(input),value=TRUE),'GeneID')] 
input_24 <- input[,c(grep('24',names(input),value=TRUE),'GeneID')] 
input_melt <- melt(input_24, variable.name = "Column", value.name = "Value")
input_melt$Group <- sub("_[^_]*$", "", input_melt$Column)
input_agg <- aggregate(Value ~ GeneID + Group, data = input_melt, FUN = mean)
input_agg$Drug <- sub("_.+$", "", gsub("^[^_]+_", "", input_agg$Group))
input_agg$Concentration <- sub(".*_", "", input_agg$Group)
input_agg$Concentration[grepl("DMSO", input_agg$Drug)] <- '0'
input_agg$Drug<- as.factor(input_agg$Drug)
input_agg$GeneID <- as.factor(input_agg$GeneID)
input_agg$Concentration <- as.numeric(input_agg$Concentration)
input_DMSO <- subset(input_agg, Drug == 'DMSO')
input_DMSO_1ae <- input_DMSO
input_DMSO_EPI <- input_DMSO
input_DMSO_1ae$Drug <- '1ae'
input_DMSO_EPI$Drug <- 'EPI'
input_remain <- subset(input_agg, Drug != 'DMSO')
input_clean <- rbind(input_remain, input_DMSO_1ae, input_DMSO_EPI)
input_clean$Drug <- as.factor(input_clean$Drug)


# Get the hallmark gene sets for humans
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

# Extract genes lists from msigdb
AR_genes <- subset(hallmark_sets, gs_name == "HALLMARK_ANDROGEN_RESPONSE")$ensembl_gene
Inflamm_genes <- subset(hallmark_sets, gs_name == "HALLMARK_INFLAMMATORY_RESPONSE")$ensembl_gene
Sperm <- subset(hallmark_sets, gs_name == "HALLMARK_SPERMATOGENESIS")$ensembl_gene
p53 <- subset(hallmark_sets, gs_name == "HALLMARK_P53_PATHWAY")$ensembl_gene
E2F_genes <- subset(hallmark_sets, gs_name == "HALLMARK_E2F_TARGETS")$ensembl_gene
MYC1_genes <- subset(hallmark_sets, gs_name == "HALLMARK_MYC_TARGETS_V1")$ensembl_gene
MYC2_genes <- subset(hallmark_sets, gs_name == "HALLMARK_MYC_TARGETS_V2")$ensembl_gene
G2M_genes <- subset(hallmark_sets, gs_name == "HALLMARK_G2M_CHECKPOINT")$ensembl_gene

# Set List var to pathway of interest
List <- AR_genes #HALLMARK_ANDROGEN_RESPONSE
Target <- input_clean[input_clean$GeneID %in% List, ]
Target$Drug <- factor(Target$Drug, levels = c('1ae','EPI'))

# Line plot for pathway
ggplot(Target, aes(x = Concentration, y = Value,  color = Drug, group = interaction(GeneID, Drug))) +
  geom_line(alpha = 0.5) + 
  labs(title = "Gene Expression by Drug and Dose", x = "Dose", y = "Expression Value") +
  theme_minimal() +
  theme(legend.position = 'none')
  

