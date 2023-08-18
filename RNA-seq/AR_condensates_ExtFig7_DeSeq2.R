################################################################################
## DESeq2 loading, transformation, and pca to produce Extended Data Figure 7B
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
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("apeglm")
library(ggplot2)
library(readxl)
library(DESeq2)
library(EnhancedVolcano)
library(apeglm)

# read in count data from GSE206853
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206853
input = 
  read_excel("GSE206853_SupplementTable3_ARcondensates_220531a.xlsx",
                   sheet = 1)

# create count matrix
countmatrix = as.matrix(input)
countmatrix <- countmatrix[, -c(1,2)]
original_colnames <- colnames(countmatrix)
countmatrix <- matrix(as.numeric(countmatrix), ncol=ncol(countmatrix))
colnames(countmatrix) <- original_colnames
rownames(countmatrix) <- input$GeneID

# Read in design matrix (deposited to Git as 'GSE206853_designmatrix.xlsx') 
# Create sample annotation
coldata = as.data.frame(read_excel("GSE206853_designmatrix.xlsx"))
coldata$Sample <- factor(coldata$Sample)
coldata$Condition <- factor(coldata$Condition)
rownames(coldata) <- coldata$Sample
colnames(countmatrix) <- rownames(coldata)

# Load and execute DESeq statistical model
dds <- DESeqDataSetFromMatrix(countData = countmatrix,
                              colData = coldata,
                              design = ~ Condition) 
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

# Perform principal component analysis
vsd <- vst(dds) # variance stabilizing transformation 
plotPCA(vsd, intgroup = "Condition") + theme_bw() # PCA
get_pca(vsd)
