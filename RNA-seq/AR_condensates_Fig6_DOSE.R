################################################################################
## GSEA code to produce Figure 6E
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
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("apeglm")
BiocManager::install('fgsea')
BiocManager::install('org.Hs.eg.db')
BiocManager::install('msigdbr', force = TRUE)
library(ggplot2)
library(readxl)
library(DESeq2)
library(fgsea)
library(tidyverse)
library(msigdbr)
library(org.Hs.eg.db)

# GSEA command
gsea <- function(x)
{
results <- x
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=results$GeneID, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")

res <- left_join(results, ens2symbol, by=c("GeneID"="ENSEMBL"))

res2 <- res %>% 
  dplyr::select(SYMBOL, log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(log2FoldChange))
res2
ranks <- deframe(res2)

msigdbr_df <- msigdbr(species = "human", category = "H")
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

fgseaRes <- fgsea(pathwaysH, ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
}
 
# read in differential expression lists from Supplemental Table S3, peform GSEA

DE_24h_1ae_5uM_vsDMSO <- gsea(read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 1))
DE_24h_EPI_25uM_vsDMSO  <- gsea(read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 2))
DE_24h_1ae_0p5uM_vsDMSO  <- gsea(read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 3))
DE_24h_EPI_2p5uM_vsDMSO  <- gsea(read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 4))
DE_6h_1ae_5uM_vsDMSO <- gsea(read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 5))
DE_6h_EPI_25uM_vsDMSO <- gsea(read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 6))
DE_6h_1ae_0p5uM_vsDMSO <- gsea(read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 7))
DE_6h_EPI_2p5uM_vsDMSO <- gsea(read_excel("SupplementTable3_ARcondensates_220531a.xlsx", sheet = 8))


# merge adjusted p-values and normalized enrichment scores 
rbind_with_name <- function(df_name) {
  df <- get(df_name)
  df$origin <- df_name
  return(df)
}
df_names <- c('DE_24h_1ae_5uM_vsDMSO','DE_24h_EPI_25uM_vsDMSO','DE_24h_1ae_0p5uM_vsDMSO','DE_24h_EPI_2p5uM_vsDMSO',
              'DE_6h_1ae_5uM_vsDMSO','DE_6h_EPI_25uM_vsDMSO','DE_6h_1ae_0p5uM_vsDMSO','DE_6h_EPI_2p5uM_vsDMSO')
df <- do.call(rbind, lapply(df_names, rbind_with_name))
df$pathway <- as.factor(df$pathway)
df$origin <- as.factor(df$origin)


selectedPathways <- c('HALLMARK_E2F_TARGETS','HALLMARK_MYC_TARGETS_V1','HALLMARK_ANDROGEN_RESPONSE','HALLMARK_G2M_CHECKPOINT',
                      'HALLMARK_MYC_TARGETS_V2','HALLMARK_MITOTIC_SPINDLE','HALLMARK_SPERMATOGENESIS','HALLMARK_NOTCH_SIGNALING',
                      'HALLMARK_COMPLEMENT','HALLMARK_UV_RESPONSE_UP','HALLMARK_P53_PATHWAY','HALLMARK_APOPTOSIS','HALLMARK_COAGULATION',
                      'HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_HEME_METABOLISM','HALLMARK_XENOBIOTIC_METABOLISM',
                      'HALLMARK_HEDGEHOG_SIGNALING','HALLMARK_ALLOGRAFT_REJECTION','HALLMARK_INFLAMMATORY_RESPONSE')

subset_df <- subset(df, pathway %in% selectedPathways)
subset_df$pathway <- factor(subset_df$pathway, levels = rev(selectedPathways))
subset_df$origin <- factor(subset_df$origin, levels = df_names)

#dotplot to show enriched and depleted pathways by GSEA stat
ggplot(subset_df, aes(x = pathway, y = origin, size = -log10(padj), color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = 'blue', mid = "white", high = "red", midpoint = 0) +
  theme_minimal() + 
  labs(size = "-log10(p-value)", color = "NES") +
  theme(legend.position = "right") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


