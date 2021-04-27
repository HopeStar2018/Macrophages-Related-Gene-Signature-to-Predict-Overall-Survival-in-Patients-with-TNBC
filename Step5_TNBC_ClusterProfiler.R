rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step5_TNBC_ClusterProfiler')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step4_TNBC_WGCNA/TCGA_BRCA_CIBERSORT.ABS/geneInfo.RData")

color_name <- 'blue'

library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(GOplot)
library(ggplot2)
library(org.Hs.eg.db)
library(patchwork)

genes <- geneInfo$geneSymbol[geneInfo$moduleColor == color_name] %>% as.character()

entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)

gene <- as.character(entrezIDs)

## BP ##
kk <- enrichGO(gene = gene,OrgDb = org.Hs.eg.db,ont = "BP", pvalueCutoff =0.05, qvalueCutoff = 0.05)

p_bp_dot <- dotplot(kk,showCategory = 10)

p_bp_dot

### CC
kk <- enrichGO(gene = gene,OrgDb = org.Hs.eg.db,ont = "CC", pvalueCutoff =0.05, qvalueCutoff = 0.05)

p_CC_dot <- dotplot(kk,showCategory = 10)

p_CC_dot

### MF
kk <- enrichGO(gene = gene,OrgDb = org.Hs.eg.db,ont = "MF", pvalueCutoff =0.05, qvalueCutoff = 0.05)

p_MF_dot <- dotplot(kk,showCategory = 10)

p_MF_dot

patchwork <- p_bp_dot + p_CC_dot + p_MF_dot + plot_layout(ncol = 1)

Plot_All <- patchwork + plot_annotation(tag_levels = 'A')

ggsave(filename = 'Plot_Blue.pdf',Plot_All,width = 12,height = 14)






















