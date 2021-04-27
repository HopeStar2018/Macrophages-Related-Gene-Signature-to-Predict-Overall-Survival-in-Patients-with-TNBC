rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr5_ClusterProfile_GSE103091')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr4_GSE103091_WGCNA/geneInfo.RData")

color_name <- c('turquoise')

library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(GOplot)
library(ggplot2)
library(org.Hs.eg.db)
library(patchwork)

color_sel <- color_name[1]

genes <- geneInfo$geneSymbol[geneInfo$moduleColor == color_sel] %>% as.character()

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

# patchwork <- p_bp_dot + p_MF_dot + plot_layout(ncol = 1)

Plot_All <- patchwork + plot_annotation(tag_levels = 'A')

ggsave(filename = paste0('Plot_',color_sel,'.pdf'),Plot_All,width = 12,height = 14)

#### Rectome ####
rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr4_ClusterProfile_GSE103091')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr3_GSE103091_WGCNA/geneInfo.RData")

color_name <- c('magenta','blue')

color_sel <- color_name[1]

library(tidyverse)
library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(GOplot)
library(ggplot2)
library(org.Hs.eg.db)
library(ReactomePA)

Hub_genes <- geneInfo$geneSymbol[geneInfo$moduleColor == color_sel] %>% as.character()

eg = bitr(Hub_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

de <- eg$ENTREZID

x <- enrichPathway(gene = de,pvalueCutoff = 0.05,readable = T,organism = "human")

p <- dotplot(x,showCategory = 10)

p

pdf(file=paste0('Result5_',color_name[i],'_dotplot_BP_hubgenes_reactome.pdf'),width = 12,height = 8)

p + scale_x_discrete(labels=function(x) str_wrap(x, width=40))

dev.off()








