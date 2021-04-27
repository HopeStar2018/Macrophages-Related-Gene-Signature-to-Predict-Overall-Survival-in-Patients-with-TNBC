rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr6_VennPlot')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr4_GSE103091_WGCNA/geneInfo.RData")

color_name <- c('turquoise')

color_sel <- color_name[1]

GSE103091_genes <- geneInfo$geneSymbol[geneInfo$moduleColor == color_sel] %>% as.character()

##

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step4_TNBC_WGCNA/TCGA_BRCA_CIBERSORT/geneInfo.RData")

TCGA_genes <- geneInfo$geneSymbol[geneInfo$moduleColor == 'blue'] %>% as.character()

library(VennDiagram)

pdf("TCGA_GSE103091_venn.pdf", width = 8, height = 8)
grid.newpage()
venn.plot <- venn.diagram(list(GSE103091 = GSE103091_genes, #画图
                               TCGA = TCGA_genes), NULL, 
                          fill = c("#E31A1C","#E7B800"), 
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names = c("GSE103091", "TCGA"), 
                          main = "Overlap")
grid.draw(venn.plot)

dev.off()

tmp_gene <- intersect(TCGA_genes,GSE103091_genes)

save(tmp_gene,file = 'intersect_genes.RData')






















