rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step8_Subset_GSEA')

library(tidyverse)

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step4_TNBC_WGCNA/Result4_rawdata.RData")

rm(Immunity_List)

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group_5gene.RData")

colnames(TCGA_BRCA_tpm_symbol_mRNA) <- substring(colnames(TCGA_BRCA_tpm_symbol_mRNA),1,12)

df_uniq <- subset(TCGA_BRCA_tpm_symbol_mRNA,select = rownames(tmp_df))

rm(TCGA_BRCA_tpm_symbol_mRNA)

NAME <- rownames(df_uniq) %>% as.character()

Description <- rep("NA",nrow(df_uniq))

input <- cbind(NAME=NAME,Description=Description,df_uniq)

# colnames(input) <- gsub("[.]","-",colnames(input))
colnames(input)[1:5]

dim(df_uniq)

write.table(input,file = "input.gct",row.names = F,sep = "\t",quote = F)

##

tmp_df$Group <- ifelse(tmp_df$Group == 1,'ClusterA','ClusterB')

Group_singlegene <- tmp_df$Group
tmp <- unique(Group_singlegene)

## output ##
sink(paste0('input.cls'),append = T)

line1 <- c(ncol(df_uniq),2,1)

cat(paste(line1,sep = '\t'),'\n')

line2 <- c('#',tmp)

cat(paste(line2,sep = '\t'),'\n')

line3 <- Group_singlegene

cat(paste(line3,sep = '\t'))

sink()


















