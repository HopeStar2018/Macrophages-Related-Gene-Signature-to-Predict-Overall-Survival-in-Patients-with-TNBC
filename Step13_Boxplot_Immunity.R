rm(list = ls())

setwd("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step13_Boxplot_immunity")

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step3_TNBC_Timer/TCGA_BRCA_Immunity.RData")

rownames(TCGA_BRCA_CIBERSORT.ABS)

group=sapply(strsplit(rownames(TCGA_BRCA_CIBERSORT.ABS),"\\-"),"[",4)

TCGA_BRCA_CIBERSORT.ABS <- TCGA_BRCA_CIBERSORT.ABS[group == '01', ]

rownames(TCGA_BRCA_CIBERSORT.ABS) <- substring(rownames(TCGA_BRCA_CIBERSORT.ABS),1,12)

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group.RData")

TCGA_BRCA_CIBERSORT.ABS$sample <- rownames(TCGA_BRCA_CIBERSORT.ABS)

df <- merge(tmp_df,TCGA_BRCA_CIBERSORT.ABS)

rownames(df) <- df$sample

df <- df[,-1]

df$Group <- ifelse(df$Group == 1,'TNBC_ClusterA','TNBC_ClusterB')

library(tidyr)
library(reshape2)

colnames(df) <- gsub('_CIBERSORT.ABS','',colnames(df))

colnames(df)

df_melt <- melt(df)

jco <- c("#EABF00", "#2874C5", "red")

library(ggplot2)
library(ggpubr)

p <- ggplot(data = df_melt, aes(x=variable, y=value))

p2 <- p + geom_boxplot(aes(fill = Group)) + 
  
  scale_fill_manual(values = jco[1:2]) + #自定义box的配色
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        plot.title = element_text(size = 12, hjust = 0.5),
        panel.grid = element_blank()) +
  
  xlab("") + ylab("")

p3 <- p2 + stat_compare_means(aes(group = Group ,label=..p.signif..))

ggsave('Immunity_boxplot.pdf',p3,width = 24,height = 12)













