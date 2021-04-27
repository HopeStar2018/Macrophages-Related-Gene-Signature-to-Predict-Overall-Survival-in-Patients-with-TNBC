rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step12_TMB')

#source("https://bioconductor.org/biocLite.R")
#biocLite("TCGAbiolinks")
#biocLite("maftools")
require(TCGAbiolinks)
require(maftools)

#下载突变数据
BRCA_mutect2 <- GDCquery_Maf(tumor = "BRCA", pipelines = "mutect2")

#确认是否都是somatic mutation
levels(factor(BRCA_mutect2$Mutation_Status))

save(BRCA_mutect2,file = 'BRCA_mutect2.RData')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step12_TMB/BRCA_mutect2.RData")

require(dplyr)

#用maftools里的函数读取maf文件
var_maf <- read.maf(maf = BRCA_mutect2, isTCGA = T)
str(var_maf)

#并且已经统计好了每个sample里的variants数量
variants_per_sample <- var_maf@variants.per.sample
dim(variants_per_sample)
#保存到文件
write.csv(variants_per_sample, "easy_input_mut1.csv", quote = F, row.names = F)

library(stringr)

#突变

myMut <- read.csv("easy_input_mut1.csv")
colnames(myMut) <- c("Tumor_Sample_Barcode","Variants")
#保留barcode的前三个label
myMut$Tumor_Sample_Barcode <- str_sub(myMut$Tumor_Sample_Barcode,1, 12)
head(myMut)

#分组
load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group_5gene.RData")

TMB_per_sample <- myMut
TMB_per_sample$TMB <- myMut$Variants%/%38

#把TMB值保存到文件，自己设定阈值，就可以用高低TMB分组进行生存分析
write.csv(TMB_per_sample, "TMB_output.csv", quote = F, row.names = F)

colnames(TMB_per_sample)
colnames(tmp_df)[2] <- 'Tumor_Sample_Barcode'

TMB_clinical_mRNA <- merge(TMB_per_sample, tmp_df, by="Tumor_Sample_Barcode")
head(TMB_clinical_mRNA)

#可以先用‘Mann-Whitney’ test看一下pvalue
res <- wilcox.test(TMB ~ Group, data = TMB_clinical_mRNA,
                   paired = FALSE, #whether you want a paired test
                   exact = FALSE)
(pvalue <- res$p.value)

require(ggplot2)
require(ggpubr)
require(ggsignif)

jco <- c("#EABF00", "#2874C5", "red")

TMB_clinical_mRNA$Group <- ifelse(TMB_clinical_mRNA$Group == 1,'TNBC_ClusterA','TNBC_ClusterB')

p <- ggplot(data = TMB_clinical_mRNA, aes(x=Group, y=TMB))

p <- p + geom_boxplot(aes(fill = Group)) + 
  
  scale_fill_manual(values = jco[1:length(unique(TMB_clinical_mRNA$Group))]) + #自定义box的配色
  
  theme(legend.position="none") + # 倾斜字体
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
  
  xlab("") + ylab("TMB") + 
  
  geom_signif(comparisons = list(unique(TMB_clinical_mRNA$Group)),test = "t.test")

p





















