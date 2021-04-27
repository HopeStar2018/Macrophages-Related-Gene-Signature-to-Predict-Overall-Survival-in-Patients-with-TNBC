rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step9_GDSC')

library(pRRophetic)
library(ggplot2)
library(cowplot)
library(tidyverse)

# install.packages("ggsignif")

library(ggsignif)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

load("G:/TCGA_database/TCGA_BRCA/BRCA_exp/TPM/TCGA_BRCA_tpm_ensg_mRNA.Rdata")

#### Step 1 ####

TCGA_BRCA_tpm_ensg_mRNA[1:6,1:6]

group=sapply(strsplit(colnames(TCGA_BRCA_tpm_ensg_mRNA),"[.]"),"[",4)

BRCA_tpm_normal <- TCGA_BRCA_tpm_ensg_mRNA[,group=='11A']
BRCA_tmp_cancer <- TCGA_BRCA_tpm_ensg_mRNA[,group=='01A']

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group_5gene.RData")

tmp_group_ClusterA <- tmp_df[tmp_df$Group == 1, ]
tmp_group_ClusterA <- tmp_group_ClusterA$sample %>% as.character()

tmp_group_ClusterB <- tmp_df[tmp_df$Group == 2, ]
tmp_group_ClusterB <- tmp_group_ClusterB$sample %>% as.character()

colnames(BRCA_tmp_cancer) <- substring(colnames(BRCA_tmp_cancer),1,12)

colnames(BRCA_tmp_cancer) <- gsub('[.]','-',colnames(BRCA_tmp_cancer))

BRCA_tmp_cancer_C1 <- subset(BRCA_tmp_cancer,select = tmp_group_ClusterA)
BRCA_tmp_cancer_C2 <- subset(BRCA_tmp_cancer,select = tmp_group_ClusterB)

colnames(BRCA_tpm_normal) <- substring(colnames(BRCA_tpm_normal),1,16)

colnames(BRCA_tpm_normal) <- gsub('[.]','-',colnames(BRCA_tpm_normal))

dat_C1_C2 <- cbind(BRCA_tmp_cancer_C1,BRCA_tmp_cancer_C2)
dat_C1_Nor <- cbind(BRCA_tmp_cancer_C1,BRCA_tpm_normal)
dat_C2_Nor <- cbind(BRCA_tmp_cancer_C2,BRCA_tpm_normal)

dat <- cbind(dat_C1_C2,BRCA_tpm_normal)

#### Step 1 anno data ####

load("H:/PHD_database/Commonly_used_statistical_analysis/AAA_used_RData_Function_code/Annotation_GPL/AnnoTCGA.RData")

colnames(AnnoTCGA)

dat <- dat_C1_C2

dat$ENSG <- rownames(dat)

dat <- merge(AnnoTCGA,dat)

dat[1:6,1:6]

dat <- dat[,-c(1,3)]

dat <- aggregate(.~GeneName,dat,median)

rownames(dat) <- dat$GeneName

dat <- dat[,-1]

colnames(dat) <- gsub('[.]','-',colnames(dat))

dat_C1_C2 <- dat

save(dat_C1_C2,file = 'Result9_dat_5gene_C1_C2.RData')

########################

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step9_GDSC/Result9_dat_5gene_All.RData")

ann <- data.frame(ImmClust = c(rep('TNBC_ClusterA',ncol(BRCA_tmp_cancer_C1)),
                               rep('TNBC_ClusterB',ncol(BRCA_tmp_cancer_C2)),
                                rep('Normal',ncol(BRCA_tpm_normal))),
                   sample = c(colnames(BRCA_tmp_cancer_C1),
                              colnames(BRCA_tmp_cancer_C2),
                              colnames(BRCA_tpm_normal)))

rownames(ann) <- ann$sample

#药物名字
GCP.drug <- read.table("drug.txt") #如果要例文的两种药物，就换成drug_eg.txt
GCP.drug <- GCP.drug$V1

# 自定义足够多的box的颜色，颜色数量至少等于分组数量
jco <- c("#EABF00", "#2874C5", "red")

### 药物敏感性预测 ###、
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表
plotp <- list()
Result_Ttest <- data.frame()

# dat <- dat_C2_Nor

for (drug in GCP.drug) {
  
  drug <- GCP.drug[64]
  
  set.seed(20200715) # 因为预测过程默认10-fold CV，所以设置种子以便结果可重复
  
  cat(drug," starts!\n") # 提示当前药物已开始分析
  
  # 预测IC50值，采用默认参数，详细可参考??pRRopheticPredict参数列表
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = as.matrix(dat[,rownames(ann)]),
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              selection = 1) # 1表示若有重复基因取均值处理
  
  if(!all(names(predictedPtype[[drug]])==rownames(ann))) {stop("Name mismatched!\n")} # 若名字不匹配则报错退出
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        
                                        "ImmClust"=ann$ImmClust, # 这里我修改了C1和C2的名字
                                        
                                        row.names = names(predictedPtype[[drug]])) 
  
  predictedBoxdat[[drug]]$ImmClust <- factor(predictedBoxdat[[drug]]$ImmClust,levels = unique(ann$ImmClust),ordered = T) # 把类改成因子变量
  
  pvalue <- t.test(predictedBoxdat[[drug]]$est.ic50 ~ predictedBoxdat[[drug]]$ImmClust)
  
  tmp <- data.frame(Drug_name = drug,
                    
                    pvalue = pvalue$p.value)
  
  Result_Ttest <- rbind(Result_Ttest,tmp)
  
  rm(tmp)
  
  p_value <- pvalue$p.value
  
  if(p_value >= 0.05){
    
    next
    
  }else{
    
    # 绘图
    p <- ggplot(data = predictedBoxdat[[drug]], aes(x=ImmClust, y=est.ic50))
    
    p <- p + geom_boxplot(aes(fill = ImmClust)) + 
      
      scale_fill_manual(values = jco[1:length(unique(ann$ImmClust))]) + #自定义box的配色
      
      theme(legend.position="none") + # 倾斜字体
      
      theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
      
      xlab("") + ylab("Estimated IC50") + 
      
      ggtitle(drug)+ 
    
      geom_signif(comparisons = list(c("TNBC_ClusterA","TNBC_ClusterB"),
                                     c("TNBC_ClusterA","Normal"),
                                     c("Normal","TNBC_ClusterB")),test = "wilcox.test")
    
    ###
    pdf(file = 'Vinorelbine_C2_Nor.pdf',height = 12,width = 6,onefile = F)
    p
    dev.off()
    
    ###
    
    
    plotp[[drug]] <- p # 保存在列表里供合并图片用
    
    cat(drug," has been finished!\n") # 提示当前药物已分析结束
    
  }
  
  
}

save(plotp,Result_Ttest,file = 'Result9_GDSC_5gene_C1_C2.RData')

length(plotp)
names(plotp)

# 合并图片
#适合展示两种药物
# p1 <- plot_grid(plotp[[1]],plotp[[2]],labels = c("A","B"),nrow = 1) # title可以AI下拉到合适位置，就如例文所示
# ggsave("boxplot of predicted IC50.pdf", width = 6, height = 5)

# 适合展示多种药物
# p2 <- plot_grid(plotlist=plotp, ncol=5)
# ggsave("boxplot of predicted IC50_multiple.pdf", width = 12, height = 6)


#  11 阿西替尼；17 贝沙罗汀；20 比卡鲁胺；23 博来霉素；26 苔藓抑素；36 阿糖胞苷；38 多西他赛；39 阿霉素；43 埃罗替尼；
# 45 吉西他滨；49 伊马替尼；55 拉帕替尼；56 来那度胺；57 二甲双胍；58甲氨蝶呤；61丝裂霉素；63 尼罗替尼；80 选择性elF2α去磷酸化抑制剂；
#  83索拉非尼；85 舒尼替尼；86西罗莫司脂化物；87 替吡法尼；89长春瑞滨

rm(list = ls())

setwd("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step9_GDSC")

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step9_GDSC/Result9_GDSC.RData")

Result_Ttest <- Result_Ttest[order(Result_Ttest$pvalue), ]

plotlist_mean <- plotp[c(1:19)]

plotlist_mean[[1]]

names(plotlist_mean)

#####

plotlist_mean <- plotp[Result_Ttest$Drug_name[1:19]]

library(cowplot)

p2 <- plot_grid(plotlist=plotlist_mean, ncol=5)

p2

ggsave("C1 C2 boxplot of predicted IC50_multiple.pdf", width = 12, height = 36)





