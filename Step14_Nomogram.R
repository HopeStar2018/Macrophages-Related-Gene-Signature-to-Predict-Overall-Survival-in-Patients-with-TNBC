rm(list = ls())

setwd("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step14_Nomogram")

library(tidyverse)
library(pheatmap)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

Clincal_BRCA <- read.csv('/Volumes/HS_SSD/Other_analysis/XL/Analysis_20200512/20200529_Analysis/Step1_Rawdata/BRCA_Clinical.csv',row.names = 1,na.strings = '')
Clincal_BRCA[1:5,1:5]

sampleInfo <- Clincal_BRCA[Clincal_BRCA$sample_type == 'Primary Tumor',]

rm(Clincal_BRCA)

sampleInfo <- na.omit(sampleInfo)

sampleInfo <- sampleInfo[sampleInfo$futime_os > 30, ]

sampleInfo <- sampleInfo[sampleInfo$Gender != 'MALE', ]

rownames(sampleInfo) <- substring(rownames(sampleInfo),1,12)

table(sampleInfo$Stage_M)

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/Result6_GetFactors.Rdata")

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step4_TNBC_WGCNA/Result4_tidyrawdata_exp.RData")

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group.RData")

sampleInfo <- sampleInfo[rownames(tmp_df), ]

## Table Baseline output ##

source(file = '/Volumes/HS_SSD/Other_analysis/JLX_Nomogram/Lesson1/Part1/Lesson1_Function/GetTable1.R')

sampleInfo$AJCC_Status[which(sampleInfo$AJCC_Status == "Stage I" | sampleInfo$AJCC_Status == "Stage IA" | sampleInfo$AJCC_Status == "Stage IB")] <- 'StageI'
sampleInfo$AJCC_Status[which(sampleInfo$AJCC_Status == "Stage IIA" | sampleInfo$AJCC_Status == "Stage IIB"| sampleInfo$AJCC_Status == "Stage II")] <- 'StageII'
sampleInfo$AJCC_Status[which(sampleInfo$AJCC_Status == "Stage IIIA" | sampleInfo$AJCC_Status == "Stage IIIB" | sampleInfo$AJCC_Status == "Stage IIIC" | sampleInfo$AJCC_Status == "Stage III")] <- 'StageIII'
sampleInfo$AJCC_Status[which(sampleInfo$AJCC_Status == "Stage IV")] <- 'StageIV'
sampleInfo$AJCC_Status[which(sampleInfo$AJCC_Status == "Stage X")] <- ''
table(sampleInfo$AJCC_Status)

GetTable1(sampleInfo,'Table_TCGA_TNBC_Baseline')

datExpr <- subset(datExpr,select = GetFactors)

# save.image('Chr1_tidyverse.RData')

##########

#### Step 1 Lasso Cluster ####

datExpr$sample <- rownames(datExpr)
sampleInfo$sample <- rownames(sampleInfo)

df <- merge(tmp_df,datExpr)
rownames(df) <- df$sample
df <- df[,-1]

df$Group <- ifelse(df$Group == 1, 'TNBC_ClusterA','TNBC_ClusterB')

i <- 1

cv.times <- 1000

cv.list_snowfall <- list()

library(snowfall)
library(parallel)
library(glmnet)

n_cores <- detectCores()/2
sfInit(parallel = TRUE,cpus = n_cores)

sfExport('df')
sfExport('i')
sfExport('cv.times')
sfExport('cv.list_snowfall')

sfLibrary(survival)
sfLibrary(survivalROC)
sfLibrary(penalized)
sfLibrary(org.Hs.eg.db)
sfLibrary(glmnet)

cv.list_snowfall <- sfLapply(1:1000,function(i){
  
  outcome <- df$Group
  
  xx <- as.matrix(df[,2:ncol(df)])
  
  cv.fit1=cv.glmnet(xx,outcome,family="binomial")
  
  coef.cv <- as.matrix(coef(cv.fit1, s = 'lambda.min'))
  
  coef.cv <- coef.cv[coef.cv[,1] != 0, ]
  
  # print(coef.cv)
  
  nameCoefCv <- names(coef.cv)
  
  if(is.null(nameCoefCv)) {nameCoefCv <- 'isNA'}
  
  return(nameCoefCv)
  
})

sfStop()

save.image(file = 'lasso_TNBC_snowfall_result.RData')

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step14_Nomogram/lasso_TNBC_snowfall_result.RData")


# build a df, whose row = lasso model; column = gene occured or not
# sum by row and column

allGenes <- unique(unlist(cv.list_snowfall))

xxx.df <- as.data.frame(matrix(NA, ncol = length(allGenes)))

colnames(xxx.df) <- c(allGenes)

for (i in 1:length(cv.list_snowfall)) {
  if(is.null(cv.list_snowfall[[i]])) {
    xxx.df[i, ] <- F
    xxx.df[i, which(colnames(xxx.df) == 'isNA')] <- T}
  else {
    xxx.df[i, ] <- F
    xxx.df[i, which(colnames(xxx.df) %in% cv.list_snowfall[[i]])] <- T
  }
}

geneSum <- apply(xxx.df[,1:ncol(xxx.df)], 2, sum)

geneSum <- geneSum[order(geneSum, decreasing = T)]

xxx.df$SUM <- apply(xxx.df, 1, sum)


gModel.sum <- c()

for (i in 1: length(cv.list_snowfall)) {
  gModel <- paste0(cv.list_snowfall[[i]], collapse = ' ')
  gModel.sum <- c(gModel.sum, gModel)
}

gModel.sum[nchar(gModel.sum) == 0] <- 'isNA'

sort(table(gModel.sum, useNA = 'ifany'), decreasing = T)

xxx.df$gModel <- gModel.sum

xxx.df <- xxx.df[,c(ncol(xxx.df), 1:(ncol(xxx.df) - 1))]


geneSum <- sort(geneSum, decreasing = T)
View(geneSum)
#geneSum <- geneSum[-1] ## del 'SUM'

xxA.gene <- sort(names(geneSum)[geneSum == 1000])

select_lncRNA <- data.frame(select_feasures = xxA.gene)

# write.table(select_lncRNA,"select_feasures_LASSO.txt")

Lasso_Cluster_feature <- select_lncRNA$select_feasures[-1]

save(Lasso_Cluster_feature,file = 'Lasso_Cluster_feature.RData')


yyy <- names(table(xxx.df$gModel)[order(table(xxx.df$gModel), decreasing = T)])

yyy.df <- as.data.frame(matrix(NA, ncol = length(allGenes), nrow = length(yyy)))

rNames_yyyDf <- c()

gm_count_sum <- c()

colnames(yyy.df) <- allGenes

for (i in 1:length(yyy)) {
  gm <- strsplit(yyy[i], ' ')
  gm <- gm [[1]]
  gm_count <- length(gm)
  gm_count_sum <- c(gm_count_sum, gm_count)
  
  yyy.df[i, ] <- F
  for(j in 1:length(allGenes)) {
    if(colnames(yyy.df)[j] %in% gm) yyy.df[i, j] <- T
  }
  
  rNames_yyyDf[i] <- paste0(gm_count, " feasures")
}

rNames_yyyDf

for (i in 2:length(rNames_yyyDf)) {
  if (rNames_yyyDf[i] %in% rNames_yyyDf[1:(i-1)]) {
    rNames_yyyDf[i] <- paste0(rNames_yyyDf[i], '*')
  }
}

rNames_yyyDf

rownames(yyy.df) <- rNames_yyyDf

yyy.df <- yyy.df[, order(colnames(yyy.df))]

yyy.df <- as.matrix(yyy.df)

y.df <- apply(yyy.df, 2, as.numeric)

rownames(y.df) <- rNames_yyyDf

library(ComplexHeatmap)

hp1 <- Heatmap(y.df, col = c('white', 'yellow'), cluster_rows = F, cluster_columns = F, rect_gp = gpar(col = "grey", lty = 1, lwd = 2), row_names_side = 'left', show_heatmap_legend = F)

gModel_Count <- table(xxx.df$gModel)[order(table(xxx.df$gModel), decreasing = T)]

gModel_Count <- as.matrix(sort(gModel_Count, decreasing = T))

rownames(gModel_Count) <- paste0(gModel_Count[,1], ' times')

hp2 <- Heatmap(gModel_Count, col = c('pink', 'red'), cluster_rows = F, cluster_columns = F, row_names_side = 'right', show_column_names = F, show_heatmap_legend = F)

hp_list <- hp1 + hp2

pdf('OP_Pt01_01_lassoModel.pdf',width = 12,height = 12, onefile=FALSE)

draw(hp_list)

dev.off()


pdf('OP_Pt01_02_lassoGenes.pdf',width = 12,height = 12, onefile=FALSE)

barplot(height = geneSum, col = 'cyan', las = 2,angle = 60,ylim = c(0,1000)) 

# abline(h = 900, lty = 2)

dev.off()

#### Step 2 PCoA ####

rm(list = ls())

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step14_Nomogram/Lasso_Cluster_feature.RData")

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step4_TNBC_WGCNA/Result4_tidyrawdata_exp.RData")

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group.RData")

#定义足够多的颜色，用于展示分组
mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

# install.packages('ade4')
# install.packages('vegan')

library(ade4)
library(ggplot2)
library(RColorBrewer)
library(vegan)#用于计算距离

tab <- subset(datExpr,select = Lasso_Cluster_feature)

tmp_df$Group <- ifelse(tmp_df$Group == 1,'TNBC_ClusterA','TNBC_ClusterB')

#### PLS-DA ####

library(mixOmics)

plsda_result <-plsda(tab, tmp_df$Group, ncomp = 2)

plsda_result

plotIndiv(plsda_result, ind.names = TRUE, style = 'ggplot2')

# plotIndiv(plsda_result, ind.names = TRUE, style = '3d')

plsda_result_eig <- {plsda_result$explained_variance$X}[1:2]

sample_site <- data.frame(plsda_result$variates)[1:2]

#为样本点坐标添加分组信息
sample_site$Sample_ID <- rownames(sample_site)

names(sample_site)[1:2] <- c('plsda1', 'plsda2')

colnames(tmp_df)[2] <- 'Sample_ID'

sample_site <- merge(sample_site, tmp_df, by = 'Sample_ID', all.x = TRUE)

library(ggplot2)

#使用 ggplot2 简单绘制 PLS-DA 结果图
plsda_plot <- ggplot(sample_site, aes(plsda1, plsda2, color = Group, label = Sample_ID)) +
  geom_point(size = 1.5, alpha = 0.6) + 
  stat_ellipse(show.legend = FALSE) +    #添加 95% 置信椭圆
  scale_color_manual(values = c('#1D7ACC', '#F67433', '#00815F')) +
  theme(panel.grid = element_line(color = 'grey50'), panel.background = element_rect(color = 'black', fill = 'transparent')) + 
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) +
  labs(x = paste('PLS-DA axis1 ( explained variance ', round(100 * plsda_result_eig[1], 2), '% )', sep = ''), y = paste('PLS-DA axis2 ( explained variance ', round(100 * plsda_result_eig[2], 2), '% )', sep = ''))

ggsave('Result2_plsda_plot.pdf', plsda_plot, width = 12, height = 12)

#### PCoA ####
library(vegan)
library(ade4)
library(ggplot2)
library(RColorBrewer)

tab.dist<-vegdist(tab,method='euclidean')#基于euclidean距离

pcoa<- dudi.pco(tab.dist, scan = FALSE,nf=3)

#坐标轴解释量（前两轴）
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

#提取样本点坐标（前两轴）
sample_site <- data.frame({pcoa$li})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')

colnames(sample_site)[3] <- 'Sample_ID'

sample_site <- merge(sample_site,tmp_df)

colnames(sample_site)[4] <- 'level'

#以最终成绩作为分组
sample_site$level <- factor(sample_site$level)

library(ggplot2)

pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2,color=level)) +
  theme_classic()+#去掉背景框
  geom_point(size = 3.5)+  #可在这里修改点的透明度、大小
  scale_color_manual(values = brewer.pal(6,"Set1")) + #可在这里修改点的颜色
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank()
  )+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%')) 

pcoa_plot

ggsave(paste0('Result2_PCoA.pdf'),pcoa_plot,width = 8, height = 8)

#### Step 3 ROC ####

rm(list = ls())

setwd('/Volumes/HS_SSD\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step14_Nomogram')

library(tidyverse)
library(pheatmap)
library("pROC")

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

Clincal_BRCA <- read.csv('/Volumes/HS_SSD/Other_analysis/XL/Analysis_20200512/20200529_Analysis/Step1_Rawdata/BRCA_Clinical.csv',row.names = 1,na.strings = '')
Clincal_BRCA[1:5,1:5]

sampleInfo <- Clincal_BRCA[Clincal_BRCA$sample_type == 'Primary Tumor',]

rm(Clincal_BRCA)

sampleInfo <- na.omit(sampleInfo)

sampleInfo <- sampleInfo[sampleInfo$futime_os > 30, ]

sampleInfo <- sampleInfo[sampleInfo$Gender != 'MALE', ]

rownames(sampleInfo) <- substring(rownames(sampleInfo),1,12)

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step4_TNBC_WGCNA/Result4_tidyrawdata_exp.RData")

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group.RData")

sampleInfo <- sampleInfo[rownames(tmp_df), ]

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step14_Nomogram/Lasso_Cluster_feature.RData")

datExpr <- subset(datExpr,select = Lasso_Cluster_feature)

datExpr$sample <- rownames(datExpr)
sampleInfo$sample <- rownames(sampleInfo)

df <- merge(tmp_df,datExpr)
rownames(df) <- df$sample
df <- df[,-1]

df$Group <- factor(df$Group)

# df$Group <- ifelse(df$Group == 1, 'TNBC_ClusterA','TNBC_ClusterB')

mod1 <- glm(Group~., family = binomial(link = 'logit'),data = df)

index.min <- mod1$coefficients[-1] %>% as.numeric()

signature <- as.matrix(df[,c(-1)]) %*% as.matrix(index.min)
signature <- as.numeric(signature)

df$signature <- signature

library(pROC)

roc1 <- roc(df$Group,df$signature)

pdf('Result3_ROC.pdf',width = 12,height = 12, onefile=FALSE)

plot(roc1,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"),
     
     max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)

dev.off()

#### Step 4 Normal ####

rm(list = ls())

setwd('/Volumes/HS_SSD\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step14_Nomogram')

library(tidyverse)
library(pheatmap)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

Clincal_BRCA <- read.csv('/Volumes/HS_SSD/Other_analysis/XL/Analysis_20200512/20200529_Analysis/Step1_Rawdata/BRCA_Clinical.csv',row.names = 1,na.strings = '')
Clincal_BRCA[1:5,1:5]

sampleInfo <- Clincal_BRCA[Clincal_BRCA$sample_type == 'Primary Tumor',]

rm(Clincal_BRCA)

sampleInfo <- na.omit(sampleInfo)

sampleInfo <- sampleInfo[sampleInfo$futime_os > 30, ]

sampleInfo <- sampleInfo[sampleInfo$Gender != 'MALE', ]

rownames(sampleInfo) <- substring(rownames(sampleInfo),1,12)

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group_5gene.RData")

sampleInfo <- sampleInfo[tmp_df$sample, ]

df <- cbind(tmp_df,sampleInfo)

df <- df[,-2]

df$Group <- ifelse(df$Group == 1,'TNBC_ClusterA','TNBC_ClusterB')

library(ggplot2)
library(survminer)
library(survival)

fit <- survfit(Surv(futime_os, fustatus_os)~Group,data = df)

p <- ggsurvplot(fit, conf.int=F, pval=T, risk.table=T,
                legend.title='Risk Score', palette = 'jco',
                risk.table.height = 0.3)
p

df <- df[,-c(4:7,12,13)]

library(dplyr)

df_nom <- df %>% dplyr::select(fustatus_os,futime_os,everything())

colnames(df_nom)[1] <- 'OS_status'
colnames(df_nom)[2] <- 'OS_time'

table(df_nom$AJCC_Status)

df_nom$AJCC_Status[which(df_nom$AJCC_Status == "Stage I" | df_nom$AJCC_Status == "Stage IA" | df_nom$AJCC_Status == "Stage IB")] <- 'StageI'
df_nom$AJCC_Status[which(df_nom$AJCC_Status == "Stage IIA" | df_nom$AJCC_Status == "Stage IIB"| df_nom$AJCC_Status == "Stage II")] <- 'StageII'
df_nom$AJCC_Status[which(df_nom$AJCC_Status == "Stage IIIA" | df_nom$AJCC_Status == "Stage IIIB" | df_nom$AJCC_Status == "Stage IIIC" | df_nom$AJCC_Status == "Stage III")] <- 'StageIII'
df_nom$AJCC_Status[which(df_nom$AJCC_Status == "Stage IV")] <- 'StageIV'
df_nom$AJCC_Status[which(df_nom$AJCC_Status == "Stage X")] <- ''
table(df_nom$AJCC_Status)

df_nom$AJCC_Status <- factor(df_nom$AJCC_Status)
df_nom$PAM50 <- factor(df_nom$PAM50)
df_nom$Stage_M <- factor(df_nom$Stage_M)
df_nom$Stage_N <- factor(df_nom$Stage_N)
df_nom$Stage_T <- factor(df_nom$Stage_T)
df_nom$Group <- factor(df_nom$Group)

for (i in 1:ncol(df_nom)) {
  
  df_nom[,i] <- as.numeric(df_nom[,i])
  
}

BaSurv <- Surv(time = df_nom$OS_time,event = df_nom$OS_status)

## 构建一个单因素Cox分析Function，便于批量操作时调用
UniCox <- function(x){
  FML <- as.formula(paste0('BaSurv~',x))
  GCox <- coxph(FML, data = df_nom)
  GSum <- summary(GCox)
  HR <- round(GSum$coefficients[,2],2)
  PValue <- round(GSum$coefficients[,5],3)
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-")
  Unicox <- data.frame("characteristics" = x,
                       "Hazard Ratio" = HR,
                       "CI95" = CI,
                       "P Value" = PValue)
  return(Unicox)
}


UniCox(colnames(df_nom)[5])

library(plyr)

VarNames <- colnames(df_nom)[3:ncol(df_nom)] ## 选择需要纳入单因素Cox回归的变量
UniVar <- lapply(VarNames,UniCox) ## 批量单因素Cox分析
UniVar <- ldply(UniVar,data.frame) ## 将结果整理成数据框

## 选择单因素Cox结果中，P值＜0.2的变量纳入多因素Cox分析。
(GetFactors_uni <- UniVar$characteristics[which(UniVar$P.Value < 0.05)] %>% as.character() )
## 多因素分析 ##
fml <-  as.formula(paste0('BaSurv~',paste0(GetFactors_uni,collapse = '+')))
MultiCox <- coxph(fml, data = df_nom)

## 交互作用
library(car)
vif(MultiCox,digits = 3)

MultiSum <- summary(MultiCox)

MultiName <- as.character(GetFactors_uni)
MHR <- round(MultiSum$coefficients[,2],2)
MPV <- round(MultiSum$coefficients[,5],3)
MCIL <- round(MultiSum$conf.int[,3],2)
MCIU <- round(MultiSum$conf.int[,4],2)
MCI <- paste0(MCIL,'-',MCIU)
MulCox <- data.frame('characteristics' = MultiName,
                     'Hazard Ratio' = MHR,
                     'CI95' = MCI,
                     'P Value' = MPV)

## merge data frame ##
Final <- merge.data.frame(UniVar,MulCox, by = "characteristics", all = T, sort = T)
View(Final)
write.csv(Final,'Result4_CoxReg_new.csv')

GetFactors_mul <- Final$characteristics[which(Final$P.Value.y < 0.05)] ## 确定纳入预测模型的因子（多因素Cox结果，P值小于0.05）


library(dplyr)

df_nom <- df %>% dplyr::select(fustatus_os,futime_os,everything())

colnames(df_nom)[1] <- 'OS_status'
colnames(df_nom)[2] <- 'OS_time'

table(df_nom$AJCC_Status)

df_nom$AJCC_Status[which(df_nom$AJCC_Status == "Stage I" | df_nom$AJCC_Status == "Stage IA" | df_nom$AJCC_Status == "Stage IB")] <- 'StageI'
df_nom$AJCC_Status[which(df_nom$AJCC_Status == "Stage IIA" | df_nom$AJCC_Status == "Stage IIB"| df_nom$AJCC_Status == "Stage II")] <- 'StageII'
df_nom$AJCC_Status[which(df_nom$AJCC_Status == "Stage IIIA" | df_nom$AJCC_Status == "Stage IIIB" | df_nom$AJCC_Status == "Stage IIIC" | df_nom$AJCC_Status == "Stage III")] <- 'StageIII'
df_nom$AJCC_Status[which(df_nom$AJCC_Status == "Stage IV")] <- 'StageIV'
df_nom$AJCC_Status[which(df_nom$AJCC_Status == "Stage X")] <- ''
table(df_nom$AJCC_Status)

library(rms)

bc <- na.omit(df_nom) ## 去掉NA

dd <- datadist(bc) ## 封装数据

options(datadist="dd") ## 全局变量设置

BaSurv <- Surv(time = bc$OS_time,event = bc$OS_status)

fml <- as.formula(paste0('BaSurv~',paste0(GetFactors_mul,collapse = '+'))) ## 构建多因素Cox公式

f <- cph(fml, x=T, y=T, surv=T, data=bc)

surv <- Survival(f)

nom <- nomogram(f, ##生存函数
                
                fun=list(function(x) surv(365, x), ## 我们想预测患者的1年、3年和5年的生存率
                         function(x) surv(365*3, x), 
                         function(x) surv(365*5, x)), 
                
                lp=F, funlabel=c("1-year survival","3-year survival", "5-year survival"), ## label注释
                
                maxscale=100, ## 评分最大100分
                
                fun.at=seq(0.1,0.9,0.1)) ## 设置刻度

## Nomogram 可视化 ##

pdf(paste0('Plot5_Result_normgram.pdf'),width = 24,height = 12, onefile=FALSE)

plot(nom,cex.var = 2,cex.axis = 1.5,lwd = 10,xfrac = 0.5,tcl = 0.5)

dev.off()


##

library(rpart)
library(rpart.plot)
library(tidyverse)
library(survival)

ct <- rpart.control(xval=10, minbucket=30, cp=0.01)  
fit <- rpart(fml,data = bc,method = 'exp',control = ct)

print(fit)

pdf('Plot12_emt_tree_os.pdf',width = 12,height = 12)
rpart.plot(fit,split.col="black",main='Decision Tree', type=4)
dev.off()

#### Chr 5 Discrimination & Calibration ####

## C index 指数 ##
##注意，我们这里是通过一个模型来预测1年、3年和5年的生存率，有的审稿人不太明白会提出各种sx问题

validate(f, method="boot", B=1000, dxy=T)

c_index <- rcorrcens(Surv(OS_time, OS_status) ~ predict(f), data = bc)

index <- 1- c_index[1] %>% round(.,3) ## 提取C-index数值，保留小数后3位

low95CI <- c(index - c_index[4]/2) %>% round(.,3)

up95CI <- c(index + c_index[4]/2) %>% round(.,3)

## 输出C-index 结果 ##
sink(paste0('Chr4_Result_C-index.txt'))

(Cindex_df <- data.frame(c_index = index,low95CI = low95CI,up95CI = up95CI))

sink()

## 1年、3年、5年的校准度分析 ##

f1 <- cph(fml, x=T, y=T, surv=T, data=bc, time.inc=365)

cal1 <- calibrate(f1, 
                  cmethod="KM", 
                  method="boot", 
                  u=365, ## 预测时间
                  m=nrow(bc)/3, ## 界值选择
                  B=1000)

pdf(paste0('Plot6 calibration 1 year cal.pdf'),width = 12,height = 12, onefile=FALSE)

plot(cal1,
     lwd = 4,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("#2166AC"),#error bar的颜色
     # xlim = c(0,1),ylim= c(0,1),
     xlab="Nomogram-Predicted Probability of 1-Year OS",
     ylab="Actual 1-Year OS (proportion)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小
lines(cal1[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 4, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("#2166AC")) #连线的颜色
mtext("")
box(lwd = 2) #边框粗细
abline(0,1,lty = 3, #对角线为虚线
       lwd = 4, #对角线的粗细
       col = c("#224444")#对角线的颜色
) 
dev.off()


f3 <- cph(fml, x=T, y=T, surv=T, data=bc, time.inc=365*3)

cal3 <- calibrate(f3, 
                  cmethod="KM", 
                  method="boot",
                  u=365*3, 
                  m=nrow(bc)/3, 
                  B=1000)

pdf(paste0('Plot6 calibration 3 year cal.pdf'),width = 12,height = 12, onefile=FALSE)

plot(cal3,
     lwd = 4,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("#2166AC"),#error bar的颜色
     # xlim = c(0,1),ylim= c(0,1),
     xlab="Nomogram-Predicted Probability of 3-Year OS",
     ylab="Actual 3-Year OS (proportion)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小
lines(cal3[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 4, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("#2166AC")) #连线的颜色
mtext("")
box(lwd = 2) #边框粗细
abline(0,1,lty = 3, #对角线为虚线
       lwd = 4, #对角线的粗细
       col = c("#224444")#对角线的颜色
) 
dev.off()

f5 <- cph(fml, x=T, y=T, surv=T, data=bc, time.inc=365*5)
cal5 <- calibrate(f5, cmethod="KM", method="boot", u=365*5, m=nrow(bc)/3, B=1000)

pdf(paste0('Plot6 calibration 5 year cal.pdf'),width = 12,height = 12, onefile=FALSE)

plot(cal5,
     lwd = 4,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("#2166AC"),#error bar的颜色
     # xlim = c(0,1),ylim= c(0,1),
     xlab="Nomogram-Predicted Probability of 5-Year OS",
     ylab="Actual 5-Year OS (proportion)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小
lines(cal5[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 4, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("#2166AC")) #连线的颜色
mtext("")
box(lwd = 2) #边框粗细
abline(0,1,lty = 3, #对角线为虚线
       lwd = 4, #对角线的粗细
       col = c("#224444")#对角线的颜色
) 
dev.off()

#### Chr 6 Time-ROC进行验证 ####
library(timeROC)

cox.timp2 <- coxph(Surv(bc$OS_time, bc$OS_status == 1) ~ bc$Group)

lpFit <- cox.timp2$linear.predictors

roc.fit <-timeROC(T = bc$OS_time, # 结局时间
                  
                  delta = bc$OS_status, # 生存结局
                  
                  marker = lpFit, # 预测变量
                  
                  cause = 1, # 阳性结局赋值，比如死亡，复发的赋值
                  
                  weighting = "marginal", ## 'marginal' = Kaplan-Meier estimator of the censoring distribution 
                  ## or 'cox' = censoring by the Cox model
                  
                  times = c(365*c(1,3,5,7,10)), ## 时间点，选取1年、3年和5年生存率
                  
                  ROC = T,
                  
                  iid = TRUE) 

library(tidyverse)

CI95 <- confint(roc.fit)[1] %>% unlist()

ROC_1 <- paste0(' : ',roc.fit$AUC[1] %>% as.numeric() %>% round(.,2),'(',CI95[1],'-',CI95[length(CI95)/2 + 1],')')
ROC_3 <- paste0(' : ',roc.fit$AUC[2] %>% as.numeric() %>% round(.,2),'(',CI95[2],'-',CI95[length(CI95)/2 + 2],')')
ROC_5 <- paste0(' : ',roc.fit$AUC[3] %>% as.numeric() %>% round(.,2),'(',CI95[3],'-',CI95[length(CI95)/2 + 3],')')
ROC_7 <- paste0(' : ',roc.fit$AUC[4] %>% as.numeric() %>% round(.,2),'(',CI95[4],'-',CI95[length(CI95)/2 + 4],')')
ROC_10 <- paste0(' : ',roc.fit$AUC[5] %>% as.numeric() %>% round(.,2),'(',CI95[5],'-',CI95[length(CI95)/2 + 5],')')

col_all <- rainbow(10)

filename <- paste0(paste0('Plot7_OS_AUC_t_group_TypeI.pdf'))

pdf(filename,onefile = F)

plot(roc.fit,time=365,col = col_all[1],add =FALSE,title = F)#time 是时间点，col是线条颜色，add指是否添加在上一张图中

plot(roc.fit,time=365*3,col = col_all[2],add = TRUE)

plot(roc.fit,time=365*5,col = col_all[3],add = TRUE)

plot(roc.fit,time=365*7,col = col_all[4],add = TRUE)

plot(roc.fit,time=365*10,col = col_all[5],add = TRUE)

title(main = 'Time-dependent ROC curve')

legend("bottomright",
       
       legend = c(paste0('AUC at 1 year',ROC_1), 
                  paste0('AUC at 3 years',ROC_3),
                  paste0('AUC at 5 years',ROC_5),
                  paste0('AUC at 7 years',ROC_7),
                  paste0('AUC at 10 years',ROC_10)),
       
       lty = c("solid","solid","solid","solid","solid"),
       col = col_all[1:5],
       bty = 'n')

dev.off()

#### Chr 7 DCA 临床决策曲线 ####

bc <- na.omit(bc) ## 去掉NA

dd <- datadist(bc) ## 封装数据

BaSurv <- Surv(time = bc$OS_time,event = bc$OS_status)

fml <- as.formula(paste0('BaSurv~',paste0(GetFactors_mul,collapse = '+'))) ## 构建多因素Cox公式

#Run the cox model
coxmod = coxph(fml, data=bc)

source("/Volumes/HS_SSD/Other_analysis/JLX_Nomogram/Lesson3/Lesson3_Rfunction/stdca.R")

#the probability of failure is calculated by subtracting the probability of 

#survival from 1. 

bc$Risk_Signature = c(1- (summary(survfit(coxmod,newdata=bc), times=365)$surv))

#Run the decision curve analysis (with a smoother)

km <- stdca(data=bc, ## 数据集
            outcome="OS_status", ## 生存状态
            ttoutcome="OS_time", ## 生存时间
            timepoint=365, ## 预测的事件，1年（365天）
            predictors="Risk_Signature", 
            xstop=0.4,##
            ymin = -0.01,
            loess.span=0.50, 
            smooth=TRUE)

#Plotting the curves

pdf(paste0('Plot8_Result_DCA_Single.pdf'),width = 12,height = 8, onefile=FALSE)

plot(km$net.benefit.threshold, km$net.benefit.none, type = "l", lwd=2, 
     
     xlim=c(0,.5), ylim=c(-.02, .07), ## x轴和y轴的宽度
     
     xlab = "Threshold Probability", ylab = "Net Benefit") ## abel

lines(km$net.benefit$threshold, km$net.benefit$all, type="l", col=3, lwd=2) 

lines(km$net.benefit$threshold, km$net.benefit$Risk_Signature, type="l", col=10 , lwd = 4,lty = 2)

abline(h = 0, lwd =2,col  ="black")

legend("topright", cex=1, legend=c("None", "All", "Risk Signature"), 
       
       col=c(17, 3, 10), ## 线的颜色
       
       lwd=c(1.5, 1.5, 2.5), ## 线的宽度
       
       lty=c(1, 1, 2), ## 线的type
       
       bty = "n") ## 不要边框

dev.off()

#### 5 gene in Cluster Boxplot ####

rm(list = ls())

setwd('/Volumes/HS_SSD\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step14_Nomogram')

library(tidyverse)
library(pheatmap)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

Clincal_BRCA <- read.csv('/Volumes/HS_SSD/Other_analysis/XL/Analysis_20200512/20200529_Analysis/Step1_Rawdata/BRCA_Clinical.csv',row.names = 1,na.strings = '')
Clincal_BRCA[1:5,1:5]

sampleInfo <- Clincal_BRCA[Clincal_BRCA$sample_type == 'Primary Tumor',]

rm(Clincal_BRCA)

sampleInfo <- na.omit(sampleInfo)

sampleInfo <- sampleInfo[sampleInfo$futime_os > 30, ]

sampleInfo <- sampleInfo[sampleInfo$Gender != 'MALE', ]

rownames(sampleInfo) <- substring(rownames(sampleInfo),1,12)

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step4_TNBC_WGCNA/Result4_tidyrawdata_exp.RData")

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group_5gene.RData")

sampleInfo <- sampleInfo[rownames(tmp_df), ]

load("/Volumes/HS_SSD/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step14_Nomogram/Lasso_Cluster_feature.RData")

datExpr <- subset(datExpr,select = Lasso_Cluster_feature)

datExpr$sample <- rownames(datExpr)
sampleInfo$sample <- rownames(sampleInfo)

df <- merge(tmp_df,datExpr)
rownames(df) <- df$sample
df <- df[,-1]

df$Group <- ifelse(df$Group == 1,'TNBC_ClusterA','TNBC_ClusterB')

library(reshape2)

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

ggsave('Gene_Expression_boxplot.pdf',p3,width = 24,height = 12)

























