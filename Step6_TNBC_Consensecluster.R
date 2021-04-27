rm(list = ls())

library(survival)

library(survminer)

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step6_TNBC_Consensecluster')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step4_TNBC_WGCNA/Result4_tidyrawdata_exp.RData")

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step4_TNBC_WGCNA/TCGA_BRCA_CIBERSORT.ABS/geneInfo.RData")

color_pick <- 'blue'

hub_gene <- geneInfo[geneInfo$moduleColor == color_pick, ]
# hub_gene <- hub_gene$geneSymbol %>% as.character()
# hub_gene <- hub_gene$geneSymbol[hub_gene$GS.Macrophage.M1 > 0.2 & hub_gene$MM.blue > 0.8]
hub_gene <- hub_gene$geneSymbol[hub_gene$p.GS.Macrophage.M1 < 0.001]

#### uni survival ####

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN/Step2_ssGSEA_Score/Step2_TNBC_datTraits.RData")

TNBC_datTraits$sample <- rownames(TNBC_datTraits)

datExpr$sample <- rownames(datExpr)

df <- merge(TNBC_datTraits,datExpr)

df[1:6,1:15]

rownames(df) <- df$sample

df <- df[,-c(1:10)]

BaSurv <- Surv(time = df$futime_os, ## 生存时间
               event = df$fustatus_os) ## 生存状态

UniCox <- function(x){ ## 构建一个R function 便于后期调用
  FML <- as.formula(paste0('BaSurv~',x)) ## 构建生存分析公式
  GCox <- coxph(FML, data = df) ## Cox分析
  GSum <- summary(GCox) ## 输出结果
  HR <- round(GSum$coefficients[,2],2) ## 输出HR值
  PValue <- round(GSum$coefficients[,5],3) ## 输出P值
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-") ## 输出HR的执行区间
  Unicox <- data.frame("characteristics" = x, ## 返回结果，并构建数据框
                       "Hazard Ratio" = HR,
                       "CI95" = CI,
                       "P Value" = PValue)
  return(Unicox)
}

VarNames <- hub_gene ## 输出需要分析的变量名字
UniVar <- lapply(VarNames,UniCox) ## 批量做Cox分析

library(plyr)

UniVar <- ldply(UniVar,data.frame) ## 将结果整理为一个数据框

## 筛选其中P值<0.0001的变量纳入多因素cox分析。
(GetFactors <- UniVar$characteristics[which(UniVar$P.Value < 0.05)] %>% as.character()) 

# save(GetFactors,file = 'Result6_GetFactors.Rdata')

#######

# load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step14_Nomogram/Lasso_Cluster_feature.RData")


data <- subset(datExpr,select = GetFactors)

naSum=function(dt){return(sum(is.na(dt)))}

sum(apply(data, 1, naSum)==0)
samples=which(apply(data, 2, naSum)==0)
data=data[,samples]
# sampleInfo=sampleInfo[samples,]
length(data[,1])

std.data.by.row <- scale(t(data), scale=T)

library(limma)
library(ConsensusClusterPlus)
library(modeest)

method_all <- c("hc",'km','pam')

(method <- method_all[2])

pathway <- getwd()

dir.create(paste0('limma.dir.',method))
title=paste0(pathway,'/limma.dir.',method)
results = ConsensusClusterPlus(as.matrix(std.data.by.row),maxK=6,reps=1000,pItem=0.8,pFeature=1,
                               title=title,clusterAlg=method,distance="euclidean",seed=1000,plot="png")

save(results,file = paste0(title,'/',method,'.RData'))

maxK <- 6
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)}#end for i# The optimal K
(Kvec[which.min(PAC)])

clsCut=2
Iteration <- 0;
k3_run1 <- results[[clsCut]][["consensusClass"]]
table(k3_run1)

tmp_df <- data.frame(k3_run1)
tmp_df$sample <- rownames(tmp_df)

colnames(tmp_df)[1] <- 'Group'

save(tmp_df,file = 'tmp_group_5gene.RData')

#### Step 2 KM analysis ####
rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step6_TNBC_Consensecluster')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN/Step2_ssGSEA_Score/Step2_TNBC_datTraits.RData")

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group_5gene.RData")

TNBC_datTraits$sample <- rownames(TNBC_datTraits)

rt <- merge(tmp_df,TNBC_datTraits)

rm(tmp_df)
rm(TNBC_datTraits)

library(survival)
library("survminer")

my.surv <- Surv(rt$futime_os, rt$fustatus_os)

fit <- survfit(my.surv~Group,data = rt)

tiffFile=paste0('Blue_Module',"_SurvivalOS_5gene.pdf")

p <- ggsurvplot(fit, conf.int=F, pval=T, risk.table=T, legend.labs = c("TNBC_ClusterA","TNBC_ClusterB"), 
                legend.title='Risk Factor', palette = 'jco',
                risk.table.height = 0.3,axes.offset = F)
p

pdf(file = tiffFile,width = 8,height = 6,onefile = F)
print(p)
dev.off()













