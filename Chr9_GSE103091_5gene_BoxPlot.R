rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr9_GSE103091_5gene_BoxPlot')

load("H:/PHD_database/Paper_writing/A_Paper_XIII_PP_TNBC/Step1_TNBC_ssGSEA/Chr3_download_cli/GSE103091_cli_OS.RData")

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr1_GSE103091_reanno/GSE103091_exprdf_uniq.RData")

GSE103091_exprdf_uniq[1:6,1:6]
GSE103091_cli_OS[1:6,1:5]

GSE103091_cli_OS <- GSE103091_cli_OS[,-1]

GSE103091_cli_OS$Age <- as.numeric(GSE103091_cli_OS$Age)
GSE103091_cli_OS$fustat <- as.numeric(GSE103091_cli_OS$fustat)
GSE103091_cli_OS$Metastasis <- as.numeric(GSE103091_cli_OS$Metastasis)
GSE103091_cli_OS$futime <- as.numeric(GSE103091_cli_OS$futime)

GSE103091_cli_OS <- GSE103091_cli_OS[GSE103091_cli_OS$futime > 30, ]

data <- subset(GSE103091_exprdf_uniq,select = rownames(GSE103091_cli_OS))

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step14_Nomogram/Lasso_Cluster_feature.RData")

data <- data[Lasso_Cluster_feature, ]

data <- log2(data + 1)

data <- data.frame(t(data))

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
                               title=title,clusterAlg=method,distance="euclidean",seed=20,plot="png")

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

save(tmp_df,file = 'GSE103091_tmp_group_5gene.RData')

#### Step 2 KM analysis ####

GSE103091_cli_OS$sample <- rownames(GSE103091_cli_OS)

rt <- merge(tmp_df,GSE103091_cli_OS)

rm(tmp_df)
rm(GSE103091_cli_OS)

library(survival)
library("survminer")

my.surv <- Surv(rt$futime, rt$fustat)

fit <- survfit(my.surv~Group,data = rt)

tiffFile=paste0('GSE103091',"_SurvivalOS_5gene.pdf")

p <- ggsurvplot(fit, conf.int=F, pval=T, risk.table=T, legend.labs = c("ClusterC","ClusterD"), 
                legend.title='GSE103091 Cluster', palette = 'jco',
                risk.table.height = 0.3,axes.offset = F)
p

pdf(file = tiffFile,width = 8,height = 8,onefile = F)
print(p)
dev.off()

######################

rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr9_GSE103091_5gene_BoxPlot')

load("H:/PHD_database/Paper_writing/A_Paper_XIII_PP_TNBC/Step1_TNBC_ssGSEA/Chr3_download_cli/GSE103091_cli_OS.RData")

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr1_GSE103091_reanno/GSE103091_exprdf_uniq.RData")

GSE103091_exprdf_uniq[1:6,1:6]
GSE103091_cli_OS[1:6,1:5]

GSE103091_cli_OS <- GSE103091_cli_OS[,-1]

GSE103091_cli_OS$Age <- as.numeric(GSE103091_cli_OS$Age)
GSE103091_cli_OS$fustat <- as.numeric(GSE103091_cli_OS$fustat)
GSE103091_cli_OS$Metastasis <- as.numeric(GSE103091_cli_OS$Metastasis)
GSE103091_cli_OS$futime <- as.numeric(GSE103091_cli_OS$futime)

GSE103091_cli_OS <- GSE103091_cli_OS[GSE103091_cli_OS$futime > 30, ]

data <- subset(GSE103091_exprdf_uniq,select = rownames(GSE103091_cli_OS))

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step14_Nomogram/Lasso_Cluster_feature.RData")

data <- data[Lasso_Cluster_feature, ]

data <- log2(data + 1)

data <- data.frame(t(data))

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr9_GSE103091_5gene_BoxPlot/GSE103091_tmp_group_5gene.RData")

library(reshape2)

df <- cbind(tmp_df,data)

rownames(df) <- df$sample

df <- df[,-2]

df$Group <- ifelse(df$Group == 1,'TNBC_ClusterC','TNBC_ClusterD')

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

ggsave('GSE103091_Gene_Expression_boxplot.pdf',p3,width = 24,height = 12)
















