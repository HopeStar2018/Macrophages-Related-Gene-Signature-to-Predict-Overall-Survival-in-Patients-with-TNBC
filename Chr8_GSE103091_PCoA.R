rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr8_GSE103091_PCoA')

#### 12 hub genes ####

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/Result6_GetFactors.Rdata")

GetFactors[4] <- "HLA-DQB2"
GetFactors[7] <- "HLA-DQA2"

#### 5 hub genes ####

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step14_Nomogram/Lasso_Cluster_feature.RData")

####

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr1_GSE103091_reanno/GSE103091_exprdf_uniq.RData")

load("H:/PHD_database/Paper_writing/A_Paper_XIII_PP_TNBC/Step1_TNBC_ssGSEA/Chr3_download_cli/GSE103091_cli_OS.RData")

GSE103091_exprdf_uniq <- subset(GSE103091_exprdf_uniq,select = rownames(GSE103091_cli_OS))

#### 12 hub genes ####

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr2_GSE103091_ConsenseCluster/tmp_group.RData")

#### 5 hub genes ####

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr2_GSE103091_ConsenseCluster/GSE103091_tmp_group_5genes.RData")

tab <- GSE103091_exprdf_uniq[Lasso_Cluster_feature, ]

tab <- data.frame(t(tab))

tab <- tab[rownames(tmp_df), ]

tmp_df$Group <- ifelse(tmp_df$Group == 1,'TNBC_ClusterC','TNBC_ClusterD')

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

ggsave('plsda_plot_5genes.pdf', plsda_plot, width = 12, height = 12)

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

ggsave(paste0('Plot4_PCoA_5genes.pdf'),pcoa_plot,width = 8, height = 8)

########

##

library(pca3d)

pca.results <- prcomp(tab, center = TRUE, scale. = FALSE)

pca3d(pca.results, group= sample_site$level)

# install.packages("devtools")
# library(devtools)
# install_github("GuangchuangYu/yyplot")
# install_github('fawda123/ggord')
library(ggplot2)
library(plyr)
library(ggord)
library(yyplot)

#调用yyplot包里的geom_ord_ellipse函数
source('./geom_ord_ellipse.R') #该文件位于当前文件夹

#用ggord画基本PCA图
pcoa_plot_pca <- ggord(pca.results, grp_in = sample_site$level, repel=TRUE,
                       ellipse = FALSE, #不显示置信区间背景色
                       size = 2, #样本的点大小
                       alpha=0.5, #设置点为半透明，出现叠加的效果
                       #如果用自定义的颜色，就运行下面这行
                       cols = mycol[1:length(unique(sample_site$level))],
                       arrow = NULL,txt = NULL) + #不画箭头和箭头上的文字
  theme(panel.grid =element_blank()) + #去除网格线
  
  #用yyplot添加置信区间圆圈
  geom_ord_ellipse(ellipse_pro = .95, #设置置信区间
                   size=1.5, #线的粗细
                   lty=1 ) #实线

pcoa_plot_pca

ggsave(paste0('Plot4_typeII_PCA.pdf'),pcoa_plot_pca,width = 8, height = 8)




