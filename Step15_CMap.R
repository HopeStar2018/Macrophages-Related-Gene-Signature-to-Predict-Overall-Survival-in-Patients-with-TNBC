rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step15_CMap')

load("G:/TCGA_database/TCGA_BRCA/TCGA_BRCA/TCGA_BRCA_mRNA_count.Rdata")

rawcount <- data.frame(rawcount)

colnames(rawcount)[1:5]

group=sapply(strsplit(colnames(rawcount),"\\."),"[",4)

TCGA_BRCA_Ca_Count <- rawcount[,group == '01A']
TCGA_BRCA_Normal_Count <- rawcount[,group == '11A']

load("H:/PHD_database/Commonly_used_statistical_analysis/AAA_used_RData_Function_code/Annotation_GPL/AnnoTCGA.RData")

TCGA_BRCA_Ca_Count$ENSG <- rownames(TCGA_BRCA_Ca_Count)
TCGA_BRCA_Normal_Count$ENSG <- rownames(TCGA_BRCA_Normal_Count)

TCGA_BRCA_Ca_Count <- merge(AnnoTCGA,TCGA_BRCA_Ca_Count)
TCGA_BRCA_Normal_Count <- merge(AnnoTCGA,TCGA_BRCA_Normal_Count)

colnames(TCGA_BRCA_Ca_Count)[1:6]

TCGA_BRCA_Ca_Count <- TCGA_BRCA_Ca_Count[,-c(1,3)]
TCGA_BRCA_Normal_Count <- TCGA_BRCA_Normal_Count[,-c(1,3)]

TCGA_BRCA_Ca_Count <- aggregate(.~ GeneName,TCGA_BRCA_Ca_Count,median)
TCGA_BRCA_Normal_Count <- aggregate(.~ GeneName,TCGA_BRCA_Normal_Count,median)

# save(TCGA_BRCA_Ca_Count,TCGA_BRCA_Normal_Count,file = 'TCGA_BRCA_Ca_Count_genes.RData')

rm(list = ls())

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step15_CMap/TCGA_BRCA_Ca_Count_genes.RData")

rownames(TCGA_BRCA_Ca_Count) <- TCGA_BRCA_Ca_Count$GeneName
rownames(TCGA_BRCA_Normal_Count) <- TCGA_BRCA_Normal_Count$GeneName

TCGA_BRCA_Ca_Count[1:5,1:5]
TCGA_BRCA_Ca_Count <- TCGA_BRCA_Ca_Count[,-1]
TCGA_BRCA_Normal_Count <- TCGA_BRCA_Normal_Count[,-1]

rm(group)
rm(AnnoTCGA)
rm(rawcount)

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/tmp_group_5gene.RData")

tmp_name <- tmp_df$sample[tmp_df$Group == 1]

rm(tmp_df)
 
colnames(TCGA_BRCA_Ca_Count) <- gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*","\\1\\-\\2\\-\\3",colnames(TCGA_BRCA_Ca_Count))

TCGA_BRCA_Ca_Count <- subset(TCGA_BRCA_Ca_Count,select = tmp_name)

rm(tmp_name)

exp <- cbind(TCGA_BRCA_Ca_Count,TCGA_BRCA_Normal_Count)

#### EdgeR ####

library("edgeR")

foldChange=1

padj=0.05

dimnames=list(rownames(exp),colnames(exp))

data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

data=avereps(data)

data=data[rowMeans(data)>1,]

group=c(rep("TNBC_ClusterA",ncol(TCGA_BRCA_Ca_Count)),rep("Normal",ncol(TCGA_BRCA_Normal_Count)))

design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("TNBC_ClusterA","Normal"))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]

save(diff,diffSig,file = 'TNBC_ClusterA&Normal_edgeR.RData')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step15_CMap/TNBC_ClusterA&Normal_edgeR.RData")

foldChange=1

padj=0.05

diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]

diffUp <- diffUp[order(diffUp$logFC,decreasing = T),]

diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]

diffUp_gene <- rownames(diffUp)[1:2000]

diffDown_gene <- rownames(diffDown)

write.table(diffUp_gene,file = 'String_diffup.txt', quote = F, col.names = F,row.names = F)

write.table(diffDown_gene,file = 'String_diffDown.txt', quote = F, col.names = F,row.names = F)

#### CMap 2.0 ####

load("H:/PHD_database/Commonly_used_statistical_analysis/AAA_used_RData_Function_code/Annotation_GPL/anno_GPL96.RData")

head(anno)

diffUp$Gene.Symbol <- rownames(diffUp)

ID_diffUp <- merge(anno,diffUp)

diffUp_ID <- ID_diffUp[order(ID_diffUp$logFC,decreasing = T), ]

diffUp_ID <- diffUp_ID$ID[1:500]


diffDown$Gene.Symbol <- rownames(diffDown)

ID_diffDown <- merge(anno,diffDown)

diffDown_ID <- ID_diffDown[order(ID_diffDown$logFC,decreasing = T), ]

diffDown_ID <- diffDown_ID$ID[1:500]

write.table(diffUp_ID, paste0("TNBC_ID_up500.grp"),
            
            row.names = F, sep = "\t", quote = F, col.names = F)

write.table(diffDown_ID, paste0("TNBC_ID_down500.grp"),
            
            row.names = F, sep = "\t", quote = F, col.names = F)


#### ico ####

diffUp_CMap <- diffUp[order(diffUp$logFC,decreasing = T), ]

diffUp_CMap <- rownames(diffUp_CMap)[1:150]

diffDown_CMap <- diffDown[order(diffDown$logFC), ]

diffDown_CMap <- rownames(diffDown_CMap)[1:150]

write.table(diffUp_CMap, paste0("TNBC_up150.grp"),
            
            row.names = F, sep = "\t", quote = F, col.names = F)

write.table(diffDown_CMap, paste0("TNBC_down150.grp"),
            
            row.names = F, sep = "\t", quote = F, col.names = F)

#### DEGs ####

foldChange=1

padj=0.05

library("DESeq")
library("limma")

dimnames=list(rownames(exp),colnames(exp))

data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

data=avereps(data)

data=data[rowMeans(data)>1,] 

data=round(data,0)

group=c(rep("TNBC_ClusterA",ncol(TCGA_BRCA_Ca_Count)),rep("Normal",ncol(TCGA_BRCA_Normal_Count)))

design = factor(group)
newTab = newCountDataSet( data, design )
newTab = estimateSizeFactors(newTab)
newData=counts(newTab, normalized=TRUE )

#have replicates
newTab = estimateDispersions( newTab, fitType = "local")
diff = nbinomTest( newTab, "TNBC_ClusterA", "Normal")
diff = diff[is.na(diff$padj)==FALSE,]
diff = diff[order(diff$pval),]

#########

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/Result6_GetFactors.Rdata")

diff_genes <- c(diffDown_CMap,diffUp_CMap)

intersect(GetFactors,diff_genes)

#######





































