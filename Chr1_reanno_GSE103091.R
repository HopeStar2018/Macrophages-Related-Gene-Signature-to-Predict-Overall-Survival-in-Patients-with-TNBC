rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr1_GSE103091_reanno')

load("H:\\PHD_database\\Paper_writing\\A_Paper_XIII_PP_TNBC\\Step1_TNBC_ssGSEA/Chr3_download_cli/GSE103091_cli_OS.RData")

library(tidyverse)

#### Step 0 load rawdata ####

library(affy)
library(limma)

setwd('G:\\GEO_database\\TNBC\\GSE103091_RAW')

# ## Method I ##
# 
# filenames <- dir('G:\\GEO_database\\TNBC\\GSE103091_RAW')
# 
# eset.rma <- justRMA(filenames=filenames, celfile.path='G:\\GEO_database\\TNBC\\GSE103091_RAW')
# 
# datExpr = exprs(eset.rma)

## Method II ##

mydata <- ReadAffy()

mydata1 <- rma(mydata)

# mydata2 <- mas5(mydata)

datExpr <- exprs(mydata1)

datExpr <- data.frame(datExpr)

datExpr[1:6,1:6]

load("H:/PHD_database/Commonly_used_statistical_analysis/AAA_used_RData_Function_code/Annotation_GPL/anno_GPL570.RData")

anno$Gene.Symbol <- as.character(anno$Gene.Symbol)

# anno[15165,2] <- 'PCDHGA1'
# anno[18494,2] <- 'PCDHGA1'
# anno[20431,2] <- 'PCDHGA1'
# anno[25130,2] <- 'PCDHGA1'

anno[21977,2] <- 'HLA-DQA2'

# save(anno,file = 'H:/PHD_database/Commonly_used_statistical_analysis/AAA_used_RData_Function_code/Annotation_GPL/anno_GPL570.RData')

datExpr$ID <- rownames(datExpr)

anno <- anno[anno$ID != '' ,]

df <- merge(anno,datExpr)

df[1:6,1:6]

df <- df[,-1]

temp <- sapply(df$Gene.Symbol, function(x){
  name <- unlist(strsplit(x %>% as.character(),'[ /// ]'))
  x <- name[1]
})

df$Gene.Symbol <- temp

df <- na.omit(df)

GSE103091_exprdf_uniq<-aggregate(.~Gene.Symbol,df,median)

rownames(GSE103091_exprdf_uniq) <- GSE103091_exprdf_uniq$Gene.Symbol %>% as.character()

GSE103091_exprdf_uniq <- GSE103091_exprdf_uniq[,-1]

GSE103091_exprdf_uniq[1:6,1:6]

colnames(GSE103091_exprdf_uniq) <- substring(colnames(GSE103091_exprdf_uniq),1,10)

save(GSE103091_exprdf_uniq,file = 'H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr1_GSE103091_reanno\\GSE103091_exprdf_uniq.RData')

#### Step 1 re-anno ####

rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr1_GSE103091_reanno')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr1_GSE103091_reanno/GSE103091_exprdf_uniq.RData")

GSE103091_exp <- GSE103091_exprdf_uniq

rm(GSE103091_exprdf_uniq)

rownames(GSE103091_exp)[1:4]

exprdf_uniq <- GSE103091_exp

library(biomaRt)

#用下面这行查看ensembl基因组版本跟host的对应关系
#listEnsemblArchives()

#运行下面两行，查看基因组
#mart = useMart('ensembl')
#listDatasets(mart)
#你需要哪个物种，就复制它在dataset列里的词，放在下面这行的`dataset = `参数里
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl", #人
                   host = "www.ensembl.org") 

feature_info <- getBM(attributes = c("gene_biotype",
                                     #"transcript_biotype",#还可以提取transcript_biotype
                                     #如果基因名是gene symbol，就运行下面这行
                                     "hgnc_symbol"), 
                      #如果基因名是ensembl ID，就运行下面这行
                      #"ensembl_gene_id"),
                      #如果基因名是gene symbol，就运行下面这行
                      filters = "hgnc_symbol", #小鼠是mgi_symbol，大鼠是mgi_symbol
                      #如果基因名是ensembl ID，就运行下面这行
                      #filters = "ensembl_gene_id",
                      values = rownames(exprdf_uniq), mart = ensembl)

#有些芯片注释的gene symbol跟最新版本ensembl的基因名不一致，需要返回上一步，换比较老的版本。
#TCGA数据的ensembl ID跟最新版ensembl一致
if (nrow(exprdf_uniq) != nrow(feature_info)){
  #查看哪些基因名不一致
  library(dplyr)
  diffName<-setdiff(rownames(exprdf_uniq),feature_info[,2])
  length(diffName)
  head(diffName)
}

length(unique(feature_info$hgnc_symbol))

table(feature_info$gene_biotype)

#此处定义protein_coding作为mRNA
mRNA <-"protein_coding"
#根据实际研究目的，调整定义为lncRNA的gene_biotype，此处根据Vega定义如下8种biotype为lncRNA
lncRNA <- paste('lncRNA',"non_coding","3prime_overlapping_ncRNA","antisense","lincRNA","sense_intronic","sense_overlapping","macro_lncRNA","bidirectional_promoter_lncRNA",sep = "|")
#还可以定义miRNA
miRNA <-"miRNA"

#下面就是表达矩阵里的mRNA、lncRNA、miRNA及其数量，把每类基因的基因名保存到相应的文件里。
mRNA.list<-feature_info[grepl(mRNA, feature_info$gene_biotype),]
# write.table(mRNA.list,"mRNA.list.txt",quote = F,row.names = F, col.names = F)
nrow(mRNA.list)

lncRNA.list<-feature_info[grepl(lncRNA, feature_info$gene_biotype),]
# write.table(lncRNA.list,"lncRNA.list.txt",quote = F,row.names = F, col.names = F)
nrow(lncRNA.list)

miRNA.list<-feature_info[grepl(miRNA, feature_info$gene_biotype),]
# write.table(miRNA.list,"miRNA.list.txt",quote = F,row.names = F, col.names = F)
nrow(miRNA.list)

GSE103091_mRNA <- exprdf_uniq[rownames(exprdf_uniq) %in% mRNA.list$hgnc_symbol, ]

GSE103091_miRNA <- exprdf_uniq[rownames(exprdf_uniq) %in% miRNA.list$hgnc_symbol, ]

GSE103091_lncRNA <- exprdf_uniq[rownames(exprdf_uniq) %in% lncRNA.list$hgnc_symbol, ]

save(GSE103091_mRNA,GSE103091_miRNA,GSE103091_lncRNA,file = 'GSE103091_expr_alltype.RData')

#####

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/Result6_GetFactors.Rdata")

GetFactors[4] <- "HLA-DQB2"
GetFactors[7] <- "HLA-DQA2"

intersect(rownames(GSE103091_exprdf_uniq),GetFactors)

















