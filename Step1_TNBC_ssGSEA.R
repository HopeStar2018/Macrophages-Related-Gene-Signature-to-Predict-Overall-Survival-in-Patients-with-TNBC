rm(list = ls())
setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step1_TNBC_ssGSEA')

#options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
library(GSVA)
library(TCGAbiolinks)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(rtracklayer)
library(SummarizedExperiment)
library(clusterProfiler)
library(RColorBrewer)
library(maftools)
library(circlize)
library(matrixStats)
library(GetoptLong)
library(GenomicRanges)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

#### Step 1 download TCGA-BRCA Counts ####

load("G:/TCGA_database/TCGA_BRCA/TCGA_BRCA/TCGA_BRCA_mRNA_count.Rdata")

rawcount[1:6,1:6]

##### Step 2 read count转TPM ####

expMatrix <- rawcount

rm(rawcount)

expMatrix <- data.frame(expMatrix)

expMatrix[1:5,1:5]

eff_length2 <-read.csv('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step1_TNBC_ssGSEA\\eff_length.csv', row.names = 1, header = T)

eff_length2$gene_id <- rownames(eff_length2)

rownames(eff_length2) <- do.call(rbind,strsplit(eff_length2$gene_id,'\\.'))[,1]

feature_ids <- rownames(expMatrix)

if (! all(feature_ids %in% rownames(eff_length2))){
  tbl <- table(feature_ids %in% rownames(eff_length2))
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]],tbl[[1]])
  warning(msg1)
  
} 

if (! identical(feature_ids, rownames(eff_length2))){
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", nrow(eff_length2), nrow(expMatrix))
  warning(msg2)
}

# trim the expression matrix and effetive gene length
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]
mm <- match(rownames(expMatrix), rownames(eff_length2))
eff_length2 <- eff_length2[mm, ]

if (identical(rownames(eff_length2), rownames(expMatrix))){
  print("GTF and expression matix now have the same gene and gene in same order")
}

x <- expMatrix / eff_length2$eff_length

expMatrix_tpm <- t( t(x) / colSums(x) ) * 1e6 

TCGA_BRCA_tpm_ensg <- expMatrix_tpm

TCGA_BRCA_tpm_ensg[1:6,1:6]

save(TCGA_BRCA_tpm_ensg,file = 'G:\\TCGA_database\\TCGA_BRCA\\BRCA_exp\\TPM\\TCGA_BRCA_tpm_ensg.Rdata')

#### Step 3 ssGSEA Immunity ####
rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step1_TNBC_ssGSEA')

load("G:/TCGA_database/TCGA_BRCA/BRCA_exp/TPM/TCGA_BRCA_tpm_ensg.Rdata")

gtf_v22 <-read.csv('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step1_TNBC_ssGSEA\\gtf_v22.csv', header = T)

tcga_expr <- data.frame(TCGA_BRCA_tpm_ensg)

rm(TCGA_BRCA_tpm_ensg)

library(data.table)

tcga_expr$gene_id <- rownames(tcga_expr)

tcga_expr_BRCA <- merge(gtf_v22,tcga_expr)

rm(gtf_v22)
rm(tcga_expr)

tcga_expr_BRCA[1:6,1:6]

tcga_expr_BRCA <- tcga_expr_BRCA[,-c(1,2,4)]

tcga_expr_BRCA <- aggregate(.~ENTREZID, tcga_expr_BRCA, median)

rownames(tcga_expr_BRCA)<- tcga_expr_BRCA$ENTREZID

tcga_expr_BRCA <- tcga_expr_BRCA[, -1]

tcga_expr_BRCA[1:3,1:2]

(load("easy_input_immunity.rdata"))

immunity[1]

tcga_gsva_BRCA <- as.data.frame(t(gsva(as.matrix(tcga_expr_BRCA), immunity, method = "ssgsea")))

str(tcga_gsva_BRCA)

library(pheatmap)

pheatmap(t(tcga_gsva_BRCA), scale = "row", show_colnames = F, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

save(tcga_gsva_BRCA,file = 'tcga_gsva_BRCA_ssGSEA.Rdata')

#### Step 4 mRNA & lncRNA ####

library(biomaRt)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl", #人
                   host = "www.ensembl.org") 

feature_info <- getBM(attributes = c("gene_biotype",
                                     #"transcript_biotype",#还可以提取transcript_biotype
                                     #如果基因名是gene symbol，就运行下面这行
                                     # "hgnc_symbol"), 
                      #如果基因名是ensembl ID，就运行下面这行
                      "ensembl_gene_id"),
                      #如果基因名是gene symbol，就运行下面这行
                      # filters = "hgnc_symbol", #小鼠是mgi_symbol，大鼠是mgi_symbol
                      #如果基因名是ensembl ID，就运行下面这行
                      filters = "ensembl_gene_id",
                      values = rownames(TCGA_BRCA_tpm_ensg), mart = ensembl)

#TCGA数据的ensembl ID跟最新版ensembl一致
if (nrow(TCGA_BRCA_tpm_ensg) != nrow(feature_info)){
  #查看哪些基因名不一致
  library(dplyr)
  diffName<-setdiff(rownames(TCGA_BRCA_tpm_ensg),feature_info[,2])
  length(diffName)
  head(diffName)
}

length(unique(feature_info$ensembl_gene_id))

unique(feature_info$gene_biotype)

#此处定义protein_coding作为mRNA
mRNA <-"protein_coding"
#根据实际研究目的，调整定义为lncRNA的gene_biotype，此处根据Vega定义如下8种biotype为lncRNA
lncRNA <- paste("lncRNA","non_coding","3prime_overlapping_ncRNA","antisense","lincRNA","sense_intronic","sense_overlapping","macro_lncRNA","bidirectional_promoter_lncRNA",sep = "|")
#还可以定义miRNA
miRNA <-"miRNA"

#下面就是表达矩阵里的mRNA、lncRNA、miRNA及其数量，把每类基因的基因名保存到相应的文件里。
mRNA.list<-feature_info[grepl(mRNA, feature_info$gene_biotype),]
write.table(mRNA.list,"mRNA.list.txt",quote = F,row.names = F, col.names = F)
nrow(mRNA.list)

lncRNA.list<-feature_info[grepl(lncRNA, feature_info$gene_biotype),]
write.table(lncRNA.list,"lncRNA.list.txt",quote = F,row.names = F, col.names = F)
nrow(lncRNA.list)

miRNA.list<-feature_info[grepl(miRNA, feature_info$gene_biotype),]
write.table(miRNA.list,"miRNA.list.txt",quote = F,row.names = F, col.names = F)
nrow(miRNA.list)

TCGA_BRCA_tpm_ensg <- data.frame(TCGA_BRCA_tpm_ensg)
TCGA_BRCA_tpm_ensg[1:6,1:6]

TCGA_BRCA_tpm_ensg_mRNA <- TCGA_BRCA_tpm_ensg[mRNA.list$ensembl_gene_id, ]
TCGA_BRCA_tpm_ensg_miRNA <- TCGA_BRCA_tpm_ensg[miRNA.list$ensembl_gene_id, ]
TCGA_BRCA_tpm_ensg_lncRNA <- TCGA_BRCA_tpm_ensg[lncRNA.list$ensembl_gene_id, ]

save(TCGA_BRCA_tpm_ensg_mRNA,file = 'G:\\TCGA_database\\TCGA_BRCA\\BRCA_exp\\TPM\\TCGA_BRCA_tpm_ensg_mRNA.Rdata')
save(TCGA_BRCA_tpm_ensg_miRNA,file = 'G:\\TCGA_database\\TCGA_BRCA\\BRCA_exp\\TPM\\TCGA_BRCA_tpm_ensg_miRNA.Rdata')
save(TCGA_BRCA_tpm_ensg_lncRNA,file = 'G:\\TCGA_database\\TCGA_BRCA\\BRCA_exp\\TPM\\TCGA_BRCA_tpm_ensg_lncRNA.Rdata')






















