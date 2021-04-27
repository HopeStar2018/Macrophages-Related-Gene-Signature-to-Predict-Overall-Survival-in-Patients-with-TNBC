rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr3_GSE103091_Cibersort')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr1_GSE103091_reanno/GSE103091_exprdf_uniq.RData")

library(parallel)

library(doParallel)

source("CIBERSORT.R")

out=rbind(ID=colnames(GSE103091_exprdf_uniq),GSE103091_exprdf_uniq)

write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)

stopCluster(cl)

#### Method II ####

library(affy)
library(annotate)
library(hgu133plus2hsentrezgcdf)
library(org.Hs.eg.db)

setwd('G:\\GEO_database\\TNBC\\GSE103091_RAW')

Data<-ReadAffy(cdfname = "hgu133plus2hsentrezgcdf")
eset<-rma(Data)
# eset<-mas5(Data)
ID<-featureNames(eset)
ID2<-sub("_at","",ID)
GS <- as.matrix(getSYMBOL(ID2, 'org.Hs.eg'))
ematrix<-exprs(eset)
rows <- GS
cols =  c("GeneSymbol",colnames(ematrix))
ematrix <- cbind(rows,ematrix)
ematrix <- ematrix[which(ematrix[,1] != "NA"),] #remove NAs
ematrix <- ematrix[order(ematrix[,1]),] #sort by gene name 
ematrix <- rbind(cols, ematrix)
write.table(ematrix,file="NormalizedExpressionArray.customCDF.txt",sep="\t", col.names=F, row.names=F,quote=FALSE)




















