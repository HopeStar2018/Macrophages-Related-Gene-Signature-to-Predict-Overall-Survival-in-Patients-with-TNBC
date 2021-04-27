rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step3_TNBC_Timer')

TCGA_BRCA_Timer <- read.csv('.\\infiltration_estimation_for_tcga.csv')

load("G:/TCGA_database/TCGA_BRCA/BRCA_exp/TPM/TCGA_BRCA_tpm_ensg.Rdata")

df <- data.frame(TCGA_BRCA_tpm_ensg)

rm(TCGA_BRCA_tpm_ensg)

colnames(df)
TCGA_BRCA_Timer$cell_type

BRCA_samples <- gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*","\\1\\-\\2\\-\\3\\-\\4",colnames(df))
BRCA_samples <- unique(BRCA_samples)
BRCA_samples <- substring(BRCA_samples,1,15)

TCGA_BRCA_Timer <- TCGA_BRCA_Timer[TCGA_BRCA_Timer$cell_type %in% BRCA_samples, ]
colnames(TCGA_BRCA_Timer)

rownames(TCGA_BRCA_Timer) <- TCGA_BRCA_Timer$cell_type %>% as.character()
TCGA_BRCA_Timer <- TCGA_BRCA_Timer[,-1]
colnames(TCGA_BRCA_Timer)

TCGA_BRCA_TIMER <- TCGA_BRCA_Timer[,c(1:6)]
TCGA_BRCA_CIBERSORT <- TCGA_BRCA_Timer[,c(7:28)]
TCGA_BRCA_CIBERSORT.ABS <- TCGA_BRCA_Timer[,c(29:50)]
TCGA_BRCA_QUANTISEQ <- TCGA_BRCA_Timer[,c(51:61)]
TCGA_BRCA_MCPCOUNTER <- TCGA_BRCA_Timer[,c(62:72)]
TCGA_BRCA_XCELL <- TCGA_BRCA_Timer[,c(73:111)]
TCGA_BRCA_EPIC <- TCGA_BRCA_Timer[,c(112:119)]

save(TCGA_BRCA_TIMER,TCGA_BRCA_CIBERSORT,TCGA_BRCA_CIBERSORT.ABS,TCGA_BRCA_QUANTISEQ,TCGA_BRCA_MCPCOUNTER,
     TCGA_BRCA_XCELL,TCGA_BRCA_EPIC,file = 'TCGA_BRCA_Immunity.RData')

############################################################################













































