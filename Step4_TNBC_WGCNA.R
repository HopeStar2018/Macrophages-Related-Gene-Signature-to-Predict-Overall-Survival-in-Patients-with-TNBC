rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step4_TNBC_WGCNA')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step3_TNBC_Timer/TCGA_BRCA_Immunity.RData")

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step1_TNBC_ssGSEA/tcga_gsva_BRCA_ssGSEA.Rdata")

load("G:/TCGA_database/TCGA_BRCA/BRCA_exp/TPM/TCGA_BRCA_tpm_ensg_mRNA.Rdata")

#### Step 1 anno BRCA ####

TCGA_BRCA_tpm_ensg_mRNA[1:6,1:6]

library(annoE)

TCGA_BRCA_tpm_symbol_mRNA <- AnnoENSG(data = TCGA_BRCA_tpm_ensg_mRNA) 

rm(TCGA_BRCA_tpm_ensg_mRNA)

TCGA_BRCA_tpm_symbol_mRNA[1:6,1:6]

rownames(TCGA_BRCA_tpm_symbol_mRNA) <- TCGA_BRCA_tpm_symbol_mRNA$SYMBOL

TCGA_BRCA_tpm_symbol_mRNA <- TCGA_BRCA_tpm_symbol_mRNA[,-c(1:5)]

TCGA_BRCA_tpm_symbol_mRNA[1:6,1:6]

group=sapply(strsplit(colnames(TCGA_BRCA_tpm_symbol_mRNA),"\\."),"[",4)

TCGA_BRCA_tpm_symbol_mRNA <- TCGA_BRCA_tpm_symbol_mRNA[,group=='01A']

colnames(TCGA_BRCA_tpm_symbol_mRNA) <- gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*","\\1\\-\\2\\-\\3\\-\\4",colnames(TCGA_BRCA_tpm_symbol_mRNA))

colnames(TCGA_BRCA_tpm_symbol_mRNA) <- substring(colnames(TCGA_BRCA_tpm_symbol_mRNA),1,15)

rm(group)

tmp_name <- intersect(rownames(TCGA_BRCA_CIBERSORT),colnames(TCGA_BRCA_tpm_symbol_mRNA))

TCGA_BRCA_SSGSEA <- tcga_gsva_BRCA

rm(tcga_gsva_BRCA)

TCGA_BRCA_SSGSEA$sample <- rownames(TCGA_BRCA_SSGSEA)

group=sapply(strsplit(TCGA_BRCA_SSGSEA$sample,"\\."),"[",4)

TCGA_BRCA_SSGSEA <- TCGA_BRCA_SSGSEA[group == '01A', ]

TCGA_BRCA_SSGSEA$sample <- gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*","\\1\\-\\2\\-\\3\\-\\4",TCGA_BRCA_SSGSEA$sample)

TCGA_BRCA_SSGSEA$sample <- substring(TCGA_BRCA_SSGSEA$sample,1,15)

TCGA_BRCA_SSGSEA <- aggregate(.~sample,TCGA_BRCA_SSGSEA,median)

rownames(TCGA_BRCA_SSGSEA) <- TCGA_BRCA_SSGSEA$sample

TCGA_BRCA_SSGSEA <- TCGA_BRCA_SSGSEA[,-1]

(tmp_b <- ls())
tmp_b <- tmp_b[-c(1,9,11)]

Immunity_List  <- ls()

Immunity_List <- lapply(tmp_b, function(x){
  
  tmp <- get(x)
  
  tmp <- tmp[rownames(tmp) %in% tmp_name, ]
  
  return(tmp)
  
})

names(Immunity_List) <- tmp_b

TCGA_BRCA_tpm_symbol_mRNA <- subset(TCGA_BRCA_tpm_symbol_mRNA,select = tmp_name)

save(Immunity_List,TCGA_BRCA_tpm_symbol_mRNA,file = 'Result4_rawdata.RData')

#### Step 2 WGCNA ####
rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step4_TNBC_WGCNA')

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step4_TNBC_WGCNA/Result4_rawdata.RData")

## TNBC ##
load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN/Step2_ssGSEA_Score/Step2_TNBC_datTraits.RData")

colnames(TCGA_BRCA_tpm_symbol_mRNA) <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(TCGA_BRCA_tpm_symbol_mRNA))

tmp_samples <- intersect(colnames(TCGA_BRCA_tpm_symbol_mRNA),rownames(TNBC_datTraits))

TCGA_BRCA_tpm_symbol_mRNA <- subset(TCGA_BRCA_tpm_symbol_mRNA,select = tmp_samples)

rm(tmp_samples)

datExpr0 = TCGA_BRCA_tpm_symbol_mRNA[order(apply(TCGA_BRCA_tpm_symbol_mRNA,1,mad), decreasing = T)[1:5000],]
datExpr0 <- data.frame(t(datExpr0))
datExpr0[1:6,1:6]

rm(TCGA_BRCA_tpm_symbol_mRNA)

datExpr <- log2(datExpr0 + 1)
rm(datExpr0)

library(WGCNA)

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

sampleTree = hclust(dist(datExpr), method = "average")

pdf('Plot4_1_Sample clustering to detect outliers.pdf',width = 12,height = 12)

par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
# abline(h = 120, col = "red");

dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 120, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

save(datExpr,file = 'Result4_tidyrawdata_exp.RData')

powers = c(c(1:10), seq(from = 12, to=20, by=1))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("Plot4_2_Threshold.pdf",width = 10, height = 8)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

sft$powerEstimate

softPower = sft$powerEstimate;
adjacency = adjacency(datExpr, power = softPower);

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)

pdf(file = "Plot4_3_geneDendro.pdf", wi = 12, he = 8)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

save(MEs, moduleLabels, moduleColors, geneTree, file = "SingleGene-02-networkConstruction-stepByStep.RData")

Immunity_List <- lapply(Immunity_List, function(x){
  
  rownames(x) <- gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(x))
  
  tmp <- x[rownames(datExpr), ]
  
  return(tmp)
  
})

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#####

i <- 6

pw_tmp <- names(Immunity_List)[i]

dir.create(pw_tmp)

setwd(paste0('./',pw_tmp))

datTraits <- Immunity_List[[i]]

colnames(datTraits)

# colnames(datTraits) <- gsub('_XCELL','',colnames(datTraits))

moduleTraitCor = cor(MEs, datTraits, use = "p");

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf(paste0('Plot4_',names(Immunity_List)[i],'_Module_trait.pdf'),width = 40,height = 40)
# sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");

# textMatrix[c(82)] <- c("0.87\n(2.0e-320)")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.5,
               zlim = c(-1,1),
               main = paste(paste0(names(Immunity_List)[i],' Module-trait relationships')))
dev.off()


# Define variable OS containing the OS column of datTrait
Macrophages = as.data.frame(datTraits$Macrophages);
names(Macrophages) = "Macrophages"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Macrophages, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(Macrophages), sep="");
names(GSPvalue) = paste("p.GS.", names(Macrophages), sep="");

module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;

pdf(paste("Plot4 Module Membership in", module, "module.pdf"),width = 12,height = 8)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Macrophages",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'black')
dev.off()

#### get Anno ####
# BiocManager::install("tidyverse")
library(tidyverse)
substanceBXH <- colnames(datExpr) %>% as.character()
gene_symbol <- substanceBXH

annot <- data.frame(substanceBXH = substanceBXH,gene_symbol = gene_symbol )
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       # LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for Macrophages
modOrder = order(-abs(cor(MEs, Macrophages, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.));
geneInfo = geneInfo0[geneOrder, ]

save(geneInfo,file = 'geneInfo.RData')

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step4_TNBC_WGCNA')











