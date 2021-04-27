rm(list = ls())

setwd('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr4_GSE103091_WGCNA')

library(tidyverse)

df <- read.csv('H:\\PHD_database\\Paper_writing\\A_Paper_XVIII_BreastCancer_LXN_ver2\\Step7_OtherGSE\\Chr3_GSE103091_Cibersort\\CIBERSORT.Output_Job18.csv')

load("H:/PHD_database/Paper_writing/A_Paper_XIII_PP_TNBC/Step1_TNBC_ssGSEA/Chr3_download_cli/GSE103091_cli_OS.RData")

Samples_ID <- GSE103091_cli_OS$Samples_ID %>% as.character()

df <- df[df$P.value < 0.05, ]

colnames(df)

data <- df[,-c(24,25,26)]

data$Input.Sample <- substring(data$Input.Sample,1,10)

data <- data[data$Input.Sample %in% Samples_ID, ]

rownames(data) <- data$Input.Sample

data <- data[,-1]

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr1_GSE103091_reanno/GSE103091_exprdf_uniq.RData")

GSE103091_exprdf_uniq <- subset(GSE103091_exprdf_uniq,select = rownames(data))

datExpr0 = GSE103091_exprdf_uniq[order(apply(GSE103091_exprdf_uniq,1,mad), decreasing = T)[1:4995],]

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/Result6_GetFactors.Rdata")

GetFactors[4] <- "HLA-DQB2"
GetFactors[7] <- "HLA-DQA2"

intersect(rownames(datExpr0),GetFactors)

tmp <- GetFactors[!GetFactors %in% rownames(datExpr0)]

tmp_df <- GSE103091_exprdf_uniq[tmp,]

datExpr0 <- rbind(datExpr0,tmp_df)

intersect(rownames(datExpr0),GetFactors)

datExpr0 <- data.frame(t(datExpr0))

datExpr0[1:6,1:6]

datExpr <- log2(datExpr0 + 1)

library(WGCNA)

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

sampleTree = hclust(dist(datExpr), method = "average")

pdf('Chr3_1_Sample clustering to detect outliers.pdf',width = 12,height = 12)

par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 21, col = "red");

dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 21, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

save(datExpr,file = 'Result3_tidyrawdata_exp.RData')

powers = c(c(1:10), seq(from = 11, to=20, by=1))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("Chr3_2_Threshold.pdf",width = 10, height = 8)
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

pdf(file = "Chr3_3_geneDendro.pdf", wi = 12, he = 8)

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

datTraits <- data[rownames(datExpr),]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p");

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf(paste0('Chr3_GSE103091_Module_trait.pdf'),width = 24,height = 24)
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
               main = paste(paste0('GSE103091 Module-trait relationships')))
dev.off()


# Define variable OS containing the OS column of datTrait
Macrophages.M1 = as.data.frame(datTraits$Macrophages.M1);
names(Macrophages.M1) = "Macrophages.M1"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Macrophages.M1, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(Macrophages.M1), sep="");
names(GSPvalue) = paste("p.GS.", names(Macrophages.M1), sep="");

module_all = c('blue','green','turquoise')

for (i in 1:length(module_all)) {
  
  module <- module_all[i]
  
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  pdf(paste("Chr4 Module Membership in", module, "module.pdf"),width = 12,height = 8)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Macrophages.M1",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'black')
  dev.off()
  
}



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
# Order modules by their significance for Macrophages.M1
modOrder = order(-abs(cor(MEs, Macrophages.M1, use = "p")));
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

########

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step7_OtherGSE/Chr4_GSE103091_WGCNA/geneInfo.RData")

load("H:/PHD_database/Paper_writing/A_Paper_XVIII_BreastCancer_LXN_ver2/Step6_TNBC_Consensecluster/Result6_GetFactors.Rdata")

# GetFactors[4] <- "HLA-DQB2"
# GetFactors[7] <- "HLA-DQA2"

geneInfo$geneSymbol <- as.character(geneInfo$geneSymbol)

tmp <- geneInfo[geneInfo$geneSymbol %in% GetFactors, ]













