
##Followed and adjusted WGCNA tutorial by Steve Horvath- UCLA

####this one
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "DESeq2", "GO.db", "impute", "phyloseq", "preprocessCore"))
#####

BiocManager::install("WGCNA")
BiocManager::install("rlang")
Biocmanager::install("RSQLite")
warnings()
####loading expression data########all 12070
options(stringsAsFactors = FALSE)
CKDdata=read.table("4Callie.txt")
dim(CKDdata)
names(CKDdata)
datExpr0=as.data.frame(t(CKDdata[,-c(1)]))
names(datExpr0) = CKDdata$Gene.Name;
rownames(datExpr0) = names(CKDdata)[-c(1)]
gsg=goodSamplesGenes(datExpr0, verbose=3);gsg$a10K
head(CKDdata)
dim(datExpr0)
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
####CLUSTERING SAMPLES TO FIND OUTLIERS######
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

####loading clinical trait data############NOT USED
traitData=read.csv("CKDTraits.csv",header=FALSE)
dim(traitData)
traitData
dim(datExpr0)
CKDSamples=rownames(datExpr0)
traitRows=match(CKDSamples,traitData$V1)
datTraits=traitData[traitRows,-1]
rownames(traitData)=traitData[traitRows,1]
###expression data in datExpr0 
###clinical traits in traitData
traitData
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry

traitColors = labels2colors(traitData$V2);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitData),
                    main = "Sample dendrogram and trait heatmap")

###automatic network construction and module detection##########don't use


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.75,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")





####save module assignment (from automatic blockwise)
moduleLabels = bwnet$colors
moduleColors = labels2colors(bwnet$colors)
MEs = bwnet$MEs;
geneTree = bwnet$dendrograms[[1]];




# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels);
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)

table(bwLabels) ###modules identified and module size




####block-wise network construction and module detection#####USE

bwnet = blockwiseModules(datExpr0, maxBlockSize = 6035,
                         power = 9, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "CKDTOM-blockwise",
                         verbose = 3)
###bwnet$colors contains the module assignment, and bwnet$MEs contains the module eigengenes of the modules
bwnet$colors
bwnet$MEs
### The hierarchical clustering dendrograms (trees) used for the module identification for each block are
#returned:
bwnet$dendrograms[[1]] ###number of objects=6018
bwnet$dendrograms[[2]] ###number of objects=4863
bwnet$dendrograms[[3]] ###number of objects=1189

###dendograms displayed together with the color assignments:
# open a graphics window
sizeGrWindow(6,6)

# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2

plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 3

plotDendroAndColors(bwnet$dendrograms[[3]], bwModuleColors[bwnet$blockGenes[[3]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 3",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

####not necessary/comparing against self instead of single blockwise
singleBlockMEs = moduleEigengenes(datExpr0, moduleColors)$eigengenes;
blockwiseMEs = moduleEigengenes(datExpr0, bwModuleColors)$eigengenes;

single2blockwise = match(names(singleBlockMEs), names(blockwiseMEs))
signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)

######
sizeGrWindow(12,9)
plotDendroAndColors(geneTree,
                    cbind(moduleColors, bwModuleColors),
                    c("Single block", "2 blocks"),
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# ???Here we show a more flexible way of plotting several trees and colors on one page
sizeGrWindow(12,6)
5
#pdf(file = "Plots/BlockwiseGeneDendrosAndColors.pdf", wi = 12, he = 6);
# Use the layout function for more involved screen sectioning
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))
#layout.show(4);
nBlocks = length(bwnet$dendrograms)
# Plot the dendrogram and the module colors underneath for each block
for (block in 1:nBlocks) 
  plotDendroAndColors(bwnet$dendrograms[[1]], moduleColors[bwnet$blockGenes[[2]]],
                      "Module colors",
                      main = paste("Gene dendrogram and module colors in block", block),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)

####3 Relating modules to external clinical trait#########NOT USED
# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traitData, use = "p",method="spearman");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# DONT USE Will display correlations and their p-values DONT USE BC CATEGORICAL VARIABLE
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))




####Intramodular analysis: identifying genes with high GS and MM########NOT USED
module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for CKD",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# order the categorical variables #####NOT USED
traitData$V2 <- as.numeric(factor(traitData$V2, levels=c("CVD", "normal")))

CorLevelPlot(data = traitData,
             x = c("V2"),
             y = c(MEs),
             col = c("darkblue"),
             cexCorval = 1.5,
             colCorval = "white",
             fontCorval = 2,
             posLab = "bottomleft",
             rotLabX = 45,
             posColKey = "top",
             cexLabColKey = 1.2,
             scale = TRUE,
             main = "World Health Organization",
             colFrame = "white",
             plotRsquared = FALSE)
glm(CKDstatus ~ redModule, family=binomial(link='logit'))
glm(TumourNormalStatus ~ pinkModule, data=data, family=binomial(link='logit'))

CKDstatus = as.data.frame(traitData$V2);

names(CKDstatus) = "Menopausestatus"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr0, CKDstatus, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(CKDstatus), sep="");
names(GSPvalue) = paste("p.GS.", names(CKDstatus), sep="");
module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module
