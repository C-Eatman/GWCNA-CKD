####followed and adjusted WGCNA tutorial by Steve Horvath- UCLA


####loading expression data NORMAL
options(stringsAsFactors = FALSE)
Ref2data=read.csv("CVDReferencedata.csv")
dim(Ref2data)
names(Ref2data)
datExpr3=as.data.frame(t(Ref2data[,-c(1)]))
names(datExpr3) = Ref2data$Gene.Name;
rownames(datExpr3) = names(Ref2data)[-c(1)]
gsg3=goodSamplesGenes(datExpr3, verbose=3);gsg$a10K
head(Ref3data)
dim(datExpr3)

####CLUSTERING SAMPLES TO FIND OUTLIERS######
sampleTree3 = hclust(dist(datExpr3), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
abline(h=14,col="red")
plot(sampleTree3, main = "Sample clustering to detect outliers-Normal", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

names(Ref3data)
Ref3data<-Ref2data[,-4]

###automatic network construction and module detection


# Choose a set of soft-thresholding powers
powers3 = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr3, powerVector = powers3, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers3,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers3, cex=cex1,col="red")

net3 = blockwiseModules(datExpr3, power = 14,
                        TOMType = "unsigned", minModuleSize = 15,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "RefTOM", 
                        verbose = 3)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors3 = labels2colors(net3$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net3$dendrograms[[1]], mergedColors3[net3$blockGenes[[1]]],
                    "Module colors",
                    main=("Dendrogram of modules in normal data"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels3 = net3$colors
moduleColors3 = labels2colors(net3$colors)
MEs3 = net3$MEs;
geneTree3 = net3$dendrograms[[1]];

consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]

net3$colors #module assignments
net3$MEs #module eigangenes
####CKD vs Normal

# Rename variables to avoid conflicts
CKDLabels = moduleLabels1;
CKDColors = moduleColors1;
CKDTree = geneTree1;
CKDMEs = orderMEs(MEs1, greyName = "ME0")
####The consensus network analysis results are represented by the variables consMEs, moduleLabels, moduleColors, and
#consTree. We are now ready to relate the female modules to the consensus modules. We calculate the overlaps
#of each pair of female-consensus modules, and use the Fisher’s exact test (also known as hypergeometric test) to
#assign a p-value to each of the pairwise overlaps.

# Isolate the module labels in the order they appear in ordered module eigengenes
CKDModuleLabels = substring(names(CKDMEs), 3)
NormalModuleLabels = substring(names(MEs3), 3)
# Convert the numeric module labels to color labels
CKDModules = labels2colors(as.numeric(CKDModuleLabels))
NormalModules = labels2colors(as.numeric(NormalModuleLabels))
# Numbers of female and consensus modules
nCKDMods = length(CKDModules)
nNormalMods = length(NormalModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nCKDMods, ncol = nNormalMods);
CountTbl = matrix(0, nrow = nCKDMods, ncol = nNormalMods);
# Execute all pairwaise comparisons
for (fmod in 1:nCKDMods)
  for (nmod in 1:nNormalMods)
  {
    CKDMembers = (CKDColors == CKDModules[fmod]);
    NormalMembers = (moduleColors3 == NormalModules[nmod]);
    pTable[fmod, nmod] = -log10(fisher.test(CKDMembers, NormalMembers, alternative = "greater")$p.value);
    CountTbl[fmod, nmod] = sum(CKDColors == CKDModules[fmod] & moduleColors3 ==
                                 NormalModules[nmod])
  }



###exporting results of the network analysis
probes = names(datExpr2[[1]]$data)
probes2Refdata = match(probes, Testdata$Gene.Name)
#######ONLY USE THIS 
consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels3, excludeGrey = TRUE)
consMEs.unord
kME = list();
for (set in 1:nSets)
{
 
  kME[[set]] = corAndPvalue(multiExpr()[[set]]$data, consMEs.unord[[set]]$data);
}
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
  c("kME.set1.", "kME.set2.", "p.kME.set1.", "p.kME.set2.", "Z.kME.meta.", "p.kME.meta"),
  rep(MEnames, rep(6, nMEs)))
