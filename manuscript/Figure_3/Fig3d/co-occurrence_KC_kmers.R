### INSTALL WGCNA ###
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R") 
biocLite(c("GO.db", "preprocessCore", "impute"))
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg"); 
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6)); 
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep=""); 
biocLite(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
install.packages("/Users/chenghe/Downloads/WGCNA_1.66.tgz", repos = NULL, lib=.Library, type = "source",dependencies = TRUE)
install.packages("robust")

library(WGCNA)
allowWGCNAThreads(nThreads = 4)

### FORMATTING EXP INPUT ###
exp = read.csv("all_KC_p_2.326453e-08_kmer_RC.csv", stringsAsFactors = F,header = FALSE)
head(exp)
datexp = as.data.frame(apply(exp[-1,-1],1,as.numeric))
dim(datexp)
names(datexp) = exp$V1[-1]
rownames(datexp) = as.character(exp[1,-1])
datexp[1:10,1:10]

### FINDING SCALE INDEPENDENCE AND MEAN CONNECTIVITY - 3 ###
powers = c(c(1:10), seq(from = 12, to = 50, by = 2))
sftdatexp = pickSoftThreshold(datexp, powerVector = powers, verbose = 3) # can run correctly on server
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sftdatexp$fitIndices[,1], -sign(sftdatexp$fitIndices[,3])*sftdatexp$fitIndices[,2], 
     xlab = "soft threshold (power)", ylab = "scale free topology model fit, signed R^2", type = "n",
     main = paste("scale independence"))
text(sftdatexp$fitIndices[,1], -sign(sftdatexp$fitIndices[,3])*sftdatexp$fitIndices[,2],
     labels=powers, cex=cex1, col = "red")
plot(sftdatexp$fitIndices[,1], sftdatexp$fitIndices[,5],
     xlab = "soft threshold (power)", ylab = "mean connectivity", type = "n",
     main = paste("Mean Connectivity"))
text(sftdatexp$fitIndices[,1], sftdatexp$fitIndices[,5], labels=powers, cex=cex1, col = "red")

### GENERATING CHOSEN NETWORK FOR CYTOSCAPE ###
net_P12M10 = blockwiseModules(datexp, power = 18,
                             TOMType = "unsigned", minModuleSize = 150,
                             reassignThreshold = 0, mergeCutHeight = 0.1,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = FALSE, verbose = 3, deepSplit = 2)
mergedcolors_P12M10 = labels2colors(net_P12M10$colors)
plotDendroAndColors(net_P12M10$dendrograms[[1]], mergedcolors_P12M10[net_P12M10$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

### SAVE RESULTS ###
moduleLabels = net_P12M10$colors
moduleColors = labels2colors(net_P12M10$colors)
MEs = net_P12M10$MEs
geneTree = net_P12M10$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, file = "DTA_top1k_net_P12M10.RData")
table(moduleColors)
length(table(moduleColors))
### Exporting to Cytoscape ###
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datexp, power = 12)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("DTA_1k_P12M10_Cyt-edges_p0.01.txt", sep=""),
                               nodeFile = paste("DTA_1k_P12M10_Cyt-nodes_p0.01.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.01,
                               nodeNames = colnames(datexp),
                               nodeAttr = moduleColors)
# Export the modules 
out_module = cbind(exp[-1,1],mergedcolors_P12M10)
colnames(out_module) = c("id","module")
write.table(out_module,file="KC_p_2.326453e-08_kmer_RC_P18M150.txt",sep="\t",quote=FALSE,row.names = FALSE)