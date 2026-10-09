suppressMessages({library(limma); library(WGCNA)})
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

expFile <- "E:/尼古丁UC科研材料/222/10.WGCNA/merge.normalize.txt"
outDir  <- "E:/尼古丁UC科研材料/555/10.WGCNA_rerun"
setwd(outDir)

rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp),
               dimnames = list(rownames(exp), colnames(exp)))
data <- avereps(data)
data <- data[apply(data, 1, sd) > 0.5, ]

Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[, order(Type)]
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
conCount <- length(Type[Type == "Control"])
treatCount <- length(Type[Type == "Treat"])
datExpr0 <- t(data)

gsg <- goodSamplesGenes(datExpr0, verbose = 0)
if (!gsg$allOK) datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]

sampleTree <- hclust(dist(datExpr0), method = "average")
pdf("1_sample_cluster.pdf", width = 9, height = 6)
par(cex = 0.6); par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 20000, col = "red")
dev.off()

clust <- cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
keepSamples <- (clust == 1)
datExpr0 <- datExpr0[keepSamples, ]

traitData <- data.frame(Control = c(rep(1, conCount), rep(0, treatCount)),
                        Treat   = c(rep(0, conCount), rep(1, treatCount)))
row.names(traitData) <- colnames(data)
sameSample <- intersect(rownames(datExpr0), rownames(traitData))
datExpr0 <- datExpr0[sameSample, ]
datTraits <- traitData[sameSample, ]

sampleTree2 <- hclust(dist(datExpr0), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)
pdf("2_sample_heatmap.pdf", width = 9, height = 7)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 0)
pdf("3_scale_independence.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.8, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

fitOut <- data.frame(Power = sft$fitIndices[, 1],
                     SFT.R.sq = signif(sft$fitIndices[, 2], 4),
                     slope = signif(sft$fitIndices[, 3], 4),
                     mean.k. = signif(sft$fitIndices[, 5], 4),
                     signedR2 = signif(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 4))
write.table(fitOut, "soft_threshold_fitIndices.txt", sep = "\t", row.names = FALSE, quote = FALSE)

softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 6

adjMat <- adjacency(datExpr0, power = softPower)
TOM <- TOMsimilarity(adjMat)
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")
pdf("4_gene_clustering.pdf", width = 8, height = 6)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,
                             pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
pdf("5_Dynamic_Tree.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MEList <- moduleEigengenes(datExpr0, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
MEDissThres <- 0.25
pdf("6_Clustering_module.pdf", width = 7, height = 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = MEDissThres, col = "red")
dev.off()

mergeRes <- mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 0)
mergedColors <- mergeRes$colors
mergedMEs <- mergeRes$newMEs
pdf("7_merged_dynamic.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, mergedColors, "Merged dynamic",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors <- mergedColors
MEs <- mergedMEs

nGenes <- ncol(datExpr0)
nSamples <- nrow(datExpr0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
pdf("8_Module_trait.pdf", width = 5.5, height = 5.5)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = FALSE, cex.text = 0.75, zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

mtOut <- data.frame(Module = rownames(moduleTraitCor),
                    Control_r = signif(moduleTraitCor[, 1], 4),
                    Control_P = signif(moduleTraitPvalue[, 1], 4),
                    Treat_r = signif(moduleTraitCor[, 2], 4),
                    Treat_P = signif(moduleTraitPvalue[, 2], 4))
write.table(mtOut, "module_trait_correlation.txt", sep = "\t", row.names = FALSE, quote = FALSE)

modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")
traitNames <- names(datTraits)
geneTraitSignificance <- as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", traitNames, sep = "")
names(GSPvalue) <- paste("p.GS.", traitNames, sep = "")

y <- datTraits[, 1]
GS1 <- as.numeric(cor(y, datExpr0, use = "p"))
GeneSignificance <- abs(GS1)
ModuleSignificance <- tapply(GeneSignificance, mergedColors, mean, na.rm = TRUE)
pdf("9_GeneSignificance.pdf", width = 12.5, height = 7.5)
plotModuleSignificance(GeneSignificance, mergedColors)
dev.off()

trait <- "Treat"
traitColumn <- match(trait, traitNames)
for (module in modNames) {
  column <- match(module, modNames)
  moduleGenes <- moduleColors == module
  if (nrow(geneModuleMembership[moduleGenes, ]) > 1) {
    outPdf <- paste("10_", trait, "_", module, ".pdf", sep = "")
    pdf(file = outPdf, width = 7, height = 7)
    par(mfrow = c(1, 1))
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, traitColumn]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for ", trait),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    dev.off()
  }
}

probes <- colnames(datExpr0)
geneInfo0 <- data.frame(probes = probes, moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance)) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneTraitSignificance[, Tra], GSPvalue[, Tra])
  names(geneInfo0) <- c(oldNames, names(geneTraitSignificance)[Tra], names(GSPvalue)[Tra])
}
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[, mod], MMPvalue[, mod])
  names(geneInfo0) <- c(oldNames, names(geneModuleMembership)[mod], names(MMPvalue)[mod])
}
geneInfo <- geneInfo0[order(geneInfo0$moduleColor), ]
write.table(geneInfo, "GS_MM.xls", sep = "\t", row.names = FALSE)

cnt <- table(moduleColors)
for (mod in 1:nrow(cnt)) {
  modules <- names(cnt)[mod]
  write.table(probes[moduleColors == modules], file = paste0("module_", modules, ".txt"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

signedR2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
topIdx <- which.max(abs(moduleTraitCor[, 2]))
turIdx <- which(rownames(moduleTraitCor) == "MEturquoise")
mm <- abs(geneModuleMembership[moduleColors == "turquoise", "MMturquoise"])
gs <- abs(geneTraitSignificance[moduleColors == "turquoise", "GS.Treat"])
ct <- cor.test(mm, gs)

sink("WGCNA_summary.txt")
cat("===================== WGCNA 重跑汇总 =====================\n")
cat("运行时间      :", format(Sys.time()), "\n")
cat("输入文件      :", expFile, "\n")
cat("输出目录      :", outDir, "\n\n")
cat("样本数        :", nSamples, "  (Control", conCount, "/ Treat", treatCount, ")\n")
cat("基因数(过滤后):", nGenes, "\n\n")
cat("softPower     :", softPower, "  [pickSoftThreshold powerEstimate =",
    ifelse(is.na(sft$powerEstimate), "NA", sft$powerEstimate), "]\n")
cat("R2 @ softPower:", signif(signedR2[sft$fitIndices[, 1] == softPower], 4), "\n")
cat("mean.k @ softPower:", signif(sft$fitIndices[sft$fitIndices[, 1] == softPower, 5], 4), "\n")
cat("R2 最大值(1:20):", signif(max(signedR2), 4), "  at power",
    sft$fitIndices[which.max(signedR2), 1], "\n")
cat("首个 signedR2>=0.8 的 power:", min(sft$fitIndices[, 1][signedR2 >= 0.8], na.rm = TRUE), "\n\n")
cat("minModuleSize :", minModuleSize, "\n")
cat("deepSplit     : 2\n")
cat("MEDissThres   :", MEDissThres, "  (1-cor >", 1 - MEDissThres, " 合并)\n")
cat("动态模块数    :", length(unique(dynamicMods)), "\n")
cat("合并后模块数  :", length(cnt), "  (非灰", sum(names(cnt) != "grey"), ")\n\n")
cat("各模块基因数:\n"); print(cnt)
cat("\n与 Treat 相关性最强的模块:", rownames(moduleTraitCor)[topIdx],
    "  |r| =", signif(abs(moduleTraitCor[topIdx, 2]), 4), "\n")
cat("turquoise 与 Treat : r =", signif(moduleTraitCor[turIdx, 2], 4),
    "  P =", signif(moduleTraitPvalue[turIdx, 2], 4), "\n")
cat("turquoise 基因数   :", as.integer(cnt["turquoise"]), "\n")
cat("turquoise MM-vs-GS : n =", length(mm), "  cor =", signif(ct$estimate, 4),
    "  P =", signif(ct$p.value, 4), "\n\n")
cat("模块-性状相关性（全部模块）:\n"); print(mtOut, row.names = FALSE)
sink()

cat("DONE\n")
