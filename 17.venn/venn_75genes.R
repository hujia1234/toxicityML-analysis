library(ggvenn)
library(ggplot2)

base <- "E:/尼古丁UC科研材料"
rd <- function(p, header=FALSE, col=1) {
  x <- read.table(p, header=header, sep="\t", check.names=FALSE,
                  stringsAsFactors=FALSE, quote="", comment.char="", fill=TRUE)
  g <- as.vector(x[, col]); g <- gsub("^ +| +$", "", g)
  unique(g[nzchar(g) & !is.na(g)])
}

DEG       <- rd(file.path(base, "222/11.Disease/diff.txt"), TRUE)
darkgreen <- rd(file.path(base, "222/10.WGCNA/module_darkgreen.txt"))
genecards <- rd(file.path(base, "222/11.Disease/genecards.txt"))
OMIM      <- rd(file.path(base, "222/11.Disease/OMIM.txt"))
TTD       <- rd(file.path(base, "222/11.Disease/TTD.txt"))
Compound  <- rd(file.path(base, "555/17.venn/Compound.txt"))

Disease <- Reduce(union, list(DEG, darkgreen, genecards, OMIM, TTD))

compoundName <- "nicotine"
diseaseName  <- "UC"
geneList <- list()
geneList[[compoundName]] <- Compound
geneList[[diseaseName]]  <- Disease

interGenes <- sort(Reduce(intersect, geneList))
cat("nicotine 靶点数 :", length(Compound), "\n")
cat("UC 疾病集基因数 :", length(Disease), "\n")
cat("交集靶点数      :", length(interGenes), "\n")
cat("交集基因:\n", paste(interGenes, collapse=", "), "\n")

p <- ggvenn(geneList, show_percentage = TRUE,
            stroke_color = "white", stroke_size = 0.5,
            fill_color = c("#E41A1C", "#1E90FF"),
            set_name_color = c("#E41A1C", "#1E90FF"),
            set_name_size = 6, text_size = 4.5)

ggsave("venn_75genes.pdf", p, width = 6, height = 6)
ggsave("venn_75genes.png", p, width = 6, height = 6, dpi = 300)

write.table(interGenes, "interGenes_75.txt", sep="\t", quote=FALSE,
            col.names=FALSE, row.names=FALSE)
networkTab <- rbind(cbind(compoundName, interGenes, "Compound"),
                    cbind(diseaseName,  interGenes, "Disease"))
colnames(networkTab) <- c("Node1", "Node2", "Type")
write.table(networkTab, "net.network_75.txt", sep="\t", quote=FALSE, row.names=FALSE)
nodeTab <- rbind(cbind(compoundName, "Compound"), cbind(diseaseName, "Disease"),
                 cbind(interGenes, "Gene"))
colnames(nodeTab) <- c("Node", "Type")
write.table(nodeTab, "net.node_75.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("已输出: venn_75genes.pdf / venn_75genes.png / interGenes_75.txt\n")
