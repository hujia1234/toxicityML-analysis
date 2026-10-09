library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

pvalueFilter=0.05
p.adjustFilter=0.05

rt=read.table("interGenes.txt", header=F, sep="\t", check.names=F)

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]

kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<p.adjustFilter),]
KEGG=na.omit(KEGG)
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

showNum = ifelse(nrow(KEGG)>30, 30, nrow(KEGG))
showKEGG = KEGG[1:showNum,]
showKEGG$pathway = factor(showKEGG$Description, levels = rev(showKEGG$Description))

bar <- ggplot(data = showKEGG, aes(x = Count, y = pathway, fill = p.adjust)) +
	geom_bar(stat = "identity", width = 0.75) +
	scale_fill_distiller(palette = "Spectral", direction = 1) +
	labs(x = "Gene count", y = "", title = "Enriched KEGG Pathways") +
	theme_bw()+ theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
		 axis.title = element_text(size = 11), axis.text = element_text(size = 10),
          legend.title = element_text(size = 11), legend.text = element_text(size = 10))
pdf(file="barplot.pdf", width=6.5, height=5.5)
print(bar)
dev.off()

showKEGG$GeneRatio=sapply(showKEGG$GeneRatio, function(x) eval(parse(text = x)))
showKEGG=showKEGG[order(showKEGG$GeneRatio),]
showKEGG$pathway = factor(showKEGG$pathway, levels = showKEGG$pathway)
bubble <- ggplot(data = showKEGG, aes(x = GeneRatio, y = pathway)) +
	geom_point(aes(size = Count, color = p.adjust)) +
	scale_color_distiller(palette = "Spectral", direction = 1) +
	labs(x = "Gene ratio", y = "", title = "Enriched KEGG Pathways", size = "Count") +
	theme_bw()+ theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
		axis.title = element_text(size = 11), axis.text = element_text(size = 10),
          legend.title = element_text(size = 11), legend.text = element_text(size = 10))
pdf(file="bubble.pdf", width=6.5, height=5.5)
print(bubble)
dev.off()
