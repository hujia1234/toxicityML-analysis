if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")

library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="merge.normalize.txt"
geneFile="important.genes.txt"
gmtFile="c2.cp.kegg.Hs.symbols.gmt"

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,Type=="Treat",drop=F]

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

for(gene in as.vector(geneRT[,1])){
	dataL=data[,data[gene,]<median(data[gene,]),drop=F]
	dataH=data[,data[gene,]>=median(data[gene,]),drop=F]
	meanL=rowMeans(dataL)
	meanH=rowMeans(dataH)
	meanL[meanL<0.00001]=0.00001
	meanH[meanH<0.00001]=0.00001
	logFC=meanH-meanL
	logFC=sort(logFC, decreasing=T)
	genes=names(logFC)

	gmt=read.gmt(gmtFile)

	kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
	kkTab=as.data.frame(kk)
	kkTab=kkTab[kkTab$pvalue<0.05,]
	write.table(kkTab,file=paste0(gene,".result.txt"),sep="\t",quote=F,row.names = F)

	termNum=5
	kkUp=kkTab[kkTab$NES>0,]
	if(nrow(kkUp)>=termNum){
		showTerm=row.names(kkUp)[1:termNum]
		gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in high ", gene, " group"))
		pdf(file=paste0(gene,".highExp.pdf"), width=6.5, height=5.5)
		print(gseaplot)
		dev.off()
	}

	termNum=5
	kkDown=kkTab[kkTab$NES<0,]
	if(nrow(kkDown)>=termNum){
		showTerm=row.names(kkDown)[1:termNum]
		gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in low ", gene, " group"))
		pdf(file=paste0(gene,".lowExp.pdf"), width=6.5, height=5.5)
		print(gseaplot)
		dev.off()
	}
}
