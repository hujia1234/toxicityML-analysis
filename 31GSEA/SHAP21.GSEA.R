library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
library(GSEABase)

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

Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,Type=="Treat",drop=F]

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
geneIdsList=geneIds(geneSets)
gmt2gene=data.frame(term=rep(names(geneIdsList), times=sapply(geneIdsList, length)),
                    gene=unlist(geneIdsList))

for(gene in as.vector(geneRT[,1])){
	lowName=colnames(data)[data[gene,]<median(data[gene,])]
	highName=colnames(data)[data[gene,]>=median(data[gene,])]

	group=c(rep(1,length(lowName)), rep(2,length(highName)))
	design=model.matrix(~group)
	fit=lmFit(data[,c(lowName,highName)], design)
	fit=eBayes(fit)
	allDiff=topTable(fit, coef=2, number=Inf)
	allDiff=allDiff[order(allDiff$logFC, decreasing=TRUE),]

	genelist=allDiff$logFC
	names(genelist)=rownames(allDiff)
	genelist=genelist[!duplicated(names(genelist))]
	genelist=sort(genelist, decreasing=TRUE)

	kk=GSEA(genelist, TERM2GENE=gmt2gene, pvalueCutoff=1, minGSSize=5, maxGSSize=500, seed=123456)
	kkTab=as.data.frame(kk)
	kkTab=kkTab[kkTab$pvalue<0.05,]
	write.table(kkTab, file=paste0(gene,".gsea.result.txt"), sep="\t", quote=F, row.names=F)

	termNum=5
	kkUp=kkTab[kkTab$NES>0,]
	if(nrow(kkUp)>=1){
		showNum=min(termNum, nrow(kkUp))
		showTerm=row.names(kkUp)[1:showNum]
		gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in high ", gene, " group"))
		pdf(file=paste0(gene,".highExp.pdf"), width=6.5, height=5.5)
		print(gseaplot)
		dev.off()
	}

	kkDown=kkTab[kkTab$NES<0,]
	if(nrow(kkDown)>=1){
		showNum=min(termNum, nrow(kkDown))
		showTerm=row.names(kkDown)[1:showNum]
		gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in low ", gene, " group"))
		pdf(file=paste0(gene,".lowExp.pdf"), width=6.5, height=5.5)
		print(gseaplot)
		dev.off()
	}
}
