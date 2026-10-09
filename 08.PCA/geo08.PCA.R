library(ggplot2)
library(ggpubr)

bioPCA=function(inputFile=null, outFile=null, titleName=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=t(rt)
	Project=gsub("(.*?)\\_.*", "\\1", rownames(data))

	data.pca=prcomp(data)
	pcaPredict=predict(data.pca)
	PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)

	pdf(file=outFile, width=5.5, height=4.25)
	p1=ggscatter(data=PCA, x="PC1", y="PC2", color="Type", shape="Type",
	          ellipse=T, ellipse.type="norm", ellipse.border.remove=F, ellipse.alpha = 0.1,
	          size=2, main=titleName, legend="right")+
	          theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))
	print(p1)
	dev.off()
}

bioPCA(inputFile="merge.preNorm.txt", outFile="PCA.preNorm.pdf", titleName="Before batch correction")
bioPCA(inputFile="merge.normalize.txt", outFile="PCA.normalzie.pdf", titleName="After batch correction")
