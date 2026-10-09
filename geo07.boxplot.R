library(reshape2)
library(ggplot2)

bioBoxplot=function(inputFile=null, outFile=null, titleName=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=t(rt)
	Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
	Sample=gsub("(.+)\\_(.+)\\_(.+)", "\\2", rownames(data))
	data=cbind(as.data.frame(data), Sample, Project)

	rt1=melt(data, id.vars=c("Project", "Sample"))
	colnames(rt1)=c("Project","Sample","Gene","Expression")

	pdf(file=outFile, width=10, height=5)
	p=ggplot(rt1, mapping=aes(x=Sample, y=Expression))+
  		geom_boxplot(aes(fill=Project), notch=T, outlier.shape=NA)+
  		ggtitle(titleName)+ theme_bw()+ theme(panel.grid=element_blank())+
  		theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5,size=2), plot.title=element_text(hjust = 0.5))
	print(p)
	dev.off()
}

bioBoxplot(inputFile="merge.preNorm.txt", outFile="boxplot.preNorm.pdf", titleName="Before batch correction")
bioBoxplot(inputFile="merge.normalize.txt", outFile="boxplot.normalzie.pdf", titleName="After batch correction")
