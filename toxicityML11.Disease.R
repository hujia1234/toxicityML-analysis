library(ggvenn)
diffFile="diff.txt"
moduleFile="module_red.txt"
setwd("C:\\Users\\Administrator\\Desktop\\toxicityML\\11.Disease")
geneList=list()

rt=read.table("diff.txt", header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])
geneNames=gsub("^ | $","",geneNames)
uniqGene=unique(geneNames)
geneList[["DEG"]]=uniqGene

rt=read.table(moduleFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])
geneNames=gsub("^ | $","",geneNames)
uniqGene=unique(geneNames)
geneList[["WGCNA"]]=uniqGene

pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList,show_percentage = T,
  stroke_color = "white", stroke_size = 0.5,
  fill_color = c("#E41A1C","#1E90FF"),
  set_name_color =c("#E41A1C","#1E90FF"),
  set_name_size=6, text_size=4.5)
dev.off()

unionGenes=Reduce(union, geneList)
write.table(file="Disease.txt", unionGenes, sep="\t", quote=F, col.names=F, row.names=F)
