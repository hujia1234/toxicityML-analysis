library(dplyr)
library(ggplot2)
library(ggrepel)

logFCfilter=0.585
adj.P.Val.Filter=0.05
diffFile="all.txt"
geneFile="important.genes.txt"

rt=read.table(diffFile, header=T, sep="\t", check.names=F)
row.names(rt)=rt[,1]
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")

rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("green", "grey","red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size=16, hjust=0.5, face = "bold"))

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(rt))
showData=rt[sameGene,]
p1=p+geom_label_repel(data=showData, max.overlaps=20,
                    box.padding=0.2, point.padding=0.2, min.segment.length=0.1,
                    size=3, aes(label=id)) + theme_bw()
pdf(file="vol.pdf", width=5.25, height=4.5)
print(p1)
dev.off()
