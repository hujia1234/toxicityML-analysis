library(limma)
library(sva)

files=dir()
files=grep("normalize.txt$", files, value=T)
geneList=list()

for(file in files){
	if(file=="merge.preNorm.txt"){next}
	if(file=="merge.normalize.txt"){next}
    rt=read.table(file, header=T, sep="\t", check.names=F)
    geneNames=as.vector(rt[,1])
    uniqGene=unique(geneNames)
    header=unlist(strsplit(file, "\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}

interGenes=Reduce(intersect, geneList)

allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
    inputFile=files[i]
	if(inputFile=="merge.preNorm.txt"){next}
	if(inputFile=="merge.normalize.txt"){next}
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)
    colnames(rt)=paste0(header[1], "_", colnames(rt))

    if(i==1){
    	allTab=rt[interGenes,]
    }else{
    	allTab=cbind(allTab, rt[interGenes,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}

outTab=rbind(geneNames=colnames(allTab), allTab)
write.table(outTab, file="merge.preNorm.txt", sep="\t", quote=F, col.names=F)

outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.normalize.txt", sep="\t", quote=F, col.names=F)
