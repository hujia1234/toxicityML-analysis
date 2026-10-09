scoreFile1="score1.txt"
scoreFile2="score2.txt"

score1=read.table(scoreFile1, header=T, sep="\t",check.names =FALSE,row.names=1)
nCol=ncol(score1)

mat=data.frame()
for(i in colnames(score1)){
	print(paste0("filter1,",i,": ",median(score1[,i])))
	value=ifelse(score1[,i]>median(score1[,i]),1,0)
	mat=rbind(mat,value)
}
mat=t(mat)
colnames(mat)=colnames(score1)
row.names(mat)=row.names(score1)
geneName=row.names(mat[rowSums(mat)==nCol,])
score2=score1[geneName,]
score2out=cbind(name=row.names(score2),score2)
write.table(score2out, file=scoreFile2, sep='\t', quote=F, row.names=F)
write.table(geneName, file="score2.gene.txt", sep='\t', quote=F, row.names=F,col.names=F)
