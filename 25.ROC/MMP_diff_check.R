suppressMessages(library(limma))

checkMMP <- function(expFile, outPrefix) {
  rt = read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
  x = c("MMP1","MMP3","MMP10","TIMP1")
  x = intersect(x, rownames(rt))
  if(length(x)==0){ cat("No target genes in", expFile, "\n"); return(NULL) }
  status = ifelse(grepl("Treat$", colnames(rt)), "UC", "Con")
  cat("\n=== ", expFile, " 样本分组:", sum(status=="UC"), "UC vs", sum(status=="Con"), "Con ===\n")
  design = model.matrix(~factor(status, levels=c("Con","UC")))
  fit = eBayes(lmFit(rt, design))
  allDiff = topTable(fit, coef=2, number=Inf)
  out = data.frame(gene=rownames(allDiff), allDiff[ ,c("logFC","AveExpr","t","P.Value","adj.P.Val","B")])
  write.table(out, file=paste0(outPrefix,".allDiff.txt"), sep="\t", quote=F, row.names=F)
  res = out[out$gene %in% x, , drop=FALSE]
  res$Sig = ifelse(res$adj.P.Val<0.05 & abs(res$logFC)>0.585, ifelse(res$logFC>0.585,"Up","Down"), "Not")
  print(res, digits=4)
  write.table(res, file=paste0(outPrefix,".target.txt"), sep="\t", quote=F, row.names=F)
  pdf(file=paste0(outPrefix,".boxplot.pdf"), width=8, height=6)
  par(mfrow=c(2,2), mar=c(5,4,3,2))
  conVals = list(); ucVals = list()
  for(g in res$gene) {
    if(g %in% rownames(rt)) {
      conVals[[g]] = as.numeric(rt[g, status=="Con"])
      ucVals[[g]]  = as.numeric(rt[g, status=="UC"])
    }
  }
  for(g in res$gene) {
    if(g %in% names(conVals)) {
      boxplot(list(Control=conVals[[g]], UC=ucVals[[g]]),
        col=c("steelblue","tomato"), main=g, ylab="Expression",
        border="grey30", outline=FALSE, boxwex=0.6)
      stripchart(list(conVals[[g]], ucVals[[g]]),
        add=TRUE, vertical=TRUE, method="jitter", pch=16, cex=0.6, col=rgb(0,0,0,0.3))
    }
  }
  dev.off()
  cat("输出: ", outPrefix, ".target.txt / .boxplot.pdf\n")
  return(res)
}

checkMMP("E:/尼古丁UC科研材料/555/21.MLdata/merge.normalize.txt", "merge_MMP")
checkMMP("E:/尼古丁UC科研材料/555/21.MLdata/GSE3365.normalize.txt", "GSE3365_MMP")
checkMMP("E:/尼古丁UC科研材料/555/21.MLdata/GSE38713.normalize.txt", "GSE38713_MMP")
