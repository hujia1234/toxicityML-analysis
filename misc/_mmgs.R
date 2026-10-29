base <- Sys.getenv("TOXICITY_DATA_DIR", unset = ".")
d <- read.table(file.path(base, "10.WGCNA/GS_MM.xls"),
                header=TRUE, sep="\t", check.names=FALSE, comment.char="")
cat("列名:", paste(names(d)[1:8], collapse=" | "), "\n")
cat("总基因数:", nrow(d), "\n")
cat("darkgreen 基因数:", sum(d$moduleColor=="darkgreen"), "\n")
cat("turquoise 基因数:", sum(d$moduleColor=="turquoise"), "\n\n")

for (m in c("darkgreen","turquoise")) {
  sub <- d[d$moduleColor==m, ]
  x <- abs(as.numeric(sub[[paste0("MM", m)]]))
  y <- abs(as.numeric(sub[["GS.Treat"]]))
  ct <- cor.test(x, y)
  cat(sprintf("%s  MM-vs-GS(Treat): n=%d  cor=%.4f  P=%.4g\n",
              m, nrow(sub), ct$estimate, ct$p.value))
}
