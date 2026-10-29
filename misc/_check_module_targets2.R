base <- Sys.getenv("TOXICITY_DATA_DIR", unset = ".")
rd <- function(p, header=FALSE, col=1) {
  x <- read.table(p, header=header, sep="\t", check.names=FALSE,
                  stringsAsFactors=FALSE, quote="", comment.char="", fill=TRUE)
  g <- as.vector(x[, col]); g <- gsub("^ +| +$", "", g)
  unique(g[nzchar(g) & !is.na(g)])
}
DEG       <- rd(file.path(base, "11.Disease/diff.txt"), TRUE)
turquoise <- rd(file.path(base, "10.WGCNA/module_turquoise.txt"))
darkgreen <- rd(file.path(base, "10.WGCNA/module_darkgreen.txt"))
genecards <- rd(file.path(base, "11.Disease/genecards.txt"))
OMIM      <- rd(file.path(base, "11.Disease/OMIM.txt"))
TTD       <- rd(file.path(base, "11.Disease/TTD.txt"))
Compound  <- rd(file.path(base, "17.venn/Compound.txt"))
Disease   <- rd(file.path(base, "17.venn/Disease.txt"))
inter     <- sort(intersect(Compound, Disease))

cat("当前 73 靶点:\n"); cat(paste(inter, collapse=", "), "\n\n")

Disease_noWGCNA <- Reduce(union, list(DEG, genecards, OMIM, TTD))
cat("若不使用任何 WGCNA 模块, 交集 =", length(intersect(Compound, Disease_noWGCNA)), "\n")
cat("turquoise 独立贡献(去掉后丢失的):",
    paste(setdiff(inter, intersect(Compound, Disease_noWGCNA)), collapse=", "), "\n\n")

cat("EPHX2 在 genecards?", "EPHX2" %in% genecards,
    " 在 DEG?", "EPHX2" %in% DEG,
    " 在 turquoise?", "EPHX2" %in% turquoise, "\n")

Disease_noDG <- Reduce(union, list(DEG, turquoise, genecards, OMIM, TTD))
cat("\ndarkgreen 独立带入的新靶点:",
    paste(sort(setdiff(intersect(Compound, Reduce(union, list(Disease_noDG, darkgreen))),
                       intersect(Compound, Disease_noDG))), collapse=", "), "\n")

for (g in c("PLAT","PTAFR")) {
  cat(sprintf("%s -> DEG:%s genecards:%s turquoise:%s darkgreen:%s OMIM:%s TTD:%s\n",
              g, g %in% DEG, g %in% genecards, g %in% turquoise, g %in% darkgreen,
              g %in% OMIM, g %in% TTD))
}
