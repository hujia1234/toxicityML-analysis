base <- "E:/尼古丁UC科研材料"
rd <- function(p, header=FALSE, col=1) {
  x <- read.table(p, header=header, sep="\t", check.names=FALSE,
                  stringsAsFactors=FALSE, quote="", comment.char="", fill=TRUE)
  g <- as.vector(x[, col])
  g <- gsub("^ +| +$", "", g)
  g <- g[nzchar(g) & !is.na(g)]
  unique(g)
}

DEG       <- rd(file.path(base, "222/11.Disease/diff.txt"), header=TRUE, col=1)
turquoise <- rd(file.path(base, "222/10.WGCNA/module_turquoise.txt"))
darkgreen <- rd(file.path(base, "222/10.WGCNA/module_darkgreen.txt"))
genecards <- rd(file.path(base, "222/11.Disease/genecards.txt"))
OMIM      <- rd(file.path(base, "222/11.Disease/OMIM.txt"))
TTD       <- rd(file.path(base, "222/11.Disease/TTD.txt"))

Compound  <- rd(file.path(base, "555/17.venn/Compound.txt"))
Disease   <- rd(file.path(base, "555/17.venn/Disease.txt"))

cat("==== 各集合大小 ====\n")
cat("DEG        :", length(DEG), "\n")
cat("turquoise  :", length(turquoise), "\n")
cat("darkgreen  :", length(darkgreen), "\n")
cat("genecards  :", length(genecards), "\n")
cat("OMIM       :", length(OMIM), "\n")
cat("TTD        :", length(TTD), "\n")
cat("Compound(尼古丁靶点):", length(Compound), "\n")
cat("Disease.txt(现有)   :", length(Disease), "\n\n")

Disease_old <- Reduce(union, list(DEG, turquoise, genecards, OMIM, TTD))
cat("==== 验证 Disease.txt 构成 ====\n")
cat("DEG ∪ turquoise ∪ genecards ∪ OMIM ∪ TTD =", length(Disease_old), "\n")
cat("与 Disease.txt 完全一致? ", setequal(Disease_old, Disease), "\n")
cat("仅在 Disease.txt 而不在并集中:", length(setdiff(Disease, Disease_old)), "\n")
cat("仅在并集而不在 Disease.txt 中:", length(setdiff(Disease_old, Disease)), "\n\n")

Disease_new <- Reduce(union, list(DEG, turquoise, darkgreen, genecards, OMIM, TTD))
cat("==== 方案A: 加上 darkgreen ====\n")
cat("新并集大小:", length(Disease_new), " (增加", length(setdiff(Disease_new, Disease_old)), "个基因)\n")
cat("darkgreen 中不在旧并集的基因数:", length(setdiff(darkgreen, Disease_old)), "\n")

inter_old <- sort(intersect(Compound, Disease_old))
inter_new <- sort(intersect(Compound, Disease_new))
cat("旧交集(尼古丁∩Disease):", length(inter_old), "\n")
cat("新交集              :", length(inter_new), "\n")
cat("新增靶点:", paste(setdiff(inter_new, inter_old), collapse=", "), "\n")
cat("丢失靶点:", paste(setdiff(inter_old, inter_new), collapse=", "), "\n\n")

cat("==== darkgreen ∩ Compound(尼古丁靶点) ====\n")
dg_in_comp <- sort(intersect(darkgreen, Compound))
cat("数量:", length(dg_in_comp), "\n")
cat(paste(dg_in_comp, collapse=", "), "\n\n")

cat("其中【已在】旧 Disease 并集中的:",
    paste(sort(intersect(dg_in_comp, Disease_old)), collapse=", "), "\n")
cat("其中【不在】旧并集、即由 darkgreen 新带入的:",
    paste(sort(setdiff(dg_in_comp, Disease_old)), collapse=", "), "\n\n")

cat("==== turquoise ∩ Compound ====\n")
cat(paste(sort(intersect(turquoise, Compound)), collapse=", "), "\n")

Disease_B <- Reduce(union, list(DEG, darkgreen, genecards, OMIM, TTD))
inter_B <- sort(intersect(Compound, Disease_B))
cat("\n==== 方案B: 用 darkgreen 替换 turquoise ====\n")
cat("交集数量:", length(inter_B), "\n")
cat("相比旧交集新增:", paste(setdiff(inter_B, inter_old), collapse=", "), "\n")
cat("相比旧交集丢失:", paste(setdiff(inter_old, inter_B), collapse=", "), "\n")
