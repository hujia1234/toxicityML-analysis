library(catboost)
library(shapviz)
library(ggplot2)
library(kernelshap)

raw_train <- read.table("data.train.txt", header=T, sep="\t", check.names=F, row.names=1)
geneRT <- read.table("model.top3_intersect.txt", header=T, sep="\t", check.names=F)
consensus_genes <- as.vector(geneRT[,1])

train <- raw_train[, intersect(consensus_genes, colnames(raw_train))]
train$Type <- as.factor(raw_train$Type)
levels(train$Type) <- c("Control", "Case")

train_labels <- as.numeric(train$Type) - 1
train_features <- as.matrix(train[, -ncol(train)])
pool <- catboost.load_pool(data = train_features, label = train_labels)
model <- catboost.train(
  learn_pool = pool,
  params = list(
    loss_function = 'Logloss',
    iterations = 1000,
    thread_count = 4,
    verbose = 0,
    random_seed = 12345,
    depth = 6,
    learning_rate = 0.1
  )
)

pred_fun <- function(object, newdata) {
  p <- catboost.load_pool(data = as.matrix(newdata))
  catboost.predict(object, p, prediction_type = 'Probability')
}
fit <- permshap(model, X = train_features, pred_fun = pred_fun, m = 10)
shp <- shapviz(fit, X_pred = train_features, X = train[, -ncol(train)])

p1 <- sv_dependence(shp, v = "TIMP1", color_factor = 3) +
  theme_bw(base_size = 14) +
  labs(
    title = "TIMP1 SHAP Dependence Plot",
    subtitle = "Color: SOD1 expression (3rd most important gene)",
    x = "TIMP1 expression value",
    y = "SHAP value (impact on prediction)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    legend.position = "bottom",
    legend.title = element_text(size = 10)
  )

p1 <- p1 + geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5)

pdf("dependence_TIMP1.pdf", width = 7, height = 5.5)
print(p1)
dev.off()
cat("已保存: dependence_TIMP1.pdf\n")

p2 <- sv_dependence(shp, v = "TIMP1", color_var = "SOD1") +
  theme_bw(base_size = 14) +
  labs(
    title = "TIMP1 SHAP Dependence (colored by SOD1)",
    x = "TIMP1 expression value",
    y = "SHAP value (impact on prediction)",
    color = "SOD1 expression"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 10)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5)

pdf("dependence_TIMP1_SOD1.pdf", width = 7, height = 5.5)
print(p2)
dev.off()
cat("已保存: dependence_TIMP1_SOD1.pdf\n")

cat("\n=== TIMP1 数据摘要 ===\n")
cat("TIMP1表达量范围:", round(range(train$TIMP1), 3), "\n")
cat("TIMP1 SHAP值范围:", round(range(shp$S[, "TIMP1"]), 4), "\n")
cat("\nTIMP1相关性:\n")
cat("  表达量 vs SHAP值: ", round(cor(train$TIMP1, shp$S[, "TIMP1"]), 4), "\n")

if (cor(train$TIMP1, shp$S[, "TIMP1"]) > 0) {
  cat("  结论: TIMP1表达越高 → 毒性预测概率越高（正相关）\n")
} else {
  cat("  结论: TIMP1表达越高 → 毒性预测概率越低（负相关）\n")
}

cat("\n完成! 生成文件: dependence_TIMP1.pdf, dependence_TIMP1_SOD1.pdf\n")
