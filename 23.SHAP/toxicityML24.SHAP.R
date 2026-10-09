library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)
library(catboost)

set.seed(12345)
trainFile   <- "data.train.txt"
testFile    <- "data.test.txt"
geneFile    <- "model.top3_intersect.txt"
methodFile  <- "refer.methodLists.txt"

raw_train <- read.table(trainFile, header=T, sep="\t", check.names=F, row.names=1, na.strings=c("NA",""," "))
raw_test  <- read.table(testFile,  header=T, sep="\t", check.names=F, row.names=1, na.strings=c("NA",""," "))

geneRT <- read.table(geneFile, header=T, sep="\t", check.names=F)
consensus_genes <- as.vector(geneRT[,1])
cat("共识基因数:", length(consensus_genes), "\n")
cat("基因列表:", paste(consensus_genes, collapse=", "), "\n")

train_genes <- intersect(consensus_genes, colnames(raw_train))
cat("训练集中匹配到的基因:", length(train_genes), "\n")
train <- raw_train[, train_genes, drop=FALSE]
train$Type <- as.factor(raw_train$Type)

test_genes <- intersect(consensus_genes, colnames(raw_test))
cat("测试集中匹配到的基因:", length(test_genes), "\n")
test <- raw_test[, test_genes, drop=FALSE]
test$Type <- as.factor(raw_test$Type)

levels(test$Type) <- levels(train$Type)

levels(train$Type) <- c("Control", "Case")
levels(test$Type)  <- c("Control", "Case")

nzv <- nearZeroVar(train[, -ncol(train)], saveMetrics=TRUE)
if (sum(nzv$nzv) > 0) {
  cat("训练集移除", sum(nzv$nzv), "个方差为0的基因\n")
  train <- train[, !nzv$nzv]
  test  <- test[, colnames(train)]
}

if (sum(is.na(train)) > 0) {
  cat("训练集填充", sum(is.na(train)), "个缺失值\n")
  train[, -ncol(train)] <- apply(train[, -ncol(train)], 2,
    function(x) ifelse(is.na(x), median(x, na.rm=T), x))
}

yTestClass <- test$Type
yTest <- ifelse(yTestClass == levels(train$Type)[1], 0, 1)
test_feat <- test[, -ncol(test)]

cat("训练集样本数:", nrow(train), " 基因数:", ncol(train)-1, "\n")
cat("测试集样本数:", nrow(test), " 基因数:", ncol(test)-1, "\n")

control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  savePredictions = TRUE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter = FALSE
)

methodRT <- read.table(methodFile, header=T, sep="\t", check.names=F)
methodRT <- methodRT[!is.na(methodRT$Method) & !is.na(methodRT$Name), ]
cat("方法总数:", nrow(methodRT), "\n")

modelList <- list()
AUC <- c()
ROCcolor <- rainbow(nrow(methodRT))
first_plot <- TRUE

pdf(file="ROC.pdf", width=5.5, height=5)

for(i in 1:nrow(methodRT)) {
  name   <- methodRT[i, "Name"]
  method <- methodRT[i, "Method"]

  cat("正在训练模型:", name, "(", method, ")\n")

  tryCatch({
    if (name == "CatBoost") {
      train_labels <- as.numeric(train$Type) - 1
      train_pool <- catboost.load_pool(
        data = as.matrix(train[, -ncol(train)]),
        label = train_labels
      )
      test_pool <- catboost.load_pool(data = as.matrix(test_feat))

      model <- catboost.train(
        learn_pool = train_pool,
        params = list(
          loss_function = 'Logloss',
          iterations = 500,
          thread_count = 4,
          verbose = 0,
          random_seed = 12345
        )
      )

      pred_prob <- catboost.predict(model, test_pool, prediction_type = 'Probability')
      roc_obj <- roc(yTest, pred_prob)
      auc_value <- as.numeric(roc_obj$auc)
      modelList[[method]] <- auc_value

    } else {
      if (name == "SVM") {
        model <- train(Type ~ ., data = train, method = method,
                       prob.model = TRUE, trControl = control)
      } else {
        model <- train(Type ~ ., data = train, method = method, trControl = control)
      }

      pred <- predict(model, newdata = test_feat, type = "prob")
      pred_col <- ifelse(levels(train$Type)[2] %in% colnames(pred),
                         levels(train$Type)[2], colnames(pred)[2])
      roc_obj <- roc(yTest, as.numeric(pred[, pred_col]))
      auc_value <- as.numeric(roc_obj$auc)
      modelList[[method]] <- auc_value
    }

    AUC <- c(AUC, paste0(name, ': ', sprintf("%.03f", auc_value)))

    if (first_plot) {
      plot(roc_obj, print.auc=F, legacy.axes=T,
           main="ROC Curve - 12 Methods (Test Set: GSE3365+GSE38713)",
           col=ROCcolor[i], lwd=3)
      first_plot <- FALSE
    } else {
      plot(roc_obj, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3, add=T)
    }
  }, error = function(e) {
    cat("  模型", name, "训练失败:", e$message, "\n")
    AUC <- c(AUC, paste0(name, ': 训练失败'))
    modelList[[method]] <- NA
  })

  if (i == nrow(methodRT)) {
    legend('bottomright', AUC, col=ROCcolor[1:length(AUC)], lwd=3, bty='n', cex=0.9)
  }
}

dev.off()
cat("ROC图已保存: ROC.pdf\n")

aucValue <- unlist(modelList)
aucValue <- aucValue[!is.na(aucValue)]
if (length(aucValue) == 0) stop("所有模型训练失败！")

cat("\n=== 测试集AUC排名 ===\n")
auc_sorted <- sort(aucValue, decreasing=TRUE)
for (k in 1:min(5, length(auc_sorted))) {
  cat("  #", k, ":", names(auc_sorted)[k], "=", round(auc_sorted[k], 4), "\n")
}
cat("\n")

cat("===== 使用CatBoost进行SHAP分析 =====\n")

train_labels_all <- as.numeric(train$Type) - 1
train_features_all <- as.matrix(train[, -ncol(train)])

final_pool <- catboost.load_pool(data = train_features_all, label = train_labels_all)
final_model <- catboost.train(
  learn_pool = final_pool,
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

test_pool_final <- catboost.load_pool(data = as.matrix(test_feat))
test_prob <- catboost.predict(final_model, test_pool_final, prediction_type = 'Probability')
test_roc <- roc(yTest, test_prob)
cat("CatBoost测试集AUC:", round(as.numeric(test_roc$auc), 4), "\n")

pred_fun <- function(object, newdata) {
  pool <- catboost.load_pool(data = as.matrix(newdata))
  catboost.predict(object, pool, prediction_type = 'Probability')
}
fit <- permshap(
  object = final_model,
  X = train_features_all,
  pred_fun = pred_fun,
  m = 10
)
shp <- shapviz(fit, X_pred = train_features_all, X = train[, -ncol(train)])

important <- sort(colMeans(abs(shp$S)), decreasing = TRUE)
showVars <- names(important)
cat("SHAP特征重要性排序:\n")
for (k in 1:length(showVars)) {
  cat("  #", k, ":", showVars[k], " (", round(important[k], 4), ")\n")
}

write.table(
  data.frame(gene = names(important), SHAP_mean_abs = round(important, 4)),
  file = "important.genes.txt", sep="\t", row.names=F, quote=F
)

pdf(file="barplot.pdf", width=6, height=5)
print(sv_importance(shp, kind="bar", show_numbers=TRUE) + theme_bw() +
  labs(title = "SHAP Feature Importance (Top 3 Consensus Genes)"))
dev.off()
cat("已保存: barplot.pdf\n")

pdf(file="bee.pdf", width=7, height=5)
print(sv_importance(shp, kind="bee", show_numbers=TRUE) + theme_bw() +
  labs(title = "SHAP Bee Swarm Plot"))
dev.off()
cat("已保存: bee.pdf\n")

pdf(file="dependence.pdf", width=9, height=6)
print(sv_dependence(shp, v = showVars) + theme_bw() +
  labs(title = "SHAP Dependence Plot"))
dev.off()
cat("已保存: dependence.pdf\n")

pdf(file="waterfall.pdf", width=7, height=5)
print(sv_waterfall(shp, row_id = 1) +
  labs(title = "SHAP Waterfall Plot (Sample 1)"))
dev.off()
cat("已保存: waterfall.pdf\n")

pdf(file="force.pdf", width=9, height=5)
print(sv_force(shp, row_id = 1) +
  labs(title = "SHAP Force Plot (Sample 1)"))
dev.off()
cat("已保存: force.pdf\n")

cat("\n========================================\n")
cat("分析完成！\n")
cat("生成的文件:\n")
cat("  ROC.pdf        - 12种方法在测试集上的ROC曲线\n")
cat("  barplot.pdf    - SHAP特征重要性柱状图\n")
cat("  bee.pdf        - SHAP蜂群图\n")
cat("  dependence.pdf - SHAP依赖图\n")
cat("  waterfall.pdf  - SHAP瀑布图\n")
cat("  force.pdf      - SHAP力图\n")
cat("  important.genes.txt - SHAP基因重要性排序\n")
cat("========================================\n")
