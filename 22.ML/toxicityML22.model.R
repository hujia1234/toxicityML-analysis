library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
library(ncvreg)
library(Boruta)
library(mRMRe)
library(lightgbm)
library(catboost)

source("refer.ML.R")

Train_data <- read.table("data.train.txt", header = T, sep = "\t", check.names=F, row.names=1, stringsAsFactors=F)
Train_expr=Train_data[,1:(ncol(Train_data)-1),drop=F]
Train_class=Train_data[,ncol(Train_data),drop=F]

Test_data <- read.table("data.test.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors = F)
Test_expr=Test_data[,1:(ncol(Test_data)-1),drop=F]
Test_class=Test_data[,ncol(Test_data),drop=F]
Test_class$Cohort=gsub("(.*)\\_(.*)\\_(.*)", "\\1", row.names(Test_class))
Test_class=Test_class[,c("Cohort", "Type")]

comgene <- intersect(colnames(Train_expr), colnames(Test_expr))
Train_expr <- as.matrix(Train_expr[,comgene])
Test_expr <- as.matrix(Test_expr[,comgene])
Train_set = scaleData(data=Train_expr, cohort = NULL, centerFlags=T, scaleFlags=T)
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)

clean_numeric_data <- function(mat) {
  mat[is.infinite(mat)] <- NA
  mat <- mat[!apply(is.na(mat), 1, all), ]
  mat <- mat[, !apply(is.na(mat), 2, all)]
  mat <- apply(mat, 2, function(col) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    return(col)
  })
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  return(mat)
}

Train_set <- clean_numeric_data(Train_set)
Test_set <- clean_numeric_data(Test_set)

comgene_clean <- intersect(colnames(Train_set), colnames(Test_set))
Train_set <- Train_set[, comgene_clean]
Test_set <- Test_set[, comgene_clean]

cat("训练集异常值数量：", sum(is.na(Train_set) | is.infinite(Train_set)), "\n")
cat("测试集异常值数量：", sum(is.na(Test_set) | is.infinite(Test_set)), "\n")

methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)
methods=methodRT$Model
methods <- methods[!grepl("^#", methods)]
methods <- gsub("-| ", "", methods)
classVar = "Type"
min.selected.var = 5
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))

preTrain.var <- list()
set.seed(seed = 123)
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method,
                                 Train_set = Train_set,
                                 Train_label = Train_class,
                                 mode = "Variable",
                                 classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(Train_set)

model <- list()
set.seed(seed = 123)
Train_set_bk = Train_set
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method)
  Variable = preTrain.var[[method[1]]]
  if (is.null(Variable) || length(Variable) == 0) {
    cat("  Skip: no variables from stage 1\n")
    next
  }
  Train_set = Train_set_bk[, Variable, drop = FALSE]
  Train_label = Train_class
  model[[method_name]] <- RunML(method = method[2],
                                Train_set = Train_set,
                                Train_label = Train_label,
                                mode = "Model",
                                classVar = classVar)

  if(!is.null(model[[method_name]]) && length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
Train_set = Train_set_bk; rm(Train_set_bk)
saveRDS(model, "model.MLmodel.rds")

FinalModel <- c("panML", "multiLogistic")[2]
if (FinalModel == "multiLogistic" && length(model) > 0){
  logisticmodel <- lapply(model, function(fit){
    vars <- ExtractVar(fit)
    if (length(vars) > 0) {
      tmp <- glm(formula = Train_class[[classVar]] ~ .,
                 family = "binomial",
                 data = as.data.frame(Train_set[, vars, drop = FALSE]))
      tmp$subFeature <- vars
      return(tmp)
    } else {
      return(NULL)
    }
  })
  logisticmodel <- logisticmodel[!sapply(logisticmodel, is.null)]
  saveRDS(logisticmodel, "model.logisticmodel.rds")
}

if (file.exists("model.MLmodel.rds")) {
  model <- readRDS("model.MLmodel.rds")
  methodsValid <- names(model)

  if (length(methodsValid) > 0) {
    RS_list <- list()
    for (method in methodsValid){
      RS_list[[method]] <- CalPredictScore(fit = model[[method]], new_data = rbind.data.frame(Train_set,Test_set))
    }
    riskTab=as.data.frame(t(do.call(rbind, RS_list)))
    riskTab=cbind(id=row.names(riskTab), riskTab)
    write.table(riskTab, "model.riskMatrix.txt", sep="\t", row.names=F, quote=F)

    Class_list <- list()
    for (method in methodsValid){
      Class_list[[method]] <- PredictClass(fit = model[[method]], new_data = rbind.data.frame(Train_set,Test_set))
    }
    Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
    classTab=cbind(id=row.names(Class_mat), Class_mat)
    write.table(classTab, "model.classMatrix.txt", sep="\t", row.names=F, quote=F)

    fea_list <- list()
    for (method in methodsValid) {
      fea_list[[method]] <- ExtractVar(model[[method]])
    }
    fea_df <- lapply(model, function(fit){
      data.frame(ExtractVar(fit))
    })
    fea_df <- do.call(rbind, fea_df)
    fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
    colnames(fea_df)[1] <- "features"
    write.table(fea_df, file="model.genes.txt", sep = "\t", row.names = F, col.names = T, quote = F)

    AUC_list <- list()
    for (method in methodsValid){
      AUC_list[[method]] <- RunEval(fit = model[[method]],
                                    Test_set = Test_set,
                                    Test_label = Test_class,
                                    Train_set = Train_set,
                                    Train_label = Train_class,
                                    Train_name = "Train",
                                    cohortVar = "Cohort",
                                    classVar = classVar)
    }
    AUC_mat <- do.call(rbind, AUC_list)
    aucTab=cbind(Method=row.names(AUC_mat), AUC_mat)
    write.table(aucTab, "model.AUCmatrix.txt", sep="\t", row.names=F, quote=F)

    AUC_mat <- read.table("model.AUCmatrix.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)

    avg_AUC <- apply(AUC_mat, 1, mean)
    avg_AUC <- sort(avg_AUC, decreasing = T)
    AUC_mat <- AUC_mat[names(avg_AUC),]
    fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
    avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

    CohortCol <- c("red", "blue")
    if(ncol(AUC_mat)>2){
      CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")}
    names(CohortCol) <- colnames(AUC_mat)

    cellwidth = 1; cellheight = 0.5
    hm <- SimpleHeatmap(Cindex_mat = AUC_mat,
                        avg_Cindex = avg_AUC,
                        CohortCol = CohortCol,
                        barCol = "steelblue",
                        cellwidth = cellwidth, cellheight = cellheight,
                        cluster_columns = F, cluster_rows = F)

    pdf(file="model.AUCheatmap.pdf", width=cellwidth * ncol(AUC_mat) + 6, height=cellheight * nrow(AUC_mat) * 0.45)
    draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
    dev.off()
  } else {
    cat("警告：没有有效的机器学习模型，跳过后续分析\n")
  }
} else {
  cat("警告：未找到模型文件model.MLmodel.rds\n")
}
