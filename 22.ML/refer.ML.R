RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
  method = gsub(" ", "", method)
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method)

  if (ncol(Train_set) < 1) {
    warning(paste("Skip", method_name, ": No features available (0 columns)"))
    return(NULL)
  }
  if (method_name %in% c("Enet", "Lasso", "Ridge", "SCAD", "MCP") && ncol(Train_set) < 2) {
    warning(paste("Skip", method_name, ": Need at least 2 features (only", ncol(Train_set), "available)"))
    return(NULL)
  }

  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    NULL
  )
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_set), " Variables")

  args = list("Train_set" = Train_set,
              "Train_label" = Train_label,
              "mode" = mode,
              "classVar" = classVar)
  args = c(args, method_param)

  obj <- do.call(what = paste0("Run", method_name),
                 args = args)

  if(mode == "Variable"){
    message(length(obj), " Variables retained;\n")
  }else{message("\n")}
  return(obj)
}

RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  if (ncol(Train_set) < 2) {
    warning(paste("glmnet requires at least 2 columns (only", ncol(Train_set), "available)"))
    return(NULL)
  }
  if (sum(is.na(Train_set)) > 0) {
    Train_set[is.na(Train_set)] <- 0
  }

  cv.fit = cv.glmnet(x = Train_set,
                     y = Train_label[[classVar]],
                     family = "binomial", alpha = alpha, nfolds = 10)
  fit = glmnet(x = Train_set,
               y = Train_label[[classVar]],
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLasso <- function(Train_set, Train_label, mode, classVar){
  if (ncol(Train_set) < 2) {
    warning("Lasso needs at least 2 features")
    return(NULL)
  }
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
}

RunRidge <- function(Train_set, Train_label, mode, classVar){
  if (ncol(Train_set) < 2) {
    warning("Ridge needs at least 2 features")
    return(NULL)
  }
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
}

RunSVM <- function(Train_set, Train_label, mode, classVar){
  if (ncol(Train_set) == 0) {
    warning("wrong")
    return(NULL)
  }

  Train_set <- as.matrix(Train_set)
  Train_set[is.na(Train_set) | is.infinite(Train_set)] <- 0
  Train_set <- Train_set[, apply(Train_set, 2, function(x) sd(x) > 1e-10)]

  if (ncol(Train_set) == 0) {
    warning("wrong")
    return(NULL)
  }

  data <- as.data.frame(Train_set)
  data[[classVar]] <- Train_label[[classVar]]

  fit = svm(formula = as.formula(paste(classVar, "~", paste(colnames(Train_set), collapse = "+"))),
            data= data, probability = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunglmBoost <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial())

  cvm <- cvrisk(fit, papply = lapply,
                folds = cv(model.weights(fit), type = "kfold"))
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial(),
                  control = boost_control(mstop = max(mstop(cvm), 40)))

  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunRF <- function(Train_set, Train_label, mode, classVar){
  rf_nodesize = 5
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  fit <- rfsrc(formula = formula(paste0(classVar, "~.")),
               data = cbind(Train_set, Train_label[classVar]),
               ntree = 1000, nodesize = rf_nodesize,
               importance = T,
               proximity = T,
               forest = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunXGBoost <- function(Train_set, Train_label, mode, classVar){
  y_label <- as.integer(Train_label[[classVar]])
  indexes = createFolds(y_label, k = 5, list=T)
  CV <- unlist(lapply(indexes, function(pt){
    dtrain = xgb.DMatrix(data = Train_set[-pt, ],
                         label = y_label[-pt])
    dtest = xgb.DMatrix(data = Train_set[pt, ],
                        label = y_label[pt])
    tryCatch({
      bst <- xgb.train(params = list(max_depth = 2, eta = 1,
                                      objective = "reg:logistic"),
                       data = dtrain, nrounds = 10,
                       evals = list(train=dtrain, test=dtest), verbose = 0)
      if (!is.null(bst$evaluation_log) && nrow(bst$evaluation_log) > 0) {
        which.min(bst$evaluation_log$test_logloss)
      } else {
        NA
      }
    }, error = function(e) NA)
  }))

  CV <- CV[!is.na(CV)]
  if (length(CV) == 0) {
    warning("XGBoost CV failed, using default nrounds=10")
    nround <- 10
  } else {
    nround <- as.numeric(names(which.max(table(CV))))
  }
  dtrain_all <- xgb.DMatrix(data = Train_set, label = y_label)
  fit <- xgb.train(params = list(max_depth = 2, eta = 1,
                                  objective = "reg:logistic"),
                   data = dtrain_all, nrounds = nround, verbose = 0)
  attr(fit, "subFeature") <- colnames(Train_set)

  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunSCAD <- function(Train_set, Train_label, mode, classVar){
  if (ncol(Train_set) < 2) {
    warning("SCAD needs at least 2 features")
    return(NULL)
  }
  if (sum(is.na(Train_set)) > 0) {
    Train_set[is.na(Train_set)] <- 0
  }
  cv.fit <- cv.ncvreg(Train_set, Train_label[[classVar]],
                       family = "binomial", penalty = "SCAD")
  fit <- ncvreg(Train_set, Train_label[[classVar]],
                family = "binomial", penalty = "SCAD")
  fit$lambda.min <- cv.fit$lambda.min
  attr(fit, "subFeature") <- colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunMCP <- function(Train_set, Train_label, mode, classVar){
  if (ncol(Train_set) < 2) {
    warning("MCP needs at least 2 features")
    return(NULL)
  }
  if (sum(is.na(Train_set)) > 0) {
    Train_set[is.na(Train_set)] <- 0
  }
  cv.fit <- cv.ncvreg(Train_set, Train_label[[classVar]],
                       family = "binomial", penalty = "MCP")
  fit <- ncvreg(Train_set, Train_label[[classVar]],
                family = "binomial", penalty = "MCP")
  fit$lambda.min <- cv.fit$lambda.min
  attr(fit, "subFeature") <- colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunBoruta <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(as.data.frame(Train_set), Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  boruta.fit <- Boruta(eval(parse(text = paste(classVar, "~."))),
                        data = data, doTrace = 0)
  selected <- names(boruta.fit$finalDecision)[boruta.fit$finalDecision != "Rejected"]
  if (length(selected) == 0) selected <- colnames(Train_set)

  if (mode == "Model"){
    if (length(selected) < 1) selected <- colnames(Train_set)[1:min(5, ncol(Train_set))]
    fit <- glm(formula = Train_label[[classVar]] ~ .,
               family = "binomial",
               data = as.data.frame(Train_set[, selected, drop=FALSE]))
    fit$subFeature <- selected
    return(fit)
  }
  if (mode == "Variable") return(selected)
}

RunmRMR <- function(Train_set, Train_label, mode, classVar){
  n_select <- min(50, ncol(Train_set))
  combined <- cbind(as.data.frame(Train_set), label = as.numeric(Train_label[[classVar]]))
  dd <- mRMR.data(data = combined)
  sel <- mRMR.classic(data = dd, target_indices = ncol(combined), feature_count = n_select)
  selected <- colnames(Train_set)[solutions(sel)[[1]]]

  if (mode == "Model"){
    if (length(selected) < 1) selected <- colnames(Train_set)[1:min(5, ncol(Train_set))]
    fit <- glm(formula = Train_label[[classVar]] ~ .,
               family = "binomial",
               data = as.data.frame(Train_set[, selected, drop=FALSE]))
    fit$subFeature <- selected
    return(fit)
  }
  if (mode == "Variable") return(selected)
}

RunLightGBM <- function(Train_set, Train_label, mode, classVar){
  dtrain <- lgb.Dataset(data = as.matrix(Train_set),
                         label = Train_label[[classVar]])
  params <- list(objective = "binary",
                 metric = "auc",
                 verbose = -1,
                 num_leaves = 15,
                 learning_rate = 0.05)
  fit <- lgb.train(params = params, data = dtrain, nrounds = 100)
  attr(fit, "subFeature") <- colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunCatBoost <- function(Train_set, Train_label, mode, classVar){
  train_pool <- catboost.load_pool(data = as.matrix(Train_set),
                                    label = as.numeric(Train_label[[classVar]]))
  fit <- catboost.train(train_pool,
                        params = list(loss_function = 'Logloss',
                                     logging_level = 'Silent',
                                     iterations = 100))
  attr(fit, "subFeature") <- colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

standarize.fun <- function(indata, centerFlag, scaleFlag) {
  scale(indata, center=centerFlag, scale=scaleFlag)
}

scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)

  if (is.null(cohort)){
    data_list <- list("training" = data)
  }else{
    data_list <- split(as.data.frame(data), cohort)
  }

  if (is.null(centerFlags)){
    centerFlags = TRUE; message("No centerFlags found, set as TRUE")
  }
  if (length(centerFlags)==1){
    centerFlags = rep(centerFlags, length(data_list))
  }
  if (is.null(names(centerFlags))){
    names(centerFlags) <- names(data_list)
  }

  if (is.null(scaleFlags)){
    scaleFlags = TRUE; message("No scaleFlags found, set as TRUE")
  }
  if (length(scaleFlags)==1){
    scaleFlags = rep(scaleFlags, length(data_list))
  }
  if (is.null(names(scaleFlags))){
    names(scaleFlags) <- names(data_list)
  }

  centerFlags <- centerFlags[names(data_list)]
  scaleFlags <- scaleFlags[names(data_list)]

  outdata_list <- mapply(function(indata, centerFlag, scaleFlag) {
    indata[is.na(indata)] <- 0
    scaled <- scale(indata, center=centerFlag, scale=scaleFlag)
    scaled[is.na(scaled) | is.infinite(scaled)] <- 0
    return(scaled)
  }, indata = data_list, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = F)

  outdata <- do.call(rbind, outdata_list)
  outdata <- outdata[samplename, ]

  if (ncol(outdata) != ncol(data)) {
    missing_cols <- setdiff(colnames(data), colnames(outdata))
    if (length(missing_cols) > 0) {
      outdata[, missing_cols] <- 0
    }
  }

  return(outdata)
}

ExtractVar <- function(fit){
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet" = rownames(coef(fit))[which(coef(fit)[, 1]!=0)],
    "glm" = names(coef(fit)),
    "svm.formula" = fit$subFeature,
    "train" = fit$coefnames,
    "glmboost" = names(coef(fit)[abs(coef(fit))>0]),
    "rfsrc" = {
      imp <- fit$importance
      if (is.null(imp) || length(imp) == 0) {
        fit$subFeature
      } else if (is.matrix(imp)) {
        col_idx <- ifelse("all" %in% colnames(imp), "all", 1)
        scores <- imp[, col_idx]
        names(scores) <- rownames(imp)
        names(sort(scores, decreasing = TRUE))
      } else if (all(imp == 0)) {
        fit$subFeature
      } else {
        names(sort(imp, decreasing = TRUE))
      }
    },
    "xgb.Booster" = attr(fit, "subFeature"),
    "ncvreg" = {
      coefs <- coef(fit, lambda = fit$lambda.min)
      names(coefs)[coefs != 0]
    },
    "lgb.Booster" = {
      imp <- tryCatch(lgb.importance(fit, percentage = TRUE), error = function(e) NULL)
      if (!is.null(imp) && nrow(imp) > 0) imp$Feature else attr(fit, "subFeature")
    },
    "catboost.Model" = attr(fit, "subFeature")
  ))

  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  return(Feature)
}

getSubFeature <- function(fit) {
  if (inherits(fit, "xgb.Booster") || inherits(fit, "lgb.Booster") ||
      inherits(fit, "catboost.Model") || inherits(fit, "ncvreg")) {
    attr(fit, "subFeature")
  } else {
    fit$subFeature
  }
}

CalPredictScore <- function(fit, new_data, type = "lp"){
  new_data <- new_data[, getSubFeature(fit), drop = FALSE]
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"       = predict(fit, type = 'response', as.matrix(new_data)),
    "glm"          = predict(fit, type = 'response', as.data.frame(new_data)),
    "svm.formula"  = predict(fit, as.data.frame(new_data), probability = T),
    "train"        = predict(fit, new_data, type = "prob")[[2]],
    "glmboost"     = predict(fit, type = "response", as.data.frame(new_data)),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$predicted[, "1"],
    "xgb.Booster"  = predict(fit, as.matrix(new_data)),
    "ncvreg" = predict(fit, X = as.matrix(new_data),
                        lambda = fit$lambda.min, type = "response"),
    "lgb.Booster" = predict(fit, as.matrix(new_data)),
    "catboost.Model" = catboost.predict(fit,
                          catboost.load_pool(as.matrix(new_data)),
                          prediction_type = "Probability")
  ))
  RS = as.numeric(as.vector(RS))
  names(RS) = rownames(new_data)
  return(RS)
}

PredictClass <- function(fit, new_data){
  new_data <- new_data[, getSubFeature(fit), drop = FALSE]
  label <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"       = predict(fit, type = 'class', as.matrix(new_data)),
    "glm"          = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                            yes = "1", no = "0"),
    "svm.formula"  = predict(fit, as.data.frame(new_data), decision.values = T),
    "train"        = predict(fit, new_data, type = "raw"),
    "glmboost"     = predict(fit, type = "class", as.data.frame(new_data)),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$class,
    "xgb.Booster"  = ifelse(test = predict(fit, as.matrix(new_data))>0.5,
                            yes = "1", no = "0"),
    "ncvreg" = ifelse(test = predict(fit, X = as.matrix(new_data),
                                      lambda = fit$lambda.min, type = "response") > 0.5,
                       yes = "1", no = "0"),
    "lgb.Booster" = ifelse(test = predict(fit, as.matrix(new_data)) > 0.5,
                            yes = "1", no = "0"),
    "catboost.Model" = {
      pred <- catboost.predict(fit,
                               catboost.load_pool(as.matrix(new_data)),
                               prediction_type = "Class")
      as.character(pred)
    }
  ))
  label = as.character(as.vector(label))
  names(label) = rownames(new_data)
  return(label)
}

RunEval <- function(fit,
                    Test_set = NULL,
                    Test_label = NULL,
                    Train_set = NULL,
                    Train_label = NULL,
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    classVar){

  if(!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  }

  if((!is.null(Train_set)) & (!is.null(Train_label))) {
    sf <- getSubFeature(fit)
    new_data <- rbind.data.frame(Train_set[, sf, drop = FALSE],
                                 Test_set[, sf, drop = FALSE])

    if(!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Training"
    }
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    Test_label <- rbind.data.frame(Train_label[,c(cohortVar, classVar)],
                                  Test_label[,c(cohortVar, classVar)])
    Test_label[,1] <- factor(Test_label[,1],
                            levels = c(unique(Train_label[,cohortVar]),
                                       setdiff(unique(Test_label[,cohortVar]),
                                               unique(Train_label[,cohortVar]))))
  } else {
    new_data <- Test_set[, getSubFeature(fit), drop = FALSE]
  }

  RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))

  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS)
  Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
  unlist(lapply(Predict.out, function(data){
    as.numeric(auc(suppressMessages(roc(data[[classVar]], data$RS))))
  }))
}

SimpleHeatmap <- function(Cindex_mat, avg_Cindex,
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5,
                          cluster_columns, cluster_rows){
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)

  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = T, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)

  Heatmap(as.matrix(Cindex_mat), name = "AUC",
          right_annotation = row_ha,
          top_annotation = col_ha,
          col = c("#4195C1", "#FFFFFF", "#CB5746"),
          rect_gp = gpar(col = "black", lwd = 1),
          cluster_columns = cluster_columns, cluster_rows = cluster_rows,
          show_column_names = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)),
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) {
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 10))
          }
  )
}
