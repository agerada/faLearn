#' Low memory cross-validation wrapper for XGBoost
#'
#' @param params parameters for xgboost
#' @param data DMatrix or matrix
#' @param nrounds number of training rounds
#' @param nfold number of folds, or if < 1 then the proportion will be used
#' as the training split in a train-test split
#' @param label data labels (alternatively provide with DMatrix)
#' @param missing handling of missing data (see xgb.cv)
#' @param prediction return predictions
#' @param metrics evaluation metrics
#' @param obj custom objective function
#' @param feval custom evaluation function
#' @param stratified whether to use stratified folds
#' @param folds custom folds
#' @param train_folds custom train folds
#' @param verbose verbosity level
#' @param print_every_n print every n iterations
#' @param early_stopping_rounds early stopping rounds (applied to each fold)
#' @param maximize whether to maximize the evaluation metric
#' @param save_models whether to save the models
#' @param ... additional arguments passed to xgb.train
#'
#' @return xgb.cv.synchronous object
#'
#' @description
#' This function performs similar operations to xgboost::xgb.cv, but with the
#' operations performed in a memory efficient manner. Unlike xgboost::xgb.cv,
#' this version does not load all folds into memory from the start. Rather it
#' loads each fold into memory sequentially, and trains trains each fold using
#' xgboost::xgb.train. This allows larger datasets to be cross-validated.
#'
#' The main disadvantage of this function is that it is not possible to perform
#' early stopping based the results of all folds. The function does accept an
#' early stopping argument, but this is applied to each fold separately. This
#' means that different folds can (and should be expected to) train for a
#' different number of rounds.
#'
#' This function also allows for a train-test split (as opposed to multiple)
#' folds. This is done by providing a value of less than 1 to nfold, or a list
#' of 1 fold to folds. This is not possible with xgboost::xgb.cv, but can be
#' desirable if there is downstream processing that depends on an
#' xgb.cv.synchromous object (which is the return object of both this function
#' and xgboost::xgb.cv).
#'
#' Otherwise, where possible this function tries to return the same data
#' structure as xgboost::xgb.cv, with the exception of callbacks (not supported
#' as a field within the return object). To save models, use the save_models
#' argument, rather than the cb.cv.predict(save_models = TRUE) callback.
#'
#' @export
#' @examples
#' train <- list(data = matrix(rnorm(20), ncol = 2),
#'              label = rbinom(10, 1, 0.5))
#' dtrain <- xgboost::xgb.DMatrix(train$data, label = train$label, nthread = 1)
#' cv <- xgb.cv.lowmem(data = dtrain,
#'                    params = list(objective = "binary:logistic"),
#'                    nrounds = 2,
#'                    nfold = 3,
#'                    prediction = TRUE,
#'                    nthread = 1)
#' cv
xgb.cv.lowmem <- function(params = list(),
                          data,
                          nrounds,
                          nfold,
                          label = NULL,
                          missing = NA,
                          prediction = FALSE,
                          metrics = list(),
                          obj = NULL,
                          feval = NULL,
                          stratified = TRUE,
                          folds = NULL,
                          train_folds = NULL,
                          verbose = 1,
                          print_every_n = 1L,
                          early_stopping_rounds = NULL,
                          maximize = NULL,
                          save_models = FALSE,
                          ...) {
  output <- list()

  if (!inherits(data, "xgb.DMatrix")) {
    tryCatch({
      data <- xgboost::xgb.DMatrix(data = data, label = label, missing = missing)
    }, error = function(e) {
      stop("data must be an xgb.DMatrix or a matrix")
    })
  }

  if (!is.null(train_folds)) {
    if (length(train_folds) != nfold) stop("Length of train_folds must be equal to nfold")
  }

  if (!is.null(label)) {
    xgboost::setinfo(data, "label", label)
  }

  if (is.null(folds)) {
    if (nfold == 1) {
      stop("nfold must be greater than 1 (for cross-validation), or less than 1
           (for train-test split)")
    }

    if (stratified) {
      stratifying_labels <- xgboost::getinfo(data, "label")
    } else {
      stratifying_labels <- rep(0, nrow(data))
    }

    if (nfold < 1) {
      folds <- caret::createDataPartition(y = stratifying_labels,
                                          times = 1,
                                          p = nfold,
                                          list = TRUE)
      # createDataPartition returns train indices, so we need to invert them
      folds <- lapply(folds, \(x) setdiff(seq_len(nrow(data)), x))
    } else {
      folds <- caret::createFolds(stratifying_labels,
                                  k = nfold,
                                  list = TRUE,
                                  returnTrain = FALSE)
    }

  }

  populated_fields <- c("call",
                        "params",
                        "nfeatures")

  models <- vector("list", length(folds))

  for (i in seq_along(folds)) {
    test_indices <- as.integer(folds[[i]])
    dtest <- xgboost::slice(data, test_indices)

    if (!is.null(train_folds)) {
      train_indices <- as.integer(train_folds[[i]])
    } else {
      train_indices <- setdiff(seq_len(nrow(data)), test_indices)
    }

    dtrain <- xgboost::slice(data, train_indices)

    model <- xgboost::xgb.train(params = params,
                                data = dtrain,
                                nrounds = nrounds,
                                watchlist = list(train = dtrain, eval = dtest),
                                obj = obj,
                                feval = feval,
                                verbose = verbose,
                                print_every_n = print_every_n,
                                early_stopping_rounds = early_stopping_rounds,
                                maximize = maximize,
                                ...)
    models[[i]] <- model
    if (prediction) {
      output$pred <- c(output$pred, stats::predict(model, dtest))
    }

    for (i in seq_along(populated_fields)) {
      field <- populated_fields[[i]]
      if (is.null(output[[field]])) {
        output[[field]] <- model[[field]]
      } else {
        stopifnot(identical(output[[field]], model[[field]]))
      }
    }

  }

  # reorder predictions to match the order of the labels
  if (prediction) {
    # first, reorganise to matrix if softmax
    if (model[["params"]]["objective"] == "multi:softprob") {
      # convert to matrix of class probabilities
      output$pred <- matrix(output$pred,
                            ncol = model[["params"]][["num_class"]],
                            byrow = TRUE)
      # reorder matrix rows to match the order of the labels
      output$pred <- output$pred[order(unlist(folds)), ]
    } else {
      output$pred <- output$pred[order(unlist(folds))]
    }
  }

  output$niter <- max(sapply(models, \(x) x$niter))

  metric_names <- names(models[[1]]$evaluation_log)
  train_metric_name <- metric_names[2]
  test_metric_name <- metric_names[3]

  train_rmse <- lapply(models, \(x) {
    x[['evaluation_log']][[train_metric_name]]
    })
  train_rmse <-
    (lapply(train_rmse, "length<-", max(lengths(train_rmse))))
  train_rmse <- do.call(cbind, train_rmse)
  train_rmse_means <- rowMeans(train_rmse, na.rm = TRUE)

  train_rmse_std <- apply(train_rmse, 1, \(x) stats::sd(x, na.rm = TRUE))

  test_rmse <- lapply(models, \(x) {
    x[['evaluation_log']][[test_metric_name]]
  })

  test_rmse <-
    (lapply(test_rmse, "length<-", max(lengths(test_rmse))))
  test_rmse <- do.call(cbind, test_rmse)
  test_rmse_means <- rowMeans(test_rmse, na.rm = TRUE)

  test_rmse_std <- apply(test_rmse, 1, \(x) stats::sd(x, na.rm = TRUE))

  output$evaluation_log <- data.frame(
    iter = seq_len(output$niter)
  )
  output$evaluation_log[paste0(train_metric_name, "_mean")] <- train_rmse_means
  output$evaluation_log[paste0(train_metric_name, "_std")] <- train_rmse_std

  # replace any "eval" with "test" in the test metric name
  test_metric_name <- gsub("eval", "test", test_metric_name)
  output$evaluation_log[paste0(test_metric_name, "_mean")] <- test_rmse_means
  output$evaluation_log[paste0(test_metric_name, "_std")] <- test_rmse_std

  rownames(output$evaluation_log) <- seq_len(nrow(output$evaluation_log))

  if (save_models) {
    output$models <- models
  }

  output$folds <- folds
  output$callbacks <- list()

  class(output) <- "xgb.cv.synchronous"
  invisible(output)
}
