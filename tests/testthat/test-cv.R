test_that("test cv", {
  mat <- matrix(rnorm(1000), ncol = 10)
  label <- rnorm(100)
  dtrain <- xgboost::xgb.DMatrix(data = mat, label = label)
  model <- xgboost::xgb.cv(data = dtrain,
                           nrounds = 10,
                           nfold = 5,
                           verbose = FALSE)
  expect_s3_class(model, "xgb.cv.synchronous")

  # using the custom low memory cv function
  model_lowmem <- xgb.cv.lowmem(data = dtrain,
                                nfold = 3,
                                nrounds = 10,
                                verbose = FALSE)
  expect_s3_class(model_lowmem, "xgb.cv.synchronous")

  # with custom folds (similar to generating using caret::createFolds())
  folds <- list()
  folds$Fold1 <- 1:20
  folds$Fold2 <- 21:40
  folds$Fold3 <- 41:60
  folds$Fold4 <- 61:80
  folds$Fold5 <- 81:100
  model <- xgboost::xgb.cv(data = dtrain,
                           nrounds = 10,
                           folds = folds,
                           verbose = FALSE)
  expect_s3_class(model, "xgb.cv.synchronous")

  model_lowmem <- xgb.cv.lowmem(data = dtrain,
                                folds = folds,
                                nrounds = 10,
                                verbose = FALSE,
                                prediction = TRUE)
  expect_s3_class(model_lowmem, "xgb.cv.synchronous")

  # check members
  expected_fields <- c("call",
                       "params",
                       "callbacks",
                       "evaluation_log",
                       "niter",
                       "nfeatures",
                       "folds")

  expect_contains(names(model), expected_fields)
  expect_contains(names(model_lowmem), expected_fields)

  # using some additional callbacks, other fields are added
  model <- xgboost::xgb.cv(data = dtrain,
                           prediction = TRUE,
                           nrounds = 10,
                           folds = folds,
                           verbose = FALSE,
                           callbacks = list(
                             xgboost::cb.cv.predict(save_models = TRUE)
                           )
  )
  expect_s3_class(model, "xgb.cv.synchronous")
  expect_contains(names(model), c(expected_fields,
                                  "pred",
                                  "models"))

  # to save models in lowmem cv, use the save_models argument
  model_lowmem <- xgb.cv.lowmem(data = dtrain,
                                folds = folds,
                                nrounds = 10,
                                verbose = FALSE,
                                prediction = TRUE,
                                save_models = TRUE)

  expect_s3_class(model_lowmem, "xgb.cv.synchronous")
  expect_contains(names(model_lowmem), c(expected_fields,
                                         "pred",
                                         "models"))

  expect_equal(dim(model$pred), dim(model_lowmem$pred))
  expect_equal(dim(model$models), dim(model_lowmem$models))

  # using one fold with standard cv function
  # (i.e., similar to test train split) fails
  folds <- list()
  folds$Test <- 1:20
  expect_error(
    xgboost::xgb.cv(data = dtrain,
                             nrounds = 10,
                             folds = folds,
                             verbose = FALSE)
  )

  # using one fold with lowmem cv function works
  model_lowmem <- xgb.cv.lowmem(data = dtrain,
                                nfold = .8,
                                nrounds = 10,
                                verbose = FALSE,
                                prediction = TRUE,
                                save_models = TRUE)
  expect_s3_class(model_lowmem, "xgb.cv.synchronous")
  expect_contains(names(model_lowmem), c(expected_fields,
                                         "pred",
                                         "models"))

})
