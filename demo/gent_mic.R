library(xgboost)
paths_data <- train_test_filesystem("~/Downloads/patric/esco/9/", ".txt")
paths_ids <- train_test_filesystem("~/Downloads/patric/esco/9/", ".csv",
                               train_folder = "train_id",
                               test_folder = "test_id")
train <- xgboost::xgb.DMatrix(paths_data[["train"]])
test <- xgboost::xgb.DMatrix(paths_data[["test"]])
train_ids <- read_csv(list.files(paths_ids["train_id"], full.names = T, pattern = "*.csv"), col_types = cols(.default = "c"))
test_ids <- read_csv(list.files(paths_ids["test_id"], full.names = T, pattern = "*.csv"), col_types = cols(.default = "c"))

meta_data <- load_patric_db("tmp/PATRIC_genomes_AMR.txt")
tidy_meta <- tidy_patric_meta_data(meta_data)

# gentamicin train
labels <- get_mic(tidy_meta, train_ids$genome_id, "genome_id", "mic_GEN")
labels <- labels |>
  mic_uncensor() |>
  as.numeric() |>
  log2()

weights <- ifelse(is.na(labels), 0, 1)
labels[is.na(labels)] <- 0
xgboost::setinfo(train, "label", labels)
xgboost::setinfo(train, "weight", weights)

gent_model <- xgb.train(data = train,
                        nrounds = 100,
                        nthread = 9,
                        max.depth = 6,
                        eta = 0.3)
predictions <- 2 ** predict(gent_model, test)
test_mics <- get_mic(tidy_meta, test_ids$genome_id, "genome_id", "mic_GEN")
predictions[is.na(test_mics)] <- NA
