path_to_test_genomes <- "fixtures/genomes/"
tmp_out_dir <- file.path(tempdir(), "kmer_batch_test/")
if (!dir.exists(tmp_out_dir)) dir.create(tmp_out_dir)

test_that("libsvm from genomes", {
  unlink(tmp_out_dir, recursive = TRUE)
  suppressMessages(
    genomes_to_kmer_dataset(path_to_test_genomes,
                            output_dir = tmp_out_dir,
                            cores = 1,
                            kmers = 3,
                            split = 1,
                            anchor = FALSE,
                            simplify = FALSE,
                            drop_n = TRUE,
                            integer_index = TRUE,
                            formats = "libsvm")
  )
  libsvm_path <- file.path(tmp_out_dir, "libsvm/")
  dir.create(libsvm_path)
  libsvm_files <- list.files(tmp_out_dir, pattern = "*.txt")
  lapply(libsvm_files,
         function(x) {
           file.rename(file.path(tmp_out_dir, x),
                       file.path(libsvm_path, x))
         })
  x <- xgboost::xgb.DMatrix(libsvm_path)
  expect_s3_class(x, "xgb.DMatrix")
})
