path_to_test_libsvm <- "fixtures/genomes_libsvm/"
tmp_out_dir <- file.path(tempdir(), "data_converters_tests/")

test_that("test combine files", {
  unlink(tmp_out_dir, recursive = TRUE)
  split_and_combine_files(path_to_files = path_to_test_libsvm,
                          train_target_path = file.path(tmp_out_dir, "train.txt"),
                          test_target_path = file.path(tmp_out_dir, "test.txt"),
                          names_backup = file.path(tmp_out_dir, "names.csv"),
                          split = 1)
  expect_true(file.exists(file.path(tmp_out_dir, "train.txt")))

  # test progress bar
  unlink(tmp_out_dir, recursive = TRUE)
  progressr::with_progress({
    split_and_combine_files(path_to_files = path_to_test_libsvm,
                            train_target_path = file.path(tmp_out_dir, "train.txt"),
                            test_target_path = file.path(tmp_out_dir, "test.txt"),
                            names_backup = file.path(tmp_out_dir, "names.csv"),
                            split = 1)
  })

  # list files instead
  unlink(tmp_out_dir, recursive = TRUE)
  filenames <- list.files(path_to_test_libsvm, full.names = TRUE,
                          pattern = "*.txt")
  progressr::with_progress({
  split_and_combine_files(path_to_files = filenames,
                          train_target_path = file.path(tmp_out_dir, "train.txt"),
                          test_target_path = file.path(tmp_out_dir, "test.txt"),
                          names_backup = file.path(tmp_out_dir, "names.csv"),
                          split = 0)
  })
})
