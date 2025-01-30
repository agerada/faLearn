test_patric_path <- "fixtures/test_patric_db.txt"
test_db <- data.frame(genome_id = 1,
                      genome_name = "test",
                      antibiotic = "test",
                      measurement = 1,
                      measurement_unit = "test",
                      laboratory_typing_method = "test",
                      resistant_phenotype = "test")

test_that("local load patric works", {
  expect_s3_class(load_patric_db(test_patric_path), "patric_db")
  db_df <- read.delim(test_patric_path)
  expect_s3_class(as_patric_db(db_df), "patric_db")
})

test_that("load patric from ftp works", {
  skip_on_cran()
  expect_s3_class(load_patric_db(), "patric_db")
})

test_that("check valid patric_db works",{
  expect_error(check_valid_patric_db(data.frame()))

  expect_no_error(check_valid_patric_db(test_db))

  test_db_invalid <- test_db
  test_db_invalid$antibiotic <- NULL
  expect_error(check_valid_patric_db(test_db_invalid))
})

test_that("check can save patric db", {
  skip_on_cran()
  tmp_path_db <- tempfile(fileext = ".txt")
  download_patric_db(tmp_path_db)
  expect_s3_class(load_patric_db(tmp_path_db), "patric_db")
})

test_that("check can download genomes", {
  skip_on_cran()
  tmp_dir <- tempdir()
  unlink(list.files(tmp_dir, full.names = TRUE, pattern = "*.fna"))
  pull_PATRIC_genomes(tmp_dir,
                      n_genomes = 1)
  all_files <- list.files(tmp_dir, pattern = "*.fna")
  expect_true(length(all_files) == 1)
})
