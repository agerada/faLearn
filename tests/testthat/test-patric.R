test_patric_path <- "fixtures/test_patric_db.txt"

test_that("multiplication works", {
  expect_s3_class(load_patric_db(test_patric_path), "patric_db")
  expect_s3_class(load_patric_db(), "patric_db")
})
