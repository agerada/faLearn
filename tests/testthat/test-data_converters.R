test_that("multiplication works", {
  dummy_genomes_path <- test_path("fixtures", "genomes")
  paths <- genome_paths_from_dir(dummy_genomes_path)
  ids <- strip_filename(paths)
  expect_contains(ids, c("001", "002"))
})
