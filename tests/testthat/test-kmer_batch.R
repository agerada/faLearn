path_to_test_genomes <- "fixtures/genomes/"
tmp_out_dir <- file.path(tempdir(), "kmer_batch_test/")
if (!dir.exists(tmp_out_dir)) dir.create(tmp_out_dir)

test_that("genomes to individual libsvm", {
  unlink(tmp_out_dir, recursive = TRUE)
  genomes_to_kmer_libsvm(source_dir = path_to_test_genomes,
                         target_dir = tmp_out_dir)
  libsvm_files <- list.files(tmp_out_dir, pattern = "*.txt")
  original_genomes <- list.files(path_to_test_genomes, pattern = "*.fna")
  expect_equal(length(libsvm_files),
               length(original_genomes))
  expect_equal(basename(tools::file_path_sans_ext(libsvm_files)),
               basename(tools::file_path_sans_ext(original_genomes)))
})
