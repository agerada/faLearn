path_to_test_genomes <- "fixtures/genomes/"
tmp_out_dir <- file.path(tempdir(), "kmer_batch_test/")
if (!dir.exists(tmp_out_dir)) dir.create(tmp_out_dir)

test_that("genomes to libsvm", {
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

test_that("genome to per-file libsvm", {
  target_file <- tempfile(fileext = ".txt")
  source_file <- list.files(path_to_test_genomes, pattern = "*.fna", full.names = TRUE)[1]
  source_dna <- Biostrings::readDNAStringSet(source_file) |> as.character()

  genome_to_libsvm(source_dna,
                   target_path = target_file,
                   k = 3)
  expect_true(file.exists(target_file))

  # try again without overwriting, should raise warning
  expect_message(
    genome_to_libsvm(source_dna,
                     target_path = target_file,
                     k = 3,
                     overwrite = FALSE))

  # try again with overwriting
  expect_no_message(
    genome_to_libsvm(source_dna,
                     target_path = target_file,
                     k = 3,
                     overwrite = TRUE)
  )

})
