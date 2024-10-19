test_that("multiplication works", {
  testthat::expect_true(is.list(kmers("ACTGC")))
})

genomes_dir <- file.path("fixtures", "genomes")

test_that("non canonical kmer counting is correct", {
  genomes_filepaths <- list.files(genomes_dir, full.names = TRUE)
  genomes <- lapply(genomes_filepaths, Biostrings::readDNAStringSet)
  tmp_dir <- tempdir()
  genomes_to_kmer_libsvm(genomes_dir, tmp_dir, k = 2, canonical = FALSE)
  kmer_files <- list.files(tmp_dir, full.names = TRUE)
  # read last file and print
  kmer_data <- readLines(kmer_files[4])
  kmer_data <- unlist(strsplit(kmer_data, " "))
  assert_str <- c("0", "1:3", "2:2", "3:1", "4:2",
                  "5:1", "7:4", "8:6",
                  "9:1", "10:4",
                  "13:2", "14:6", "16:1")
  expect_true(all.equal(kmer_data, assert_str))
})

test_that("canonical kmer counting is correct", {
  genome_str <- "ACACTT"
  kmers <- kmers(genome_str, k = 2, canonical = TRUE, anchor = FALSE)
  expect_true(all.equal(kmers, list(kmer_string = c("AA", "AC", "AG", "CA"),
                                    kmer_value = c(1, 2, 1, 1))))
  })
