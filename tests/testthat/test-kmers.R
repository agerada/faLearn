test_that("multiplication works", {
  testthat::expect_true(is.list(kmers("ACTGC")))
})

genomes_dir <- file.path("fixtures", "genomes")

test_that("non canonical kmer counting is correct", {
  genomes_filepaths <- list.files(genomes_dir, full.names = TRUE)
  genomes <- lapply(genomes_filepaths, Biostrings::readDNAStringSet)
  tmp_dir <- tempdir()
  unlink(tmp_dir, recursive = TRUE)
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

test_that("test kmer squeeze", {
  for (k in 3:8) {
    all_mers <- kmers("AAAAAAAAAAAAAAA", k = k, anchor = TRUE)
    all_mers <- all_mers$kmer_string
    all_mers <- paste0(all_mers, collapse = "")
    all_counts <- kmers(all_mers, k = k, anchor = TRUE)
    filtered_mers <- all_counts$kmer_string[all_counts$kmer_value > 0]
    squeezed_mers <- kmers("AAAAAAAAAAAAAAA", k = k, anchor = TRUE, squeeze = TRUE)
    squeezed_mers <- squeezed_mers$kmer_string
    expect_true(all.equal(filtered_mers, squeezed_mers))
  }
})

test_that("test kmer to libsvm squeeze", {
  use_genome <- 6
  for (k in 3:8) {
    # calculate without squeeze first - use genomes_to_kmer_libsvm()
    genomes_filepaths <- list.files(genomes_dir, full.names = TRUE)
    tmp_dir <- tempdir()
    unlink(tmp_dir, recursive = TRUE)
    genomes_to_kmer_libsvm(genomes_dir, tmp_dir, k = k, canonical = TRUE, squeeze = FALSE)
    kmer_files <- list.files(tmp_dir, full.names = TRUE)
    kmer_data <- kmer_files[use_genome] |>
      readLines() |>
      strsplit(" ") |>
      unlist()

    # calculate without squeeze - using kmers() directly
    genome_str <- genomes_filepaths[use_genome] |>
      Biostrings::readDNAStringSet() |>
      as.character()
    kmer_data2 <- kmers(genome_str, k = k, canonical = TRUE, squeeze = FALSE)
    # get the index of the last recorded kmer:
    last_kmer <- max(which(kmer_data2$kmer_value > 0))
    # check this is the same as using the libsvm file
    last_kmer_libsvm <- max(as.numeric(strsplit(kmer_data[length(kmer_data)], ":")[[1]][1]))
    expect_true(last_kmer == last_kmer_libsvm)

    # now calculate with squeeze
    unlink(tmp_dir, recursive = TRUE)
    genomes_to_kmer_libsvm(genomes_dir, tmp_dir, k = k, canonical = TRUE, squeeze = TRUE)
    kmer_files <- list.files(tmp_dir, full.names = TRUE)
    kmer_data <- kmer_files[use_genome] |>
      readLines() |>
      strsplit(" ") |>
      unlist()

    # calculate with squeeze - using kmers() directly
    kmer_data2 <- kmers(genome_str, k = k, canonical = TRUE, squeeze = TRUE)
    # get the index of the last recorded kmer:
    last_kmer <- max(which(kmer_data2$kmer_value > 0))
    # check this is the same as using the libsvm file
    last_kmer_libsvm <- max(as.numeric(strsplit(kmer_data[length(kmer_data)], ":")[[1]][1]))
    expect_true(last_kmer == last_kmer_libsvm)
  }
})

test_that("check that can revert the squeeze", {
    use_genome <- 6
    for (k in 3:8) {
      squeezed_mers <- squeezed_mers(k = k)
      calculated_mers <- kmers("AAAAAAAAAAAAAAA", k = k, anchor = TRUE, squeeze = TRUE)
      expect_true(all.equal(squeezed_mers, calculated_mers$kmer_string))

      genomes_filepaths <- list.files(genomes_dir, full.names = TRUE)
      tmp_dir <- tempdir()
      unlink(tmp_dir, recursive = TRUE)
      genomes_to_kmer_libsvm(genomes_dir, tmp_dir, k = k, canonical = TRUE, squeeze = FALSE)
      kmer_files <- list.files(tmp_dir, full.names = TRUE)
      kmer_data_unsqueezed <- kmer_files[use_genome] |>
        readLines() |>
        strsplit(" ") |>
        unlist()

      # calculate with squeeze
      unlink(tmp_dir, recursive = TRUE)
      genomes_to_kmer_libsvm(genomes_dir, tmp_dir, k = k, canonical = TRUE, squeeze = TRUE)
      kmer_files <- list.files(tmp_dir, full.names = TRUE)
      kmer_data_squeezed <- kmer_files[use_genome] |>
        readLines() |>
        strsplit(" ") |>
        unlist()

      expect_false(isTRUE(all.equal(kmer_data_unsqueezed, kmer_data_squeezed)))

      # now unsqueeze
      squeezed_indices <- kmer_data_squeezed[2:length(kmer_data_squeezed)]
      squeezed_indices <- sapply(squeezed_indices, \(x) strsplit(x, ":")[[1]][[1]], USE.NAMES = F)
      squeezed_indices <- as.numeric(squeezed_indices)
      squeezed_str <- squeezed_index_to_str(squeezed_indices, k = k)

      # compare to unsqueezed data
      unsqueezed_indices <- kmer_data_unsqueezed[2:length(kmer_data_unsqueezed)]
      unsqueezed_indices <- sapply(unsqueezed_indices, \(x) strsplit(x, ":")[[1]][[1]], USE.NAMES = F)
      unsqueezed_indices <- as.numeric(unsqueezed_indices)
      unsqueezed_str <- unsqueezed_index_to_str(unsqueezed_indices, k = k)

      expect_true(all.equal(squeezed_str, unsqueezed_str))
    }
})

test_that("test reverse complement", {
  genome_str <- "ACACTT"
  rev_str <- reverse_complement(genome_str)
  expect_true(all.equal(rev_str, "AAGTGT"))

  genome_str <- "ttt"
  rev_str <- reverse_complement(genome_str)
  expect_true(all.equal(rev_str, "AAA"))
})
