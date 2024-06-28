test_that("test get_mic", {
  meta_dt <- data.frame(genome_id = c("a1",
                                      "b2",
                                      "f6",
                                      "c1"),
                        gent_mic = suppressWarnings(
                          AMR::as.mic(c("4", "<0.5", "1.5", "NA"))
                          )
                        )
  ids <- c("c1", "b2", "f6", "n4")
  expected_mics <- suppressWarnings(
    AMR::as.mic(c("NA", "<0.5", "1.5", "NA"))
    )
  expect_equal(get_mic(meta_dt,
                       ids,
                       ab_col = "gent_mic",
                       id_col = "genome_id",
                       simplify = TRUE),
               expected_mics)

})
