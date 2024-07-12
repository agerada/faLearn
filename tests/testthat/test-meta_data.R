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

test_that("test qc_in_range", {
  expect_true(qc_in_range(AMR::as.mic(0.5), 25922, AMR::as.ab("GEN")))
  expect_true(qc_in_range(AMR::as.mic(0.25), 25922, AMR::as.ab("GEN")))
  expect_true(qc_in_range(AMR::as.mic(1.0), 25922, AMR::as.ab("AMK")))
  expect_true(qc_in_range(AMR::as.mic(4.0), 25922, AMR::as.ab("AMK")))

  expect_false(qc_in_range(AMR::as.mic("<0.25"), 25922, AMR::as.ab("GEN")))
  expect_false(qc_in_range(AMR::as.mic(8.0), 25922, AMR::as.ab("AMK")))
  expect_false(qc_in_range(AMR::as.mic(">4.0"), 25922, AMR::as.ab("AMK")))
  expect_false(qc_in_range(AMR::as.mic(0.25), 25922, AMR::as.ab("AMK")))
})

test_that("test qc_on_target", {
  expect_true(qc_on_target(AMR::as.mic(0.5), 25922, AMR::as.ab("GEN")))
  expect_true(qc_on_target(AMR::as.mic(1), 25922, AMR::as.ab("AMK")))
  expect_true(qc_on_target(AMR::as.mic(2), 25922, AMR::as.ab("AMK")))

  expect_false(qc_on_target(AMR::as.mic("<0.5"), 25922, AMR::as.ab("GEN")))
  expect_false(qc_on_target(AMR::as.mic(4), 25922, AMR::as.ab("AMK")))
})

test_that("test standardise_mic", {
  expect_warning(
    standardise_mic(AMR::as.mic(8.0),
                    AMR::as.mic(2),
                    25922,
                    AMR::as.ab("GEN"))
    )

  expect_equal(standardise_mic(AMR::as.mic(8.0),
                               AMR::as.mic(1),
                               25922, AMR::as.ab("GEN")),
               AMR::as.mic(4.0))

  expect_equal(standardise_mic(AMR::as.mic(8.0),
                               AMR::as.mic(2),
                               25922,
                               AMR::as.ab("AMK")),
               AMR::as.mic(8.0))

  expect_equal(standardise_mic(AMR::as.mic(8.0),
                               AMR::as.mic(1),
                               25922,
                               AMR::as.ab("AMK")),
               AMR::as.mic(8.0))

  expect_equal(standardise_mic(AMR::as.mic(8.0),
                               AMR::as.mic(0.5),
                               25922,
                               AMR::as.ab("AMK"),
                               prefer_upper = TRUE),
               AMR::as.mic(32.0))

  expect_equal(standardise_mic(AMR::as.mic(8.0),
                               AMR::as.mic(0.5),
                               25922,
                               AMR::as.ab("AMK"),
                               prefer_upper = FALSE),
               AMR::as.mic(16.0))

  expect_equal(standardise_mic(AMR::as.mic(">8.0"),
                               AMR::as.mic(0.5),
                               25922,
                               AMR::as.ab("AMK"),
                               prefer_upper = T),
               AMR::as.mic(">8.0"))
})
