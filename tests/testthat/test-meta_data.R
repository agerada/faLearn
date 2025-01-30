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

  expect_true(qc_in_range(NA, 25922, "GEN"))
  expect_true(qc_in_range(AMR::as.mic(8.0), NA, "GEN"))

  expect_false(qc_in_range(AMR::as.mic("<0.0625"), 25922, "CIP"))
})

test_that("test qc_on_target", {
  expect_true(qc_on_target(AMR::as.mic(0.5), 25922, AMR::as.ab("GEN")))
  expect_true(qc_on_target(AMR::as.mic(1), 25922, AMR::as.ab("AMK")))
  expect_true(qc_on_target(AMR::as.mic(2), 25922, AMR::as.ab("AMK")))

  expect_false(qc_on_target(AMR::as.mic("<0.5"), 25922, AMR::as.ab("GEN")))
  expect_false(qc_on_target(AMR::as.mic(4), 25922, AMR::as.ab("AMK")))

  expect_true(qc_on_target(NA, 25922, AMR::as.ab("AMK")))
  expect_true(qc_on_target(AMR::as.mic("<0.5"), 25922, NA))
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

  suppressWarnings(
    expect_equal(standardise_mic(NA,
                                 AMR::as.mic(0.5),
                                 25922,
                                 AMR::as.ab("AMK"),
                                 prefer_upper = T),
                 AMR::NA_mic_)
  )

  expect_equal(standardise_mic(c(AMR::as.mic("4"),
                                 AMR::as.mic("8")),
                               c(AMR::as.mic("1"),
                                 AMR::as.mic("0.5")),
                               25922,
                               AMR::as.ab("GEN")),
               AMR::as.mic(c("2", "8")))
})

test_that("test compare mic", {
  mic1 <- AMR::as.mic(c("4","8", ">32"))
  mic2 <- AMR::as.mic(c("4", "4", "2"))
  validation <- compare_mic(mic1, mic2)
  expect_s3_class(validation, "mic_validation")
})

test_that("test essential agreement", {
  expect_true(essential_agreement(2, 2))
  expect_true(essential_agreement(2, 4))
  expect_true(essential_agreement(2, ">2"))
  expect_true(essential_agreement("2", "<2"))
  expect_true(essential_agreement(0.002, "<0.002"))
  expect_true(essential_agreement("0.001", "<0.002"))

  expect_false(essential_agreement(2, ">4"))
  expect_false(essential_agreement(256, ">512"))
  expect_false(essential_agreement(128, ">256"))
  expect_false(essential_agreement(512, "128"))
  expect_false(essential_agreement(0.002, "0.006"))

  expect_true(suppressWarnings(is.na(essential_agreement(2, ">1024"))))
  expect_warning(essential_agreement(2, ">1024"))
  expect_true(is.na(essential_agreement(2, "<0.001")))

  left_check <- c(4, 2, 7, ">32")
  right_check <- c(4, 1, 1, "32")
  expect_equal(essential_agreement(left_check, right_check),
               c(TRUE, TRUE, FALSE, TRUE))
})

test_that("test mic_censor", {
  censor_rules <- list("B_ESCHR_COLI" = list(
    "AMK" = list(min = 2, max = 32),
    "CHL" = list(min = 4, max = 64),
    "GEN" = list(min = 1, max = 16),
    "CIP" = list(min = 0.015, max = 4),
    "MEM" = list(min = 0.016, max = 16),
    "AMX" = list(min = 2, max = 64),
    "AMC" = list(min = 2, max = 64),
    "FEP" = list(min = 0.5, max = 64),
    "CAZ" = list(min = 1, max = 128),
    "TGC" = list(min = 0.25, max = 1)
  ))

  suppressMessages(
    {
      expect_s3_class(mic_censor(AMR::as.mic(0.5), "TGC", "B_ESCHR_COLI", censor_rules), "mic")

      expect_equal(mic_censor(AMR::as.mic(128), "AMK", "B_ESCHR_COLI", censor_rules),
                   AMR::as.mic(">32"))
      expect_equal(mic_censor(AMR::as.mic(">128"), "AMK", "B_ESCHR_COLI", censor_rules),
                   AMR::as.mic(">32"))
      expect_equal(mic_censor(">128", "AMK", "B_ESCHR_COLI", censor_rules),
                   AMR::as.mic(">32"))

      # test vectorized
      expect_equal(mic_censor(c(AMR::as.mic(128), AMR::as.mic(0.5)), "AMK", "B_ESCHR_COLI", censor_rules),
                   c(AMR::as.mic(">32"), AMR::as.mic("<=2")))
    }
  )
})

test_that("test compare_mic", {
  gs <- c("0.5", "4", ">8", "2")
  test <- c("0.5", "8", "2", "0.5")
  expect_s3_class(compare_mic(gs, test), "mic_validation")
  expect_equal(summary(compare_mic(gs, test))[["EA_pcent"]], 0.5)

  ab <- c("amoxicillin", "amoxicillin", "gentamicin", "gentamicin")
  mo <- "Escherichia coli"
  suppressMessages(
    {
      expect_s3_class(compare_mic(gs, test, ab, mo), "mic_validation")
      expect_equal(summary(compare_mic(gs, test, ab, mo))[[1, "EA_pcent"]], 1)
      expect_equal(summary(compare_mic(gs, test, ab, mo))[[2, "EA_pcent"]], 0)
    }
  )

})

test_that("test compare_sir", {
  gs <- c("4", "16", ">8", "0.5")
  test <- c("0.5", "8", "0.25", ">64")
  ab <- "GEN"
  mo <- "Escherichia coli"
  val <- suppressMessages(compare_mic(gs, test, ab, mo))
  expect_s3_class(val, "mic_validation")
  sum_val <- summary(val)
  expect_equal(sum_val$EA_pcent[sum_val$ab == "GEN"], 0.25)
  expect_equal(sum_val$very_major_error_pcent[sum_val$ab == "GEN"], 50)
})

test_that("test fall back to ECOFFS in compare_mic", {
  # single
  val <- suppressMessages(
    compare_mic(AMR::as.mic(0.5),
                     AMR::as.mic(0.5),
                     "CHL",
                     "Escherichia coli",
                     accept_ecoff = TRUE)
  )
  expect_equal(val$error[val$ab == "CHL"], factor(NA,
                                                  levels = c("M", "vM", "m")))

  val <- suppressMessages(
    compare_mic(AMR::as.mic(64),
                     AMR::as.mic(0.5),
                     "CHL",
                     "Escherichia coli",
                     accept_ecoff = TRUE)
  )
  expect_equal(val$error[val$ab == "CHL"], factor("vM",
                                                  levels = c("M", "vM", "m")))

  # multiple
  gs <- c("0.5", "4", ">64", "2")
  test <- c("0.5", "8", "0.25", ">64")
  ab <- "CHL"
  mo <- "Escherichia coli"
  val <- suppressMessages(
    compare_mic(gs, test, ab, mo, accept_ecoff = TRUE))
  expect_s3_class(val, "mic_validation")
  sum_val <- summary(val)
  expect_equal(sum_val$EA_pcent[sum_val$ab == "CHL"], 0.5)
  expect_equal(sum_val$very_major_error_pcent[sum_val$ab == "CHL"], 25)
  expect_equal(sum_val$minor_error_pcent[sum_val$ab == "CHL"], 0)
  expect_equal(sum_val$major_error_pcent[sum_val$ab == "CHL"], 25)
})

test_that("test mic_uncensor", {
  expect_equal(mic_uncensor(">16", method = "scale"), AMR::as.mic(32))
  expect_equal(mic_uncensor(">16", method = "simple"), AMR::as.mic(16))

  expect_equal(mic_uncensor(c("0.5", "<1"), method = "scale"), AMR::as.mic(c("0.5", "0.5")))
  expect_equal(mic_uncensor(c("0.5", "<1"), method = "simple"), AMR::as.mic(c("0.5", "1")))

  ab <- "GEN"
  mo <- "Escherichia coli"

  expect_true(as.numeric(mic_uncensor(">16", method = "bootstrap", ab = ab, mo = mo)) > 16)
  mic1 <- c("0.5", "4", ">8", "2", "<=1", ">16")
  uncensored_mic1 <- mic_uncensor(mic1, method = "bootstrap", ab = ab, mo = mo)
  expect_true(as.numeric(uncensored_mic1[3]) > 8)

  expect_equal(mic_uncensor(NA, method = "bootstrap", ab = ab, mo = mo), AMR::NA_mic_)

  expect_warning(mic_uncensor("<0.004", method = "bootstrap", ab = ab, mo = mo))
})

test_that("test force_mic", {
  mode <- "closest"
  expect_equal(force_mic("2", method = mode),
               "2")
  expect_equal(force_mic("2.1", method = mode),
               "2")
  expect_equal(force_mic("2.5", method = mode),
               "2")
  expect_equal(force_mic("2.9", method = mode),
               "2")
  expect_equal(force_mic("3", method = mode, prefer = "max"),
               "4")

  mode <- "round up"
  expect_equal(force_mic("2", method = mode),
               "2")
  expect_equal(force_mic("2.1", method = mode),
               "4")
  expect_equal(force_mic("2.5", method = mode),
               "4")
  expect_equal(force_mic("2.9", method = mode),
               "4")
  expect_equal(force_mic("3", method = mode, prefer = "max"),
               "4")

})

test_that("test mic_s_breakpoint", {
  suppressMessages(
    expect_s3_class(mic_s_breakpoint(mo = "E. coli", ab = "TZP", accept_ecoff = FALSE), "mic")
  )
  suppressMessages(
    expect_s3_class(mic_s_breakpoint(mo = "E. coli", ab = "CHL", accept_ecoff = TRUE), "mic")
  )

  expect_error(
    suppressMessages(
      mic_s_breakpoint(mo = "E. coli", ab = "CHL", accept_ecoff = FALSE)
      )
    )
})

test_that("test mic_r_breakpoint", {
  suppressMessages(
    expect_s3_class(mic_r_breakpoint(mo = "E. coli", ab = "TZP", accept_ecoff = FALSE), "mic")
  )
  suppressMessages(
    expect_s3_class(mic_r_breakpoint(mo = "E. coli", ab = "CHL", accept_ecoff = TRUE), "mic")
  )

  expect_error(
    suppressMessages(
      mic_r_breakpoint(mo = "E. coli", ab = "CHL", accept_ecoff = FALSE)
      )
    )
})

test_that("test fill dilution levels", {
  test_range <- AMR::as.mic(c("0.5", "4", "32"))
  test_range <- droplevels(test_range, as.mic = TRUE)
  filled_range <- fill_dilution_levels(test_range)
  expect_true(all.equal(as.character(test_range),
                        as.character(filled_range)))
  expect_false(all(levels(filled_range) %in% levels(test_range)))
  expect_true(length(levels(filled_range)) >= length(levels(test_range)))
})
