# MIC (development version)

* `compare_mic` is faster when only one `ab` is provided
* `subset` S3 method added for `mic_validation`
* `plot.mic_validation` now properly matches dilutions on the lower end of the
scale
* `droplevels.mic_validation` method added that allows unnecessary MIC levels
in a validation object to be dropped
* `tidy_patric_meta_data` was missing MICs when `laboratory_typing_method` was
"MIC" (see 81c69f08)
* `pull_patric_genomes` now takes an `ab` argument to only download strains
where the specified antibiotic was tested

# MIC 1.0.2

* CRAN resubmission with minor changes.

# MIC 1.0.1

* CRAN resubmission with minor changes to DESCRIPTION file.

# MIC 1.0.0

* Initial release and CRAN submission.

