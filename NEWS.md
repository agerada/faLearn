# MIC (development version)

* Future versions of `MIC` may remove the genomic functions of the package, to
focus on MIC analysis. Any removed functions will be ported to a new package
that focuses on machine learning with genomic data, so the only change for users
will be to install and load a new package.

* Functionality added to handle <= and => in MIC values, distinguishing them from
< and >. In general, inhibitory MIC assays report <= and >. The `leq` and `geq`
arguments for `force_mic` can be used to enforce these. The `essential_agreement`
and `compare_mic` functions now have `tolerate_leq` and `tolerate_geq` arguments.
When all settings are kept at default values, the functions assume typical censoring,
i.e., <= and >.

# MIC 1.1.0

*  revamp of the `essential_agreement` function to allow a more robust, flexible,
and explicit approach to dealing with censored values. Now, `essential_agreement`
(and `compare_mic`) have `tolerate_censoring` and `tolerate_matched_censoring`
arguments to control how censored values are handled. The default values should
be appropriate for most situations where the user is comparing an investigational
method to a gold standard method.
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
* kmer counting is now case insensitive
* compatibility with `AMR` v3.0

# MIC 1.0.2

* CRAN resubmission with minor changes.

# MIC 1.0.1

* CRAN resubmission with minor changes to DESCRIPTION file.

# MIC 1.0.0

* Initial release and CRAN submission.

