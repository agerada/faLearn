# MIC 1.1.0

One note for `plot.mic_validation` examples CPU >10s.
This warning remains when wrapping in /dontrun{} and even when example removed
altogether.

## Resubmission

Only change is removal of space from doi entry in description.

## Previous Resubmission

This is a resubmission with the below changes, following CRAN review by
Konstanze Lauseker. Thank you for the feedback.

* Removed single quotes from DESCRIPTION file;
* References (DOI or URL) added to the Description field of DESCRIPTION file;
* Return value of the following functions have been clarified (or function
refactored to have an explicit return value):

    - combined_file_system
    - download_patric_db
    - genomes_to_kmer_libsvm
    - move_files
    - pull_PATRIC_genomes

* \dontrun removed from:

    - genomes_to_kmer_libsvm
    - load_patric_db
    - download_patric_db
    - pull_PATRIC_genomes

* removed use of getwd() as default argument in data_converters.R

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

Functions `download_patric_db` and `pull_PATRIC_genomes` examples not run (internet access required).
