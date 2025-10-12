
<!-- README.md is generated from README.Rmd. Please edit that file -->

# faLearn

<!-- badges: start -->

[![R-CMD-check](https://github.com/agerada/faLearn/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/agerada/faLearn/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Introduction

`faLearn` provides utilities to transform whole-genome sequence data
into feature representations suitable for machine learning, with a focus
on k-mer based features and XGBoost-compatible input (libsvm). The
package is intended for workflows that start from raw or assembled
genome FASTA/FNA files and end with model-ready data.

This repository was migrated from a prior package; MIC-specific
functions have been removed — the current scope is machine learning with
genome sequence data.

## Main features

- Convert individual genomes (FASTA/FNA) into k-mer counts.
- Export k-mer counts in XGBoost-friendly libsvm format (one .txt per
  genome).
- Fast k-mer counting implemented in C++ (Rcpp) for performance.
- Helpers to process directories of genomes in parallel (via
  future.apply) and to split/combine libsvm files for
  training/validation workflows.
- Memory-efficient XGBoost cross-validation using
  [xgb.cv.lowmem](man/xgb.cv.lowmem.Rd).

## Installation

Install from GitHub (development version):

``` r
# install.packages("remotes")
remotes::install_github("agerada/faLearn")
```

## Quick examples

[`libsvm`]( "https://xgboost.readthedocs.io/en/stable/tutorials/input_format.html")
is the preferred input format for XGBoost. Data are represented as a
sparce matrix in a special text format.

To convert a directory of genomes (.fna or .fasta) directly to libsvm
files (one line per genome):

``` r
library(faLearn)
# make 5 tiny genomes for the example
tmp_dir <- tempfile("genomes")
dir.create(tmp_dir)
set.seed(1)
for (i in 1:5) {
  seq <- paste0(sample(c("A","C","G","T"), 200, replace = TRUE), collapse = "")
  writeLines(c(paste0(">gen", i), seq), file.path(tmp_dir, paste0(i, ".fna")))
}

tmp_out <- file.path(tempdir(), "kmers")
unlink(tmp_out, recursive = TRUE)

future::plan(future::sequential) # or multisession, etc.

# convert the first genome file to libsvm using the single-file helper
genome_path <- list.files(tmp_dir, pattern = "\\.fna$", full.names = TRUE)[1]
dir.create(tmp_out, recursive = TRUE, showWarnings = FALSE)
target_file <- file.path(normalizePath(tmp_out), paste0(tools::file_path_sans_ext(basename(genome_path)), ".txt"))
progressr::with_progress({
  genome_to_libsvm(as.character(Biostrings::readDNAStringSet(genome_path)),
                   target_file,
                   k = 5)
})
#> [1] TRUE

list.files(tmp_out)
#> [1] "1.txt"
readLines(list.files(tmp_out, full.names = TRUE)[1])
#> [1] "0 2:1 5:2 6:1 7:1 10:1 15:1 17:1 19:1 22:3 27:1 37:1 42:1 43:2 53:1 56:1 57:1 58:1 66:1 74:1 75:2 78:1 83:2 84:1 85:1 86:3 88:1 89:1 90:1 95:1 107:2 113:1 129:2 133:1 134:1 141:1 143:1 145:1 148:1 155:2 157:1 159:1 162:1 166:2 167:1 170:2 174:1 182:1 198:2 201:2 202:1 205:3 206:1 212:3 213:1 221:2 222:1 226:1 230:1 231:1 235:1 238:1 242:1 243:1 245:1 253:1 262:1 267:1 277:1 278:1 283:1 290:1 295:1 297:1 298:2 299:1 306:1 307:2 309:2 313:1 315:1 317:1 326:1 329:1 331:3 333:1 339:1 341:1 342:1 343:3 345:1 354:1 355:1 358:1 361:2 362:2 381:1 382:1 390:1 393:1 407:1 413:1 417:2 418:1 422:2 429:1 430:1 441:1 457:2 458:2 461:1 466:1 481:1 489:1 494:2 497:1 502:1 506:1 509:1 513:1 514:2 518:1 526:1 533:1 541:1 545:1 562:1 570:2 573:1 598:1 602:1 613:1 614:1 617:1 633:1 641:1 642:1 645:1 653:1 662:2 677:1 685:1 689:1 693:1 697:1 722:1 733:1 745:1 765:1 789:1 801:1 805:1 825:2 885:1 901:1 929:1 961:3 "
```

Count k-mers for a single sequence string using the fast C++
implementation:

``` r
kmers("ATCGATCGA", k = 3)
#> $kmer_string
#>  [1] "AAA" "AAC" "AAG" "AAT" "ACA" "ACC" "ACG" "ACT" "AGA" "AGC" "AGG" "AGT"
#> [13] "ATA" "ATC" "ATG" "ATT" "CAA" "CAC" "CAG" "CAT" "CCA" "CCC" "CCG" "CCT"
#> [25] "CGA" "CGC" "CGG" "CGT" "CTA" "CTC" "CTG" "CTT" "GAA" "GAC" "GAG" "GAT"
#> [37] "GCA" "GCC" "GCG" "GCT" "GGA" "GGC" "GGG" "GGT" "GTA" "GTC" "GTG" "GTT"
#> [49] "TAA" "TAC" "TAG" "TAT" "TCA" "TCC" "TCG" "TCT" "TGA" "TGC" "TGG" "TGT"
#> [61] "TTA" "TTC" "TTG" "TTT"
#> 
#> $kmer_value
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

## Notes and next steps

- The package focuses on feature generation. Model training (XGBoost) is
  expected to be performed outside the package using standard tools.
- If you previously used MIC functions, consult the function reference
  to map older workflows to the new `faLearn` utilities.

## Contact

Report issues or feature requests on the GitHub repository:
<https://github.com/agerada/faLearn>
