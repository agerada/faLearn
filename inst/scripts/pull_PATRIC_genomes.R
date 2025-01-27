#!/usr/bin/env Rscript

option_list <- list(
  optparse::make_option(c("-f", "--filter"),
              default = "MIC",
              help = "Which genomes to download,
              based on available meta-data. One of: all, MIC or disc.
              default [MIC]"),
  optparse::make_option(c("-d", "--database"),
              default = "ftp://ftp.bvbrc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt",
              help = "Database of genomes downloaded from PATRIC website:
              https://www.bv-brc.org/docs/quick_references/ftp.html
              default [ftp://ftp.bvbrc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt]"),
  optparse::make_option(c("-o", "--output_directory"),
              default = "data/genomes/patric/",
              help = "Folder to save downloaded genomes.
              default [data/genomes/patric]"),
  optparse::make_option(c("-n", "--n_genomes"),
              default = 0,
              help = "Number of genomes to try to download, where:
              0 = all,
              -1 = none,
              default [0]")
)

args <- optparse::parse_args(
  optparse::OptionParser(
    usage = "%script [options] taxonomic_name",
    option_list = option_list),
  positional_arguments = 1)

taxonomic_name <- args$args[[1]]
opt <- args$options

MIC::pull_PATRIC_genomes(database = opt$database,
                    taxonomic_name = taxonomic_name,
                    filter = opt$filter,
                    output_directory = opt$output_directory,
                    n_genomes = opt$n_genomes)
