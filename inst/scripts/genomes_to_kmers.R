#!/usr/bin/env Rscript

option_list <- list(
  optparse::make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Print extra output"),
  optparse::make_option(c("-c", "--cores"), type = "integer", default = 1,
              help =
                "When >1, indicates number of parallel cores to use [default 1]"),
  optparse::make_option(c("-k", "--kmers"), type = "integer", default = 3,
              help = "Kmer length [default 3]"),
  optparse::make_option(c("-n", "--n_genomes"), type = "integer",
              help = "Max number of genomes to process [default all]"),
  optparse::make_option(c("-p", "--split"), type = "integer", default = 1,
              help =
                "Split data into s chunks for parallel processing [default 1]"),
  optparse::make_option(c("-a", "--anchor"), action = "store_true", default = FALSE,
              help = "Include unobserved permutations"),
  optparse::make_option(c("-s", "--simplify"), action = "store_true", default = FALSE,
              help = "Store only the kmer counts (key: value -> value)"),
  optparse::make_option(c("-d", "--drop_n"), action = "store_true", default = TRUE,
              help = "Drop kmers that have N [default false]"),
  optparse::make_option(c("-f", "--formats"), type = "character",
              default = "libsvm",
              help =
                "Comma-separated list of formats to export, currently support:
              csv libsvm rds rdata parquet"),
  optparse::make_option(c("-i", "--integer_index"), action = "store_true",
              default = TRUE, help = "Force integer index"),
  optparse::make_option("--starting_index", type = "integer",
                        default = 1,
                        help = "Starting integer index, if integer_index TRUE"),
  optparse::make_option(c("-r", "--random_shuffle"), action = "store_true",
              default = TRUE,
              help =
                "Randomly shuffle genomes - useful for test train split with -p")
)

args <- optparse::parse_args(
  optparse::OptionParser(usage = "%script [options] input_dir output_dir",
               option_list = option_list),
  positional_arguments = 2)

opt <- args$options
dirs <- args$args
input_dir <- dirs[[1]]
output_dir <- dirs[[2]]

MIC::genomes_to_kmer_dataset(input_dir = input_dir,
                                output_dir = output_dir,
                                cores = opt$cores,
                                kmers = opt$kmers,
                                n_genomes = opt$n_genomes,
                                split = opt$split,
                                anchor = opt$anchor,
                                simplify = opt$simplify,
                                drop_n = opt$drop_n,
                                integer_index = opt$integer_index,
                                starting_index = opt$starting_index,
                                random_shuffle = opt$random_shuffle)
