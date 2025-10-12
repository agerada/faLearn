##' Example genomes dataset
##'
##' A synthetic dataset for examples and demos. The dataset is stored as a
##' single object, \code{example_genomes}, which is a list with the following
##' components:
##' \describe{
##'   \item{genomes}{Named list of length 100; each element is a
##'     \code{Biostrings::DNAStringSet} holding contigs for that synthetic genome.}
##'   \item{annotations}{A two-column logical matrix (columns \code{gene1} and \code{gene2}) indicating presence/absence of the two inserted genes. Row order matches \code{example_genomes$genomes}.}
##'   \item{labels}{Numeric vector of MIC values on the log2 scale (log2(mg/L)), in the same order as \code{names(example_genomes$genomes)}.}
##' }
##'
##' These data are synthetic and intended only for package examples and testing.
##'
##' @docType data
##' @keywords datasets
##' @name example_genomes
##' @usage data(example_genomes)
##' @examples
##' data("example_genomes", package = "faLearn")
##' names(example_genomes$genomes)[1]
##' example_genomes$labels[1]
"example_genomes"
