#' R.ComDim: ComDim Analysis in R
#'
#' Provides functions to run ComDim analysis. ComDim is an unsupervised
#' multi-block method that simultaneously considers multiple data tables to
#' find latent components common to all tables as well as those specific to
#' each table.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats var sd median rnorm qnorm setNames
#' @importFrom methods hasArg
#' @importFrom utils globalVariables
NULL

# Suppress NOTEs about Bioconductor functions used in optional conversion helpers.
# These packages (MultiAssayExperiment, SummarizedExperiment) are in Suggests.
utils::globalVariables(c(
  "ExperimentList", "MultiAssayExperiment",
  "assays", "colData", "experiments", "intersectColumns", "sampleMap"
))
