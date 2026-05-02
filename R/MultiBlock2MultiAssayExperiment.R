#' MultiBlock2MultiAssayExperiment
#'
#' Converts a MultiBlock into a MultiAssayExperiment. Each block becomes one experiment; sample
#' names, batch information, and metadata are carried over into the \code{colData} of the result.
#' @param MB A \code{MultiBlock} object.
#' @param MSEmetadata An optional list of unstructured metadata describing the overall content of
#'   the MultiAssayExperiment (stored in its \code{metadata} slot). Pass \code{NULL} (default) to
#'   omit experiment-level metadata.
#' @details
#'   Each block in \code{MB} is transposed (features x samples) and stored as a named matrix in
#'   the \code{ExperimentList}. Row names are set to the variable names and column names to the
#'   sample names of the MultiBlock. The \code{colData} is constructed from \code{MB@@Samples};
#'   any \code{Metadata} and \code{Batch} information present in the MultiBlock is appended as
#'   additional columns. A \code{sampleMap} is generated mapping every sample to every experiment
#'   using the same primary and column names.
#' @return A \code{MultiAssayExperiment} object with one experiment per block in \code{MB}.
#' @seealso \code{\link{MultiBlock}}, \code{\link{MultiAssayExperiment2MultiBlock}}
#' @examples
#' \dontrun{
#' b1 <- matrix(rnorm(50), 5, 10, dimnames = list(paste0("s", 1:5), paste0("v", 1:10)))
#' b2 <- matrix(rnorm(30), 5,  6, dimnames = list(paste0("s", 1:5), paste0("w", 1:6)))
#' mb <- MultiBlock(Data = list(block1 = b1, block2 = b2))
#' mae <- MultiBlock2MultiAssayExperiment(mb)
#' mae <- MultiBlock2MultiAssayExperiment(mb, MSEmetadata = list(study = "example"))
#' }
#' @export
MultiBlock2MultiAssayExperiment <- function(MB,
                                            MSEmetadata = NULL) { # List of unstructured metadata describing the overall content of the object.

  if (!inherits(MB, "MultiBlock")) {
    stop("'MB' is  not a MultiBlock.")
  }

  ## Create ExperimentList
  exp <- list()
  for (i in seq_along(MB@Data)) {
    exp[[names(MB@Data)[i]]] <- t(MB@Data[[i]])
    colnames(exp[[names(MB@Data)[i]]]) <- MB@Samples
    rownames(exp[[names(MB@Data)[i]]]) <- MB@Variables[[names(MB@Data)[i]]]
  }
  exp <- ExperimentList(exp)

  ## Create colData
  colDat <- matrix(MB@Samples, nrow = length(MB@Samples), ncol = 1)
  for (i in seq_along(MB@Data)) {
    if (names(MB@Data)[i] %in% names(MB@Metadata)) {
      colDat <- cbind.data.frame(colDat, MB@Metadata[[names(MB@Metadata)[i]]])
    }
    if (names(MB@Data)[i] %in% names(MB@Batch)) {
      colDat <- cbind.data.frame(colDat, MB@Batch[[names(MB@Data)[i]]])
    }
  }
  rownames(colDat) <- MB@Samples

  # SampleMap
  for (i in seq_along(MB@Data)) {
    if (i == 1) {
      SEsampleMap <- data.frame(
        assay = rep(names(MB@Data)[i], length(MB@Samples)),
        primary = MB@Samples,
        colname = MB@Samples
      )
    } else {
      SEsampleMap <- rbind.data.frame(
        SEsampleMap,
        data.frame(
          assay = rep(
            names(MB@Data)[i],
            length(MB@Samples)
          ),
          primary = MB@Samples,
          colname = MB@Samples
        )
      )
    }
  }


  if (is.null(MSEmetadata)) {
    se <- MultiAssayExperiment(
      experiments = exp,
      colData = colDat,
      sampleMap = SEsampleMap
    )
  } else {
    se <- MultiAssayExperiment(
      experiments = exp,
      colData = colDat,
      sampleMap = SEsampleMap,
      metadata = MSEmetadata
    )
  }

  return(se)
}
