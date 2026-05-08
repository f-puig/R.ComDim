#' MultiAssayExperiment2MultiBlock
#'
#' Converts a MultiAssayExperiment into a MultiBlock. Samples are first intersected across all
#' experiments so that only common samples are retained.
#' @param se A \code{MultiAssayExperiment} object.
#' @param colData_samplenames Character string giving the name of the column in \code{colData(se)}
#'   that holds the sample names used for matching. When \code{NULL}, \code{colData} metadata is
#'   not stored in the resulting MultiBlock.
#' @param Batch Character string giving the name of the column in \code{colData(se)} to use as
#'   the Batch variable. Requires \code{colData_samplenames} to be set; an error is raised if
#'   \code{Batch} is supplied but \code{colData_samplenames} is \code{NULL} or the column is not
#'   found in \code{colData}.
#' @details
#'   Columns (samples) are intersected across all experiments via \code{intersectColumns} before
#'   conversion. Each experiment in the \code{MultiAssayExperiment} becomes one block in the
#'   MultiBlock (rows = samples, columns = features). Sample names are taken from the
#'   \code{colData} row names of \code{se}. Metadata from \code{colData} is stored only for the
#'   first block; subsequent blocks share the same sample order but carry no additional metadata.
#'   If \code{Batch} is specified, the corresponding column is extracted from \code{colData},
#'   removed from the metadata, and stored as the \code{Batch} slot of the MultiBlock.
#' @return A \code{MultiBlock} object with one block per experiment in \code{se}.
#' @seealso \code{\link{MultiBlock}}, \code{\link{MultiBlock2MultiAssayExperiment}}
#' @examples
#' \donttest{
#' if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
#' library(MultiAssayExperiment)
#' mae <- MultiAssayExperiment(
#'   experiments = ExperimentList(
#'     block1 = matrix(rnorm(50), nrow = 5, dimnames = list(paste0("s", 1:5), paste0("v", 1:10))),
#'     block2 = matrix(rnorm(30), nrow = 5, dimnames = list(paste0("s", 1:5), paste0("w", 1:6)))
#'   )
#' )
#' mb <- MultiAssayExperiment2MultiBlock(mae)
#' }
#' }
#' @export
MultiAssayExperiment2MultiBlock <- function(se,
                                            colData_samplenames = NULL, # Character string with the column metadata.
                                            Batch = NULL) { # Batch is the name of the column that will be considered as Batch

  # Here the metadata of the MultiAssayExperiment2MultiBlock is stored in the first MultiBlock.

  se <- intersectColumns(se[, , names(experiments(se))])

  simplelist <- assays(se)

  samplemap <- sampleMap(se)

  for (i in seq_along(simplelist@listData)) {
    Data <- list(t(as.matrix(simplelist@listData[[i]])))
    names(Data) <- names(simplelist@listData)[i]
    Variables <- list(rownames(simplelist@listData[[i]]))
    names(Variables) <- names(Data)
    eq.names <- samplemap@listData$primary[match(rownames(Data[[1]]), samplemap@listData$colname)]
    Data[[1]] <- Data[[1]][match(se@colData@rownames, eq.names), ]
    rownames(Data[[1]]) <- NULL
    colnames(Data[[1]]) <- NULL

    if (i == 1) {
      Samples <- se@colData@rownames
      if (is.null(colData_samplenames)) {
        Metadata <- NULL
      } else {
        Metadata <- as.data.frame(colData(se))
      }
      if (is.null(Metadata) && is.null(Batch)) {
        MB <- MultiBlock(
          Samples = Samples,
          Data = Data,
          Variables = Variables
        )
      } else if (is.null(Metadata) && !is.null(Batch)) {
        stop("There is no metadata in the MultiAssayExperiment to extract the Batch")
      } else if (!is.null(Metadata) && is.null(Batch)) {
        if (!(colData_samplenames %in% colnames(Metadata))) {
          stop("The string in 'colData_samplenames' is not a column name of colData.")
        }
        Metadata <- Metadata[match(
          Metadata[[colData_samplenames]],
          se@colData@rownames
        ), ]
        if (any(is.na(match(
          Metadata[[colData_samplenames]],
          se@colData@rownames
        )))) {
          Metadata <- Metadata[!is.na(match(
            Metadata[[colData_samplenames]],
            se@colData@rownames
          )), ]
        }
        Metadata <- list(x = Metadata)
        names(Metadata) <- names(Data)

        MB <- MultiBlock(
          Samples = Samples,
          Data = Data,
          Variables = Variables,
          Metadata = Metadata
        )
      } else {
        if (Batch %in% colnames(Metadata)) {
          batch_col <- Batch
          Batch <- Metadata[[batch_col]]
          noBatch <- setdiff(colnames(Metadata), batch_col)
          Metadata <- Metadata[, noBatch, drop = FALSE]
          Metadata <- Metadata[match(
            Metadata[[colData_samplenames]],
            se@colData@rownames
          ), ]
          if (any(is.na(match(
            Metadata[[colData_samplenames]],
            se@colData@rownames
          )))) {
            Metadata <- Metadata[!is.na(match(
              Metadata[[colData_samplenames]],
              se@colData@rownames
            )), ]
          }
          Metadata <- list(x = Metadata)
          names(Metadata) <- names(Data)
          Batch <- setNames(list(Batch), names(Data))
          MB <- MultiBlock(
            Samples = Samples,
            Data = Data,
            Variables = Variables,
            Batch = Batch,
            Metadata = Metadata
          )
        } else {
          stop("The string in batch was not found in the metadata.")
        }
      }
    } else {
      MB@Data <- c(MB@Data, Data)
      MB@Variables <- c(MB@Variables, Variables)
    }
  }

  return(MB)
}
