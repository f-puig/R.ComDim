#' SummarizedExperiment2MultiBlock
#'
#' Converts a \code{SummarizedExperiment} into a MultiBlock. Each assay in the
#' \code{SummarizedExperiment} becomes one block (rows = samples, columns = features).
#' @param se A \code{SummarizedExperiment} object.
#' @param colData_samplenames Character string giving the name of the column in \code{colData(se)}
#'   that holds the sample names used for matching. When \code{NULL}, \code{colData} metadata is
#'   not stored in the resulting MultiBlock.
#' @param Batch Character string giving the name of the column in \code{colData(se)} to use as
#'   the Batch variable. When supplied, the column is extracted from \code{colData}, removed from
#'   the metadata, and stored in the \code{Batch} slot of the MultiBlock. Requires at least one
#'   of \code{colData_samplenames} or a named \code{colData}; an error is raised if the column
#'   is not found.
#' @details
#'   Each assay is transposed so that samples are in rows and features in columns. Sample order
#'   is aligned to \code{colData(se)} row names; samples not present in \code{colData} are
#'   removed. If an assay has no name, the block is labelled \code{"X"}. Metadata from
#'   \code{colData} (excluding the Batch column if specified) is stored only for the first block;
#'   subsequent assays are appended as additional blocks with no extra metadata.
#' @return A \code{MultiBlock} object with one block per assay in \code{se}.
#' @seealso \code{\link{MultiBlock}}, \code{\link{MultiAssayExperiment2MultiBlock}},
#'   \code{\link{MultiBlock2MultiAssayExperiment}}
#' @examples
#' \dontrun{
#' library(SummarizedExperiment)
#' nrows <- 20; ncols <- 6
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' rownames(counts) <- paste0("feature", seq_len(nrows))
#' colnames(counts) <- paste0("sample",  seq_len(ncols))
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' mb <- SummarizedExperiment2MultiBlock(se)
#' }
#' @export
SummarizedExperiment2MultiBlock <- function(se,
                                            colData_samplenames = NULL, # Character string with the column metadata.
                                            Batch = NULL) { # Batch is the name of the column that will be considered as Batch

  # Batch is the name inlistData with the replicate information.
  numassays <- length(assays(se))

  for (i in 1:numassays) {
    Data <- t(as.matrix(assays(se)[[i]]))
    Variables <- list(colnames(Data))
    Data <- list(Data)
    block_name <- names(assays(se))[[i]]
    if (length(block_name) > 1) {
      block_name <- "X"
    }
    names(Variables) <- names(Data) <- block_name

    Data[[block_name]] <- Data[[block_name]][match(rownames(Data[[block_name]]), se@colData@rownames), ]
    if (any(is.na(match(rownames(Data[[block_name]]), se@colData@rownames)))) {
      Data[[block_name]] <- Data[[block_name]][!is.na(match(rownames(Data[[block_name]]), se@colData@rownames)), ]
    }
    rownames(Data[[block_name]]) <- NULL
    colnames(Data[[block_name]]) <- NULL


    if (i == 1) { # First MB

      if (!is.null(Batch)) {
        if (length(Batch) > 1) {
          stop("'Batch' must be a character string of length 1.")
        }
        batch_col <- Batch
        Batch <- se@colData@listData[[batch_col]]
        noBatch <- setdiff(names(se@colData@listData), batch_col)
        if (is.null(colData_samplenames)) {
          MB <- MultiBlock(
            Samples = se@colData@rownames,
            Data = Data,
            Variables = Variables,
            Batch = setNames(list(Batch), block_name)
          )
        } else {
          Metadata <- as.data.frame(se@colData@listData)
          Metadata <- Metadata[, noBatch, drop = FALSE]
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
          Metadata <- list(Metadata)
          names(Metadata) <- block_name
          MB <- MultiBlock(
            Samples = as.character(Metadata[[block_name]][[colData_samplenames]]),
            Data = Data,
            Variables = Variables,
            Batch = setNames(list(Batch), block_name),
            Metadata = Metadata
          )
        }
      } else { # If Batch is NULL
        if (!is.null(colData_samplenames)) {
          Metadata <- as.data.frame(se@colData@listData)
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
          Metadata <- list(Metadata)
          names(Metadata) <- block_name
          MB <- MultiBlock(
            Samples = as.character(Metadata[[block_name]][[colData_samplenames]]),
            Data = Data,
            Variables = Variables,
            Metadata = Metadata
          )
        } else {
          MB <- MultiBlock(
            Samples = se@colData@rownames,
            Data = Data,
            Variables = Variables
          )
        }
      }
    } else { # Iterations i > 1
      MB@Data <- c(MB@Data, Data)
      MB@Variables <- c(MB@Variables, Variables)
    }
  }

  return(MB)
}
