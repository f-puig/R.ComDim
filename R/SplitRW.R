#' SplitRW
#'
#' Splits a multi-block into a replicate-wise (RW) structure by expanding each block along its
#' batch dimension. Each batch within each original block becomes a separate block in the output,
#' enabling replicate-wise ComDim analysis.
#' @param MB A \code{MultiBlock} object built with \code{\link{MultiBlock}()}. Must contain
#'   \code{Batch} information; blocks without batch labels are treated as a single batch.
#' @param checkSampleCorrespondence Logical. If \code{FALSE} (default), identical sample order
#'   and count are assumed across all batches and only the minimum batch size is used. If
#'   \code{TRUE}, sample names are intersected across all batches and only common samples are
#'   retained, which is safer when batches have different sizes or orderings.
#' @param batchNormalisation Logical. If \code{TRUE} (default), each replicate block is divided
#'   by its Frobenius norm and by the square root of the number of replicates in its original
#'   block, to prevent blocks with more replicates from dominating the analysis.
#' @param showSampleCorrespondence Logical. If \code{TRUE} (default), a matrix showing the
#'   sample names assigned to each replicate block is printed to the console.
#' @details
#'   Output block names follow the convention \code{<original_block>} when the original block has
#'   only one batch, or \code{<original_block>_<batch_label>} when it has multiple batches.
#'   The \code{Metadata} slot of each source block is also split and carried over to the
#'   corresponding replicate blocks. If the MultiBlock has no \code{Batch} information at all,
#'   the original object is returned unchanged with a warning.
#' @return A \code{MultiBlock} object in which each block corresponds to one batch of one
#'   original data block (a replicate-wise structure ready for \code{\link{ComDim_PCA}} or
#'   similar).
#' @seealso \code{\link{MultiBlock}}, \code{\link{ComDim_PCA}}, \code{\link{SelectFeaturesRW}}
#' @examples
#' b1 <- matrix(rnorm(1500), 30, 50)
#' b2 <- matrix(rnorm(2400), 30, 80)
#' batch_b <- c(rep(1, 10), rep(2, 10), rep(3, 10))
#' # Generate the multi-block (mb) with 3 batches of 10 samples each
#' mb <- MultiBlock(
#'   Data = list(b1 = b1, b2 = b2),
#'   Batch = list(b1 = batch_b, b2 = batch_b),
#'   ignore.names = TRUE
#' )
#' rw <- SplitRW(mb)
#' @export
SplitRW <- function(MB = MB, checkSampleCorrespondence = FALSE,
                    batchNormalisation = TRUE, showSampleCorrespondence = TRUE) {
  newMB <- list()

  if (!inherits(MB, "MultiBlock")) {
    stop("'MB' is  not a MultiBlock.")
  }

  if (!is.list(MB@Samples)) {
    samples <- MB@Samples
    MB@Samples <- list()
    xx <- names(MB@Data)
    for (i in seq_along(xx)) {
      MB@Samples[[xx[i]]] <- samples
    }
    rm(samples)
  }

  # Check that the batch information is correctly provided for all the blocks:
  give_error <- 0

  batch_names <- list()
  samples_per_batch <- list()

  if (length(MB@Batch) != 0) {
    for (i in names(MB@Data)) { # Check the batch information is consistent across RWs
      if (i %in% names(MB@Batch)) {
        batch_names[[i]] <- sort(unique(MB@Batch[[i]]))
        sampled <- as.vector(NULL)
        for (j in seq_along(batch_names[[i]])) {
          sampled[j] <- length(which(MB@Batch[[i]] == batch_names[[i]][j]))
        }
        samples_per_batch[[i]] <- sampled
        if (length(unique(sampled)) > 1) {
          warning(sprintf("The replicate blocks in %s have a different number of samples.", names(MB@Data)[i]))
        }
      } else {
        batch_names[[i]] <- "Batch 1"
        samples_per_batch[[i]] <- MB@Samples[[i]]
      }
    }
    if (length(unique(unlist(samples_per_batch))) > 1 && checkSampleCorrespondence == FALSE) {
      print("Information is missing regarding the splitting.")
      print("Using checkSampleCorrespondence as TRUE is recommended.")
      stop("The data cannot be split into replicate blocks.")
    }

    for (i in names(MB@Batch)) {
      if (checkSampleCorrespondence == TRUE) {
        for (j in seq_along(batch_names[[i]])) {
          batch_position <- which(MB@Batch[[i]] == batch_names[[i]][j])
          if (i == names(MB@Batch)[1] && j == 1) {
            replicate_names <- MB@Samples[[i]][batch_position]
          } else {
            replicate_names <- intersect(replicate_names, MB@Samples[[i]][batch_position])
          }
        }
      } else {
        replicate_names <- seq_len(min(unlist(samples_per_batch))) # number of samples in the smallest batch
      }
    }

    if (checkSampleCorrespondence == TRUE && length(replicate_names) == 0) {
      print("There are 0 samples in common across the replicate blocks")
      give_error <- 1
    }

    if (give_error) {
      stop("The data cannot be split into replicate blocks.")
    }
  } else {
    warning("The MultiBlock does not contain 'Batch' information. It cannot be split into replicate blocks.")
    return(MB)
  }


  ## PROCEED WITH THE RW SPLITTING.
  # df_SampleNames is a table printed in the console (go to the end of the script for more info)
  if (showSampleCorrespondence) {
    df_SampleNames <- matrix(, nrow = length(replicate_names), ncol = 0)
  }

  k <- 1
  for (i in names(MB@Batch)) {
    for (j in seq_along(batch_names[[i]])) {
      batch_position <- which(MB@Batch[[i]] == batch_names[[i]][j])
      replicate_position <- as.vector(NULL)
      if (checkSampleCorrespondence == TRUE) {
        for (namei in replicate_names) {
          replicate_position <- c(
            replicate_position,
            intersect(which(MB@Samples[[i]] == namei), batch_position)
          )
        }
      } else {
        replicate_position <- batch_position
      }
      # Keep only the common samples across blocks
      if (length(replicate_position) != length(replicate_names)) {
        print(sprintf("There are sample duplicates in batch %s from block %s.", as.character(batch_names[[i]][j]), as.character(i)))
        print("Duplicate samples should be removed.")
        stop("Existence of duplicate samples within one or more batches.")
      }

      sorted <- order(MB@Samples[[i]][replicate_position])

      df_SampleNames <- cbind(df_SampleNames, MB@Samples[[i]][replicate_position[sorted]])

      growingMB <- MultiBlock(
        Samples = replicate_names,
        Data = list(s1 = MB@Data[[i]][replicate_position[sorted], ]),
        Variables = list(s1 = MB@Variables[[i]]),
        Batch = list(s1 = rep(batch_names[[i]][j], length(sorted)))
      )

      if (length(batch_names[[i]]) == 1) {
        blockNames(growingMB) <- i
      } else {
        blockNames(growingMB) <- paste(i, as.character(batch_names[[i]][j]), sep = "_")
      }

      if (batchNormalisation) { # Divide each RW block by its norm and the sqrt(number_of_replicates) to avoid overrepresentation of the data in the multi-set.
        normed <- growingMB@Data[[1]] / norm(growingMB@Data[[1]], type = "F")
        growingMB@Data[[1]] <- normed / sqrt(length(batch_names[[i]]))
      }

      # if(length(names(MB[[i]])) > 5){ # In case there are additional fields
      if (length(MB@Metadata[[i]]) != 0) { # In case there is metadata

        growingMB@Metadata[[blockNames(growingMB)]] <- MB@Metadata[[i]][replicate_position[sorted], ]
      }
      k <- k + 1

      if (j == 1) {
        mbb <- growingMB
      } else {
        mbb@Data <- c(mbb@Data, growingMB@Data)
        mbb@Variables <- c(mbb@Variables, growingMB@Variables)
        mbb@Batch <- c(mbb@Batch, growingMB@Batch)
        mbb@Metadata <- c(mbb@Metadata, growingMB@Metadata)
      }
    }
    if (i == names(MB@Batch)[1]) {
      newMB <- mbb
    } else {
      newMB@Data <- c(newMB@Data, mbb@Data)
      newMB@Variables <- c(newMB@Variables, mbb@Variables)
      newMB@Batch <- c(newMB@Batch, mbb@Batch)
      newMB@Metadata <- c(newMB@Metadata, mbb@Metadata)
    }
  }


  ## SHOW A DATA-FRAME WITH ALL THE SAMPLE NAMES FOR ALL THE (REPLICATE) BLOCKS.
  ## It is wise to use it when checkSampleCorrespondence is set to FALSE.
  if (showSampleCorrespondence) {
    colnames(df_SampleNames) <- names(newMB@Data)
    print("The sample names are:")
    print(df_SampleNames)
  }

  return(newMB)
}
