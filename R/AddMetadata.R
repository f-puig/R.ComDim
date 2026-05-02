#' AddMetadata
#'
#' Adds or overwrites batch and/or metadata information for one or more blocks
#' of an existing MultiBlock object.
#' @param MB A MultiBlock object.
#' @param block The name of the block to modify (character). If \code{NULL} and
#'   the MultiBlock contains a single block, that block is used automatically.
#'   To modify several blocks at once, pass a character vector of block names.
#' @param batches A numeric vector of batch labels, one per sample in the
#'   target block. Replaces any existing Batch entry for that block.
#' @param metadata A data.frame with one row per sample in the target block.
#'   Replaces any existing Metadata entry for that block.
#' @return The updated MultiBlock.
#' @seealso \code{\link{MultiBlock}}, \code{\link{FilterSamplesMultiBlock}}
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50) # 10 samples, 50 variables
#' b2 <- matrix(rnorm(800), 10, 80) # 10 samples, 80 variables
#' batch_b1 <- rep(1, 10)
#' meta_b1 <- data.frame(condition = rep(c("A", "B"), 5))
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' # Add batch information to block 'b1':
#' mb <- AddMetadata(mb, block = "b1", batches = batch_b1)
#' # Add (or overwrite) metadata for block 'b1':
#' mb <- AddMetadata(mb, block = "b1", metadata = meta_b1)
#' @export
AddMetadata <- function(MB, block = NULL, batches = NULL, metadata = NULL) {
  if (!is(MB, "MultiBlock")) {
    stop("'MB' must be a MultiBlock object.")
  }

  if (is.null(block)) {
    if (length(MB@Data) == 1) {
      block <- names(MB@Data)
    } else {
      stop("'block' must be specified when the MultiBlock contains more than one block.")
    }
  }

  if (!all(block %in% names(MB@Data))) {
    missing_blocks <- block[!block %in% names(MB@Data)]
    stop(sprintf(
      "Block(s) not found in MultiBlock: %s",
      paste(missing_blocks, collapse = ", ")
    ))
  }

  for (blk in block) {
    n <- nrow(MB@Data[[blk]])

    if (!is.null(batches)) {
      if (!is.numeric(batches)) {
        stop("'batches' must be a numeric vector.")
      }
      if (length(batches) != n) {
        stop(sprintf(
          "'batches' must have %d elements (samples in block '%s').", n, blk
        ))
      }
      MB@Batch[[blk]] <- batches
    }

    if (!is.null(metadata)) {
      if (!is.data.frame(metadata)) {
        stop("'metadata' must be a data.frame.")
      }
      if (nrow(metadata) != n) {
        stop(sprintf(
          "'metadata' must have %d rows (samples in block '%s').", n, blk
        ))
      }
      MB@Metadata[[blk]] <- metadata
    }
  }

  return(MB)
}
