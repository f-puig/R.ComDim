#' ProcessMultiBlock
#'
#' Apply a custom function to transform a MultiBlock and/or select variables or blocks.
#' When multiple operations are supplied, the order of execution is:
#' \code{vars} subsetting, then \code{FUN}, then \code{FUN.SelectVars}, then \code{FUN.SelectBlocks}.
#' @param MB A MultiBlock object.
#' @param blocks The blocks to process. A vector of integers or block names. When omitted, all
#'   blocks are processed.
#' @param vars The variables to keep. A list with the same length as \code{blocks}. Each element
#'   contains the indices (integer) or names (character) of the variables to retain in that block.
#'   Subsetting is applied before any of the \code{FUN} arguments.
#' @param FUN A function applied to each selected block's data matrix (samples x variables). It
#'   receives the matrix as its sole argument and must return a matrix of the same dimensions.
#' @param FUN.SelectVars A function applied to each selected block's data matrix to determine which
#'   variables to keep. It receives the matrix as its sole argument and must return a logical vector
#'   of length equal to \code{ncol} of the matrix (\code{TRUE} = keep).
#' @param FUN.SelectBlocks A function applied to each selected block's data matrix to determine
#'   whether the block should be retained. It receives the matrix as its sole argument and must
#'   return a single \code{TRUE} or \code{FALSE}.
#' @return The processed MultiBlock object, with data matrices transformed and/or blocks/variables
#'   removed according to the supplied arguments. Blocks not listed in \code{blocks} are left
#'   unchanged.
#' @seealso \code{\link{MultiBlock}}, \code{\link{MultiBlock2Matrix}}
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' # Normalize each block to 0-100 range
#' BY100 <- function(x) 100 * (x - min(x)) / (max(x) - min(x))
#' mb <- ProcessMultiBlock(mb, FUN = BY100)
#' # Keep only variables with non-zero variance
#' mb <- ProcessMultiBlock(mb, FUN.SelectVars = function(x) apply(x, 2, var) > 0)
#' # Remove blocks where all values are below 0.5
#' mb <- ProcessMultiBlock(mb, FUN.SelectBlocks = function(x) max(x) >= 0.5)
#' @export
ProcessMultiBlock <- function(MB = MB,
                              blocks = NULL,
                              vars = NULL,
                              FUN = NULL,
                              FUN.SelectVars = NULL,
                              FUN.SelectBlocks = NULL) {
  if (is.null(FUN) && is.null(FUN.SelectBlocks) &&
    is.null(FUN.SelectVars) && is.null(blocks) &&
    is.null(vars)) {
    warning("'FUN','FUN.SelectBlocks' and 'FUN.SelectVars' arguments were missing. 'MB' was not processed.")
    return(MB)
  }

  if (!inherits(MB, "MultiBlock")) {
    stop("'MB' is not a MultiBlock object.")
  }

  if (!(is.null(blocks))) {
    if (is.character(blocks)) {
      block.ids <- match(blocks, names(MB@Data))
      if (any(is.na(block.ids))) {
        warning(
          "Some block names were not found in the MultiBlock: %s",
          paste0(blocks[which(is.na(block.ids))], sep = ", ")
        )
        block.ids <- block.ids[!is.na(block.ids)]
      }
    } else if (is.numeric(blocks)) {
      block.ids <- intersect(seq_along(MB@Data), blocks)
      if (length(blocks) != length(block.ids)) {
        warning(
          "Some block names were not found in the MultiBlock: %s",
          paste0(setdiff(seq_along(MB@Data), blocks), sep = ", ")
        )
      }
    }
  } else {
    block.ids <- seq(1, length(MB@Data))
  }
  block.fixed <- setdiff(seq(1, length(MB@Data)), block.ids)

  if (!is.null(vars)) {
    if (!is.list(vars)) {
      stop("'variables' must be a list object.")
    } else if (length(vars) != length(block.ids)) {
      stop("'variables' should have the same length as 'blocks'.")
    }

    for (j in seq_along(vars)) {
      i <- block.ids[j]
      if (is.character(vars[[j]])) {
        xx <- match(vars[[j]], MB@Variables[[i]])
        if (any(is.na(xx))) {
          warning(sprintf(
            "Some variable names were not found in block %s: %s", names(MB@Data)[i],
            paste0(vars[[j]][which(is.na(xx))], collapse = ", ")
          ))
          xx <- xx[!is.na(xx)]
        }
      } else if (is.numeric(vars[[j]])) {
        xx <- intersect(seq_along(MB@Variables[[i]]), vars[[j]])
        if (length(xx) != length(vars[[j]])) {
          warning(sprintf(
            "Some variable indices were out of range for block %s: %s", names(MB@Data)[i],
            paste0(setdiff(vars[[j]], seq_along(MB@Variables[[i]])), collapse = ", ")
          ))
        }
      }
      MB@Data[[i]]      <- MB@Data[[i]][, xx, drop = FALSE]
      MB@Variables[[i]] <- MB@Variables[[i]][xx]
    }
  }

  if (!is.null(FUN)) {
    for (i in block.ids) {
      MB@Data[[i]] <- FUN(MB@Data[[i]])
    }
  }

  if (!is.null(FUN.SelectVars)) {
    for (i in block.ids) {
      pos <- which(FUN.SelectVars(MB@Data[[i]]))

      if (length(pos) != 0) {
        MB@Data[[i]] <- MB@Data[[i]][, pos]
        MB@Variables[[i]] <- MB@Variables[[i]][pos]
      }
    }
  }

  if (!is.null(FUN.SelectBlocks)) {
    to_remove <- integer(0)
    for (i in block.ids) {
      if (!FUN.SelectBlocks(MB@Data[[i]])) {
        to_remove <- c(to_remove, i)
      }
    }
    if (length(to_remove) > 0) {
      to_remove <- sort(to_remove, decreasing = TRUE)
      for (i in to_remove) {
        MB@Data[[i]] <- NULL
        MB@Variables[[i]] <- NULL
      }
    }
  }

  return(MB)
}
