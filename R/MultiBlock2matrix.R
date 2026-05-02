#' MultiBlock2Matrix
#'
#' Combines the blocks of a MultiBlock into a single matrix by column-binding the selected blocks.
#' Useful for computing summary statistics across the whole MultiBlock (e.g. maximum value).
#' @param MB A \code{MultiBlock} object.
#' @param blocks The blocks to combine. A vector of integers (block indices) or a vector of block
#'   names. When omitted, all blocks are included.
#' @param vars The variables to keep. A list of the same length as \code{blocks}. Each element
#'   contains the indices (integer) or names (character) of the variables to retain from the
#'   corresponding block. When omitted, all variables in each block are included.
#' @details
#'   Blocks are column-bound in the order given by \code{blocks}. Row names of the output matrix
#'   are the sample names of \code{MB}. Column names are the variable names; if the same variable
#'   name appears in more than one block, duplicates are disambiguated by appending
#'   \code{.<block_name>} to all occurrences of the repeated name, and a warning is issued.
#' @return A numeric matrix with rows corresponding to samples and columns corresponding to the
#'   variables of the selected blocks concatenated in order.
#' @seealso \code{\link{MultiBlock}}, \code{\link{ProcessMultiBlock}}
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' # Combine all blocks into a matrix and compute a summary statistic:
#' mat <- MultiBlock2Matrix(mb)
#' max(mat)
#' # Combine only the first block, keeping variables 1-10:
#' mat <- MultiBlock2Matrix(mb, blocks = 1, vars = list(1:10))
#' @export
MultiBlock2Matrix <- function(MB = MB,
                              blocks = NULL,
                              vars = NULL) {
  if (!("MultiBlock" %in% class(MB))) {
    stop("'MB' is not a MultiBlock object.")
  }

  if (!(is.null(blocks))) {
    if (is.character(blocks)) {
      block.ids <- match(blocks, names(MB@Data))
      if (all(is.na(block.ids))) {
        stop("None of the block names were found in the MultiBlock.")
      }
      if (any(is.na(block.ids))) {
        warning(
          "Some block names were not found in the MultiBlock: %s",
          paste0(blocks[which(is.na(block.ids))], sep = ", ")
        )
        block.ids <- block.ids[!is.na(block.ids)]
      }
    } else if (is.numeric(blocks)) {
      block.ids <- intersect(seq_along(MB@Data), blocks)
      if (length(block.ids) == 0) {
        stop("None of the block names were found in the MultiBlock.")
      }
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

  if (!is.null(vars)) {
    if (!is.list(vars)) {
      stop("'variables' must be a list object.")
    } else if (length(vars) != length(block.ids)) {
      stop("'variables' should have the same length as 'blocks'.")
    }

    for (j in seq_along(block.ids)) {
      i <- block.ids[j]
      if (is.character(vars[[j]])) {
        xx <- match(vars[[j]], MB@Variables[[i]])
        if (any(is.na(xx))) {
          warning(sprintf(
            "Some variable names were not found in block %s: %s", i,
            paste0(vars[[j]][which(is.na(xx))], sep = ", ")
          ))
          vars[[j]] <- xx[!is.na(xx)]
        } else {
          vars[[j]] <- xx
        }
      } else if (is.numeric(vars[[j]])) {
        xx <- intersect(seq_along(MB@Variables[[i]]), vars[[j]])
        if (length(xx) != length(vars[[j]])) {
          warning(sprintf(
            "Some variables were not found in the MultiBlock object %s: %s", i,
            paste0(setdiff(vars[[j]], seq_along(MB@Variables[[i]])), sep = ", ")
          ))
        }
        vars[[j]] <- xx
      }
    }
  } else {
    vars <- lapply(block.ids, function(i) seq_along(MB@Variables[[i]]))
  }

  cnames <- vector()
  rep_names <- vector()
  block_names <- vector()
  common <- FALSE

  for (j in seq_along(block.ids)) {
    i <- block.ids[j]
    if (j == 1) {
      mat <- MB@Data[[i]][, vars[[j]], drop = FALSE]
    } else {
      mat <- cbind(mat, MB@Data[[i]][, vars[[j]], drop = FALSE])
    }

    names_new <- MB@Variables[[i]]
    if (any(names_new %in% cnames)) {
      common <- TRUE
      rep_names <- append(rep_names, names_new[which(names_new %in% cnames)])
      names_new[which(names_new %in% cnames)] <-
        paste0(
          names_new[which(names_new %in% cnames)],
          sprintf(".%s", names(MB@Data)[i])
        )
    }
    block_names <- append(block_names, rep(names(MB@Data)[i], length(names_new)))
    cnames <- append(cnames, names_new)
  }

  # Now replace the name of the variables that had a replicated name.
  if (length(rep_names) != 0) {
    unique_rep_names <- as.character(unique(rep_names))
    cnames <- as.character(cnames)
    for (i in unique_rep_names) {
      pos <- which(cnames == i)
      cnames[pos] <- paste0(i, ".", block_names[match(i, cnames)])
    }
  }

  if (common) {
    warning("Some variable names were repeated across blocks.
  The repeated variable names were renamed.")
  }

  rownames(mat) <- MB@Samples
  colnames(mat) <- cnames

  return(mat)
}
