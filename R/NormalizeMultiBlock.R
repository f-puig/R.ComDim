# Pareto, mean-center, auto-scale, geometric mean, ranknorm

#' NormalizeMultiBlock
#'
#' Normalize all blocks from a MultiBlock object. The ranknorm transform is
#' based on that from the RNOmni package (doi:10.1111/biom.13214).
#' @param MB The MultiBlock object.
#' @param blocks Blocks to normalize. A vector of integers or block names
#'   (optional; all blocks are processed when omitted).
#' @param method Normalization method. One of: \code{'none'} (no
#'   transformation beyond optional Inf-to-NA conversion), \code{'auto'}
#'   (auto-scaling), \code{'mean'} (mean-centering), \code{'pareto'}
#'   (mean-centered and divided by the square root of the SD), \code{'norm'}
#'   (mean-centered and divided by its Frobenius norm), \code{'geometric'}
#'   (geometric mean scaling), \code{'ranknorm'} (rank-based inverse normal
#'   transform). Default: \code{'norm'}.
#' @param infinite.as.NA If \code{TRUE}, Infinite values are converted to NA
#'   before normalization. The presence of NA or Inf values will cause ComDim
#'   to stop; use \code{\link{NAInfRemoveMultiBlock}} to handle them.
#' @param constant For \code{'geometric'}: value added to each element before
#'   computing the geometric mean (useful when data contain zeros). Default: 0.
#' @param offset For \code{'ranknorm'}: offset in the Blom-type transform.
#'   Default: 3/8 (Blom transform).
#' @param showWarning If \code{TRUE}, warn when all-NA variables are
#'   discarded. Default: \code{TRUE}.
#' @return The normalized MultiBlock object.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' b2[c(2, 3, 5), c(1, 2, 3)] <- NA
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' mb <- NormalizeMultiBlock(mb, method = "auto")
#' @export
NormalizeMultiBlock <- function(MB,
                                blocks = NULL,
                                method = c(
                                  "none", "auto", "mean", "pareto",
                                  "norm", "geometric", "ranknorm"
                                )[5],
                                infinite.as.NA = FALSE,
                                constant = 0,
                                offset = 3 / 8,
                                showWarning = TRUE) {
  if (is.null(method) || length(method) > 1) {
    stop("A method should be provided. Try 'none', 'auto', 'mean', 'pareto', 'norm', 'geometric', or 'ranknorm'.")
  }

  if (!("MultiBlock" %in% class(MB))) {
    stop("'MB' is not a MultiBlock object.")
  }

  method <- tolower(method)
  if (!(method %in% c("none", "auto", "mean", "pareto", "norm", "geometric", "ranknorm"))) {
    stop("The provided method is not in the list of available methods. Try 'none', 'auto', 'mean', 'norm', 'pareto', 'geometric', or 'ranknorm'.")
  }

  if (!is.null(blocks)) {
    if (is.character(blocks)) {
      block.ids <- match(blocks, names(MB@Data))
      if (any(is.na(block.ids))) {
        warning(sprintf(
          "Some block names were not found in the MultiBlock: %s",
          paste(blocks[is.na(block.ids)], collapse = ", ")
        ))
        block.ids <- block.ids[!is.na(block.ids)]
      }
    } else if (is.numeric(blocks)) {
      block.ids <- intersect(seq_along(MB@Data), blocks)
      if (length(blocks) != length(block.ids)) {
        warning(sprintf(
          "Some block indices were not found in the MultiBlock: %s",
          paste(setdiff(blocks, seq_along(MB@Data)), collapse = ", ")
        ))
      }
    }
  } else {
    block.ids <- seq_along(MB@Data)
  }

  if (length(block.ids) == 0) {
    warning("0 blocks were selected. Data was not modified.")
    return(MB)
  }

  to_ignore <- 0
  to_delete <- integer(0)

  for (i in block.ids) {
    if (any(is.complex(MB@Data[[i]]))) {
      stop("The multi-block contains complex values that cannot be processed with this function.")
    }

    if (any(is.logical(MB@Data[[i]]))) {
      stop("The multi-block contains logical values that cannot be processed with this function.")
    }

    if (any(is.infinite(MB@Data[[i]]))) {
      if (infinite.as.NA) {
        MB@Data[[i]][is.infinite(MB@Data[[i]])] <- NA
      } else {
        warning(sprintf(
          "Block '%s' contains Infinite values. Convert them to NA with infinite.as.NA = TRUE or using NAInfRemoveMultiBlock().",
          names(MB@Data)[i]
        ))
      }
    }

    # Discard variables that are entirely NA
    numna <- apply(MB@Data[[i]], 2, function(x) sum(is.na(x)))
    colna <- which(numna >= nrow(MB@Data[[i]]))
    if (length(colna) > 0) {
      keep_idx <- setdiff(seq_len(ncol(MB@Data[[i]])), colna)
      MB@Data[[i]] <- MB@Data[[i]][, keep_idx, drop = FALSE]
      MB@Variables[[i]] <- MB@Variables[[i]][keep_idx]
      to_ignore <- to_ignore + length(colna)
    }
    if (ncol(MB@Data[[i]]) == 0) {
      to_delete <- c(to_delete, i)
    }
  }

  # Remove blocks that became empty
  if (length(to_delete) > 0) {
    to_delete <- sort(to_delete, decreasing = TRUE)
    names_to_delete <- names(MB@Data)[to_delete]
    old_names <- names(MB@Data)
    for (i in to_delete) {
      nm <- names(MB@Data)[i]
      if (nm %in% names(MB@Batch)) MB@Batch[[nm]] <- NULL
      if (nm %in% names(MB@Metadata)) MB@Metadata[[nm]] <- NULL
      MB@Data[[i]] <- NULL
      MB@Variables[[i]] <- NULL
    }
    kept_names <- setdiff(old_names[block.ids], names_to_delete)
    if (length(kept_names) == 0) {
      stop("The selected blocks contained only NA values. They cannot be used in ComDim.")
    }
    block.ids <- match(kept_names, names(MB@Data))
  }

  for (i in block.ids) {
    if (method == "none") {
      next # Can be used solely to convert Inf to NA via infinite.as.NA.
    } else if (method == "auto") {
      MB@Data[[i]] <- scale(MB@Data[[i]], center = TRUE, scale = TRUE)
    } else if (method == "mean") {
      MB@Data[[i]] <- scale(MB@Data[[i]], center = TRUE, scale = FALSE)
    } else if (method == "pareto") {
      MB@Data[[i]] <- scale(MB@Data[[i]],
        center = TRUE,
        scale = apply(
          MB@Data[[i]], 2,
          function(x) sqrt(sd(x, na.rm = TRUE))
        )
      )
    } else if (method == "norm") {
      Xmean <- scale(MB@Data[[i]], center = TRUE, scale = FALSE)
      Norm_X <- sqrt(sum(Xmean * Xmean, na.rm = TRUE))
      MB@Data[[i]] <- Xmean / Norm_X
    } else if (method == "geometric") {
      for (nc in seq_len(ncol(MB@Data[[i]]))) {
        notna <- which(!is.na(MB@Data[[i]][, nc]))
        if (length(notna) == 0) next
        vals <- MB@Data[[i]][notna, nc] + constant
        if (any(vals <= 0)) {
          stop(sprintf(
            "Block '%s', column %d: values must be positive for the geometric mean. Use the 'constant' argument to shift them.",
            names(MB@Data)[i], nc
          ))
        }
        MB@Data[[i]][notna, nc] <- vals / prod(vals)^(1 / length(vals))
      }
    } else if (method == "ranknorm") {
      for (nc in seq_len(ncol(MB@Data[[i]]))) {
        notna <- which(!is.na(MB@Data[[i]][, nc]))
        n <- length(notna)
        r <- rank(MB@Data[[i]][notna, nc])
        MB@Data[[i]][notna, nc] <- qnorm((r - offset) / (n - 2 * offset + 1))
      }
    }
  }

  if (to_ignore != 0 && showWarning) {
    warning(sprintf("%d variable(s) contained only NAs and were deleted.", to_ignore))
  }

  return(MB)
}
