#' NAInfRemoveMultiBlock
#'
#' Remove NA and Infinite values from a MultiBlock object in a single pass.
#' Variables whose combined count of NA and Infinite values meets or exceeds
#' \code{minfrac * nrow} are discarded first. Remaining Infinite values are
#' then replaced according to \code{inf.method}; remaining NA values are
#' imputed according to \code{na.method}.
#' @param MB The MultiBlock object.
#' @param blocks Blocks to process: a vector of integers or block names
#'   (optional; all blocks are processed when omitted).
#' @param minfrac Minimum fraction of valid (non-NA, finite) values required
#'   to retain a variable. Variables at or below this threshold are discarded.
#'   Default: 0.5.
#' @param na.method Imputation method for remaining NA values after the
#'   minfrac filter. One of: \code{'none'} (keep NAs), \code{'zero'},
#'   \code{'median'}, \code{'discard'} (drop columns still containing NAs),
#'   \code{'fixed.value'}, \code{'fixed.value.all'}, \code{'fixed.noise'},
#'   \code{'random.noise'}, \code{'QRILC'}. Default: \code{'random.noise'}.
#' @param inf.method Replacement method for remaining Infinite values after
#'   the minfrac filter. One of: \code{'none'} (keep Inf),
#'   \code{'fixed.noise'}, \code{'random.noise'}.
#'   Default: \code{'random.noise'}.
#' @param constant For \code{na.method = 'fixed.value'}: value to replace NAs.
#'   For \code{na.method = 'fixed.value.all'}: value added to every element
#'   (NAs are then set to this value). Default: 0.
#' @param factor.NA Noise factor used by \code{'fixed.noise'} and
#'   \code{'random.noise'} for both NA and Inf imputation. For NA and -Inf:
#'   the replacement mean equals \code{col_min * factor.NA}. For +Inf:
#'   \code{col_max / factor.NA}. Default: 0.5.
#' @param sd.noise Standard-deviation factor for \code{'random.noise'}:
#'   \code{SD = abs(mean * sd.noise)}. Default: 0.3.
#' @param seed.number Optional integer seed for reproducible
#'   \code{'random.noise'} results. The RNG state is restored on exit.
#' @param tune.sigma Tuning scalar for \code{na.method = 'QRILC'}. Default:
#'   1 (Gaussian complete-data assumption).
#' @param showWarning If \code{TRUE}, emit a warning when variables are
#'   discarded by the minfrac filter. Default: \code{TRUE}.
#' @return The processed MultiBlock object.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' b2[c(2, 3, 5), c(1, 2, 3)] <- NA
#' b2[c(1, 4), c(4, 5)] <- Inf
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' mb <- NAInfRemoveMultiBlock(mb, na.method = "zero", inf.method = "fixed.noise")
#' @export
NAInfRemoveMultiBlock <- function(MB,
                                  blocks = NULL,
                                  minfrac = 0.5,
                                  na.method = c(
                                    "none", "zero", "median", "discard",
                                    "fixed.value", "fixed.value.all",
                                    "fixed.noise", "random.noise", "QRILC"
                                  )[8],
                                  inf.method = c("none", "fixed.noise", "random.noise")[3],
                                  constant = 0,
                                  factor.NA = 0.5,
                                  sd.noise = 0.3,
                                  seed.number = NULL,
                                  tune.sigma = 1,
                                  showWarning = TRUE) {
  # ---- Input validation -------------------------------------------------- #

  if (!("MultiBlock" %in% class(MB))) {
    stop("'MB' is not a MultiBlock object.")
  }

  if (is.null(na.method) || length(na.method) > 1) {
    stop(
      "Provide a single na.method: 'none', 'zero', 'median', 'discard', ",
      "'fixed.value', 'fixed.value.all', 'fixed.noise', 'random.noise', or 'QRILC'."
    )
  }
  na.method <- tolower(na.method)
  if (!(na.method %in% c(
    "none", "zero", "median", "discard", "fixed.value",
    "fixed.value.all", "fixed.noise", "random.noise", "qrilc"
  ))) {
    stop(
      "Invalid na.method. Choose from: 'none', 'zero', 'median', 'discard', ",
      "'fixed.value', 'fixed.value.all', 'fixed.noise', 'random.noise', 'QRILC'."
    )
  }
  if (na.method == "qrilc" && !requireNamespace("imputeLCMD", quietly = TRUE)) {
    stop("Package 'imputeLCMD' is required for na.method = 'QRILC'. ",
         "Please install it with: install.packages('imputeLCMD').")
  }

  if (is.null(inf.method) || length(inf.method) > 1) {
    stop("Provide a single inf.method: 'none', 'fixed.noise', or 'random.noise'.")
  }
  inf.method <- tolower(inf.method)
  if (!(inf.method %in% c("none", "fixed.noise", "random.noise"))) {
    stop("Invalid inf.method. Choose from: 'none', 'fixed.noise', 'random.noise'.")
  }

  # ---- Resolve block indices --------------------------------------------- #

  if (!is.null(blocks)) {
    if (is.character(blocks)) {
      block.ids <- match(blocks, names(MB@Data))
      missing_b <- blocks[is.na(block.ids)]
      if (length(missing_b) > 0) {
        warning(sprintf(
          "Block name(s) not found in the MultiBlock: %s",
          paste(missing_b, collapse = ", ")
        ))
        block.ids <- block.ids[!is.na(block.ids)]
      }
    } else if (is.numeric(blocks)) {
      block.ids <- intersect(seq_along(MB@Data), blocks)
      missing_b <- setdiff(blocks, seq_along(MB@Data))
      if (length(missing_b) > 0) {
        warning(sprintf(
          "Block index/indices not found in the MultiBlock: %s",
          paste(missing_b, collapse = ", ")
        ))
      }
    } else {
      stop("'blocks' must be a character or numeric vector.")
    }
  } else {
    block.ids <- seq_along(MB@Data)
  }

  if (length(block.ids) == 0) {
    warning("0 blocks were selected. Data was not modified.")
    return(MB)
  }

  # ---- Validate and convert minfrac to an absolute count ----------------- #

  if (is.null(minfrac) || minfrac < 0) {
    if (showWarning) warning("'minfrac' was set to 0.")
    minfrac <- 0
  } else if (minfrac > 1) {
    if (showWarning) warning("'minfrac' was set to 1.")
    minfrac <- 1
  }
  minfrac_n <- round(nrow(MB@Data[[block.ids[1]]]) * minfrac)

  # ---- Set seed once for all random steps -------------------------------- #

  if (!is.null(seed.number) &&
    (na.method == "random.noise" || inf.method == "random.noise")) {
    old_seed <- .Random.seed
    on.exit(
      {
        .Random.seed <<- old_seed
      },
      add = TRUE
    )
    set.seed(seed.number)
  }

  # ---- Step 1: discard variables with too many NA / Inf values ----------- #

  n_discarded <- 0L
  to_delete <- integer(0)

  for (i in block.ids) {
    if (any(is.complex(MB@Data[[i]]))) {
      stop(sprintf(
        "Block '%s' contains complex values that cannot be processed.",
        names(MB@Data)[i]
      ))
    }
    if (any(is.logical(MB@Data[[i]]))) {
      stop(sprintf(
        "Block '%s' contains logical values that cannot be processed.",
        names(MB@Data)[i]
      ))
    }

    bad_n <- apply(MB@Data[[i]], 2, function(x) sum(is.na(x) | is.infinite(x)))
    drop_idx <- which(bad_n >= minfrac_n)

    if (length(drop_idx) > 0) {
      keep_idx <- setdiff(seq_len(ncol(MB@Data[[i]])), drop_idx)
      MB@Data[[i]] <- MB@Data[[i]][, keep_idx, drop = FALSE]
      MB@Variables[[i]] <- MB@Variables[[i]][keep_idx]
      n_discarded <- n_discarded + length(drop_idx)
    }

    if (ncol(MB@Data[[i]]) == 0) {
      to_delete <- c(to_delete, i)
    }
  }

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
      stop("All selected blocks contained only NA/Infinite values and were removed.")
    }
    block.ids <- match(kept_names, names(MB@Data))
  }

  if (n_discarded > 0L && showWarning) {
    warning(sprintf(
      "%d variable(s) exceeded the NA/Infinite threshold and were deleted.",
      n_discarded
    ))
  }

  # ---- Step 2: replace remaining Infinite values ------------------------- #

  if (inf.method != "none") {
    for (i in block.ids) {
      for (nc in seq_len(ncol(MB@Data[[i]]))) {
        col <- MB@Data[[i]][, nc]
        if (!any(is.infinite(col))) next

        finite_vals <- col[is.finite(col) & !is.na(col)]
        col_min <- if (length(finite_vals) > 0) min(finite_vals) else 0
        col_max <- if (length(finite_vals) > 0) max(finite_vals) else 0
        neg_pos <- which(is.infinite(col) & col < 0)
        pos_pos <- which(is.infinite(col) & col > 0)

        if (inf.method == "fixed.noise") {
          if (length(neg_pos) > 0) {
            MB@Data[[i]][neg_pos, nc] <- col_min * factor.NA
          }
          if (length(pos_pos) > 0) {
            MB@Data[[i]][pos_pos, nc] <- col_max * (1 / factor.NA)
          }
        } else { # random.noise
          if (length(neg_pos) > 0) {
            m <- col_min * factor.NA
            MB@Data[[i]][neg_pos, nc] <- rnorm(length(neg_pos),
              mean = m,
              sd   = abs(m * sd.noise)
            )
          }
          if (length(pos_pos) > 0) {
            m <- col_max * (1 / factor.NA)
            MB@Data[[i]][pos_pos, nc] <- rnorm(length(pos_pos),
              mean = m,
              sd   = abs(m * sd.noise)
            )
          }
        }
      }
    }
  }

  # ---- Step 3: impute remaining NA values -------------------------------- #

  if (na.method == "none") {
    return(MB)
  }

  for (i in block.ids) {
    if (na.method == "zero") {
      MB@Data[[i]][is.na(MB@Data[[i]])] <- 0
    } else if (na.method == "median") {
      for (nc in seq_len(ncol(MB@Data[[i]]))) {
        na_pos <- which(is.na(MB@Data[[i]][, nc]))
        if (length(na_pos) > 0) {
          MB@Data[[i]][na_pos, nc] <- median(MB@Data[[i]][, nc], na.rm = TRUE)
        }
      }
    } else if (na.method == "fixed.value") {
      MB@Data[[i]][is.na(MB@Data[[i]])] <- constant
    } else if (na.method == "fixed.value.all") {
      MB@Data[[i]] <- MB@Data[[i]] + constant # adds to non-NA elements
      MB@Data[[i]][is.na(MB@Data[[i]])] <- constant # fills former NAs
    } else if (na.method == "discard") {
      has_na <- apply(MB@Data[[i]], 2, anyNA)
      keep_cols <- which(!has_na)
      if (length(keep_cols) == 0) {
        nm <- names(MB@Data)[i]
        if (nm %in% names(MB@Batch)) MB@Batch[[nm]] <- NULL
        if (nm %in% names(MB@Metadata)) MB@Metadata[[nm]] <- NULL
        MB@Data[[i]] <- NULL
        MB@Variables[[i]] <- NULL
        warning(sprintf("Block '%s' was removed because all columns still contained NA values.", nm))
      } else {
        MB@Data[[i]] <- MB@Data[[i]][, keep_cols, drop = FALSE]
        MB@Variables[[i]] <- MB@Variables[[i]][keep_cols]
      }
    } else if (na.method == "fixed.noise") {
      for (nc in seq_len(ncol(MB@Data[[i]]))) {
        na_pos <- which(is.na(MB@Data[[i]][, nc]))
        if (length(na_pos) > 0) {
          col_min <- min(MB@Data[[i]][, nc], na.rm = TRUE)
          MB@Data[[i]][na_pos, nc] <- col_min * factor.NA
        }
      }
    } else if (na.method == "random.noise") {
      for (nc in seq_len(ncol(MB@Data[[i]]))) {
        na_pos <- which(is.na(MB@Data[[i]][, nc]))
        if (length(na_pos) > 0) {
          col_min <- min(MB@Data[[i]][, nc], na.rm = TRUE)
          m <- col_min * factor.NA
          MB@Data[[i]][na_pos, nc] <- rnorm(length(na_pos),
            mean = m,
            sd   = abs(m * sd.noise)
          )
        }
      }
    } else if (na.method == "qrilc") {
      for (nc in seq_len(ncol(MB@Data[[i]]))) {
        if (anyNA(MB@Data[[i]][, nc])) {
          MB@Data[[i]][, nc] <- imputeLCMD::impute.QRILC(
            as.matrix(MB@Data[[i]][, nc]),
            tune.sigma = tune.sigma
          )[[1]]
        }
      }
    }
  }

  return(MB)
}
