#' MultiBlock
#'
#' Creates a MultiBlock object from a named list of data blocks.
#' @param Samples A vector of sample names shared across all blocks (optional).
#'   When omitted, sample names are taken from the row names of each block.
#'   If no row names exist and all blocks have the same number of rows,
#'   samples are numbered as integers. Use \code{ignore.names} or
#'   \code{ignore.size} for more flexible behaviour.
#' @param Data A named list of matrices or data.frames (one entry per block).
#' @param Variables A named list of variable-name vectors, one per block
#'   (optional). When omitted, column names are taken from each block; if
#'   absent, variables are numbered as integers.
#' @param Batch A named list of batch vectors, one per block (optional).
#' @param Metadata A named list of metadata data.frames, one per block
#'   (optional).
#' @param ignore.names If TRUE, sample names are not checked across blocks.
#'   All blocks must have the same number of rows unless \code{ignore.size}
#'   is also TRUE.
#' @param ignore.size If TRUE (only meaningful when \code{ignore.names = TRUE}),
#'   blocks with different row counts are accepted. The resulting MultiBlock
#'   stores a per-block Samples list and is not directly suitable for ComDim
#'   (use \code{SplitRW()} first).
#' @return A MultiBlock object.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50) # 10 samples, 50 variables
#' b2 <- matrix(rnorm(800), 10, 80) # 10 samples, 80 variables
#' # Minimal call: Samples and Variables are filled in automatically.
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#'
#' # With explicit sample names (enables cross-block alignment):
#' rownames(b1) <- paste0("s", 1:10)
#' rownames(b2) <- paste0("s", 1:10)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#'
#' # Blocks with different row counts (replicate design):
#' b3 <- matrix(rnorm(800), 30, 80)
#' batch_b3 <- c(rep(1, 10), rep(2, 10), rep(3, 10))
#' mb3 <- MultiBlock(
#'   Data = list(b3 = b3), Batch = list(b3 = batch_b3),
#'   ignore.names = TRUE, ignore.size = TRUE
#' )
#' @references Puig-Castellví F, Jouan-Rimbaud Bouveresse D, Mazéas L, Chapleur O, Rutledge DN (2021).
#'   Rearrangement of incomplete multi-omics datasets combined with ComDim for evaluating replicate
#'   cross-platform variability and batch influence.
#'   \emph{Chemometrics and Intelligent Laboratory Systems}, 218, 104422.
#'   \doi{10.1016/j.chemolab.2021.104422}
#' @export
MultiBlock <- function(Samples = NULL,
                       Data,
                       Variables = NULL,
                       Batch = NULL,
                       Metadata = NULL,
                       ignore.names = FALSE,
                       ignore.size = FALSE) {
  if (is.null(Batch)) Batch <- list()
  if (is.null(Metadata)) Metadata <- list()

  # -------------------------------------------------------------------------
  # Expand any MultiBlock objects passed inside Data
  # -------------------------------------------------------------------------
  if (any(vapply(Data, is, logical(1), "MultiBlock"))) {
    expanded <- list()
    for (nm in names(Data)) {
      item <- Data[[nm]]
      if (is(item, "MultiBlock")) {
        for (bn in names(item@Data)) {
          data_block <- item@Data[[bn]]
          # Restore row names from the Samples slot so sample labels survive
          # the dimname-strip in the main processing loop.
          if (!is.list(item@Samples)) {
            rownames(data_block) <- as.character(item@Samples)
          } else if (!is.null(item@Samples[[bn]])) {
            rownames(data_block) <- as.character(item@Samples[[bn]])
          }
          expanded[[bn]] <- data_block
        }
        for (bn in names(item@Batch)) {
          if (is.null(Batch[[bn]])) Batch[[bn]] <- item@Batch[[bn]]
        }
        for (bn in names(item@Metadata)) {
          if (is.null(Metadata[[bn]])) Metadata[[bn]] <- item@Metadata[[bn]]
        }
      } else {
        expanded[[nm]] <- item
      }
    }
    Data <- expanded
  }

  # -------------------------------------------------------------------------
  # Validate Data input
  # -------------------------------------------------------------------------
  if (!is.list(Data) || length(Data) == 0) {
    stop("'Data' must be a non-empty named list of matrices or data.frames.")
  }
  if (is.null(names(Data)) || any(names(Data) == "")) {
    stop("All elements of 'Data' must be named.")
  }

  # -------------------------------------------------------------------------
  # Process each block: convert to matrix, capture row/col names
  # -------------------------------------------------------------------------
  detected_samples <- vector("list", length(Data))
  detected_vars <- vector("list", length(Data))
  block_nrow <- integer(length(Data))

  for (i in seq_along(Data)) {
    block <- Data[[i]]
    nm <- names(Data)[i]

    if (is.data.frame(block)) {
      if (any(sapply(block, is.character))) {
        stop(sprintf("Block '%s' contains character columns.", nm))
      }
      if (!is.null(rownames(block))) detected_samples[[i]] <- rownames(block)
      detected_vars[[i]] <- if (!is.null(colnames(block))) colnames(block) else seq_len(ncol(block))
      Data[[i]] <- as.matrix(block)
      dimnames(Data[[i]]) <- NULL
    } else if (is.matrix(block)) {
      if (!is.null(rownames(block))) detected_samples[[i]] <- rownames(block)
      detected_vars[[i]] <- if (!is.null(colnames(block))) colnames(block) else seq_len(ncol(block))
      dimnames(Data[[i]]) <- NULL
    } else {
      stop(sprintf("Block '%s' must be a matrix or a data.frame.", nm))
    }
    block_nrow[i] <- nrow(Data[[i]])
  }

  # -------------------------------------------------------------------------
  # Build Variables list when not provided by the user
  # -------------------------------------------------------------------------
  if (is.null(Variables)) {
    Variables <- setNames(detected_vars, names(Data))
  }

  # -------------------------------------------------------------------------
  # Build Samples when not provided by the user
  # -------------------------------------------------------------------------
  if (!is.null(Samples)) {
    # User supplied Samples explicitly — use as-is.
    # The validator will catch any size mismatches.
  } else if (ignore.names) {
    # -- ignore.names = TRUE -------------------------------------------------
    # Do not use names for cross-block alignment; only check/enforce row counts.
    if (!ignore.size && length(unique(block_nrow)) != 1) {
      stop(paste(
        "Blocks do not have the same number of rows.",
        "Provide sample names for all blocks to enable matching,",
        "or set 'ignore.size = TRUE'."
      ))
    }

    if (length(unique(block_nrow)) == 1) {
      # All blocks same size: assign integer indices
      Samples <- seq_len(block_nrow[1])
    } else {
      # ignore.size = TRUE: store a per-block Samples list
      Samples <- setNames(
        lapply(seq_along(Data), function(i) {
          if (!is.null(detected_samples[[i]])) {
            detected_samples[[i]]
          } else {
            seq_len(block_nrow[i])
          }
        }),
        names(Data)
      )
      warning(paste(
        "The resulting MultiBlock is not directly suitable for ComDim because",
        "sample numbers are not common across blocks.",
        "If data contains replicates, use SplitRW().",
        "Otherwise provide sample names and use ignore.names = FALSE."
      ))
    }
  } else {
    # -- ignore.names = FALSE (default) --------------------------------------
    # ignore.size is overridden to FALSE when ignore.names is FALSE
    ignore.size <- FALSE

    has_names <- !sapply(detected_samples, is.null)

    if (!any(has_names)) {
      # No block has row names: require equal sizes, assign integers
      if (length(unique(block_nrow)) != 1) {
        stop(paste(
          "Blocks do not have the same number of rows.",
          "Provide sample names for all blocks to enable matching,",
          "or set 'ignore.names = TRUE'."
        ))
      }
      Samples <- seq_len(block_nrow[1])
    } else if (any(has_names) && !all(has_names)) {
      # Some blocks have row names and others do not: ambiguous
      stop(paste(
        "Some blocks have row names and others do not.",
        "Either provide row names for all blocks, or none."
      ))
    } else {
      # All blocks have row names: align on common samples
      common <- detected_samples[[1]]
      for (i in seq_along(detected_samples)[-1]) {
        common <- intersect(common, detected_samples[[i]])
      }
      if (length(common) == 0) {
        stop("There are no common sample names across blocks.")
      }
      if (length(common) != max(block_nrow)) {
        warning("Some samples are not common across blocks. They will be discarded.")
      }
      # Reorder / subset each block to the common sample set
      for (i in seq_along(Data)) {
        idx <- match(common, detected_samples[[i]])
        Data[[i]] <- Data[[i]][idx, , drop = FALSE]
      }
      Samples <- common
    }
  }

  # -------------------------------------------------------------------------
  # Build and validate the MultiBlock object
  # -------------------------------------------------------------------------
  MB <- new("MultiBlock",
    Samples   = Samples,
    Data      = Data,
    Variables = Variables,
    Batch     = Batch,
    Metadata  = Metadata
  )

  if (validObject(MB)) {
    return(MB)
  } else {
    stop("The 'MultiBlock' could not be built.")
  }
}
