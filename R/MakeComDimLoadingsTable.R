#' MakeComDimLoadingsTable
#'
#' Creates a long (tidy) data frame with the local P-loadings from a ComDim model, suitable for
#' use with \pkg{ggplot2}.
#' @param model The output from a ComDim analysis (a \code{ComDim} object).
#' @param blocks The blocks from which loadings will be extracted. A vector of integers (block
#'   indices) or block names. When omitted, all blocks are included.
#' @param dim Integer vector of predictive component indices to include. When omitted, all
#'   predictive components in the model are included.
#' @param dim.ort Integer vector of orthogonal component indices to include. When \code{NULL},
#'   all orthogonal components present in the model are included. When \code{0}, no orthogonal
#'   components are included.
#' @return A long data frame with one row per variable–component combination, containing the
#'   following columns:
#'   \describe{
#'     \item{\code{variable.id}}{Variable name (factor).}
#'     \item{\code{variable.id.number}}{Position of the variable across all blocks (factor).}
#'     \item{\code{block.id}}{Integer index of the block.}
#'     \item{\code{block.name}}{Name of the block (factor).}
#'     \item{\code{dim}}{Component number.}
#'     \item{\code{value}}{Loading value (from \code{P.loadings} or \code{Orthogonal$P.loadings.ort}).}
#'   }
#' @seealso \code{\link{MakeComDimScoresTable}}, \code{\link{ComDim_PCA}}
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50) # 10 rows and 50 columns
#' b2 <- matrix(rnorm(800), 10, 80) # 10 rows and 80 columns
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' model <- ComDim_PCA(mb, ndim = 2)
#' tbl <- MakeComDimLoadingsTable(model)
#' @export
MakeComDimLoadingsTable <- function(model,
                                    blocks = NULL,
                                    dim = NULL,
                                    dim.ort = NULL) {
  if (!inherits(model, "ComDim")) {
    stop("'model' is  not the output of a ComDim analysis.")
  }

  if (is.null(dim)) {
    dim <- 1:model@ndim
  }

  if (is.null(dim.ort) && length(model@Orthogonal) != 0) {
    if (model@Orthogonal$nort != 0) {
      dim.ort <- 1:model@Orthogonal$nort
    }
  } else if (is.null(dim.ort) && length(model@Orthogonal) == 0) {
    dim.ort <- 0
  }

  block_names <- names(model@T.scores)

  if (!(is.null(blocks))) {
    if (is.character(blocks)) {
      block.ids <- match(blocks, block_names)
      if (any(is.na(block.ids))) {
        warning(
          "Some block names were not found in the MultiBlock: %s",
          paste0(blocks[which(is.na(block.ids))], sep = ", ")
        )
        block.ids <- block.ids[!is.na(block.ids)]
      }
    } else if (is.numeric(blocks)) {
      block.ids <- intersect(seq_along(block_names), blocks)
      if (length(blocks) != length(block.ids)) {
        warning(
          "Some block names were not found in the MultiBlock: %s",
          paste0(setdiff(seq_along(block_names), blocks), sep = ", ")
        )
      }
    }
  } else {
    block.ids <- seq(1, length(block_names))
  }


  if (length(block.ids) == 0) {
    stop("0 blocks were selected. The table cannot be built.")
  }


  for (i in block.ids) { # For each of the dimensions

    posvars <- which(model@variable.block == block_names[i])

    for (v in dim) {
      df2 <- data.frame(
        variable.id = rownames(model@P.loadings[posvars, ]),
        variable.id.number = posvars,
        block.id = rep(i, length(posvars)),
        block.name = rep(block_names[i], length(posvars)),
        dim = rep(v, length(posvars)),
        value = model@P.loadings[posvars, v]
      )
      if (v == dim[1] && i == block.ids[1]) {
        df <- df2
      } else {
        df <- rbind.data.frame(df, df2)
      }
    }

    if (dim.ort != 0) {
      for (v in dim.ort) {
        df2 <- data.frame(
          variable.id = rownames(model@Orthogonal$P.loadings.ort[posvars, ]),
          variable.id.number = posvars,
          block.id = rep(i, length(posvars)),
          block.name = rep(block_names[i], length(posvars)),
          dim = rep(v, length(posvars)),
          value = model@Orthogonal$P.loadings.ort[posvars, v]
        )
        df <- rbind.data.frame(df, df2)
      }
    }
  }
  rownames(df) <- NULL


  # Add factors to data
  df$variable.id <- factor(df$variable.id, levels = unique(df$variable.id))
  df$variable.id.number <- factor(df$variable.id.number)
  df$block.name <- factor(df$block.name, levels = block_names)

  return(df)
}
