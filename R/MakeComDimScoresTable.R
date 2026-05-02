#' MakeComDimScoresTable
#'
#' Creates a long (tidy) data frame with the global and/or local scores from a ComDim model,
#' suitable for use with \pkg{ggplot2}.
#' @param model The output from a ComDim analysis (a \code{ComDim} object).
#' @param blocks The blocks from which local scores will be extracted. A vector of integers (block
#'   indices) or block names. When omitted, all blocks are included. Only relevant when
#'   \code{"T.scores"} or \code{"T.scores.ort"} are in \code{include}.
#' @param dim Integer vector of predictive component indices to include. When omitted, all
#'   predictive components in the model are included.
#' @param dim.ort Integer vector of orthogonal component indices to include. When \code{NULL},
#'   all orthogonal components present in the model are included. When \code{0}, no orthogonal
#'   components are included.
#' @param include Character vector selecting which score types to include. Accepted values
#'   (case-insensitive) are \code{"Q.scores"} (global predictive scores),
#'   \code{"T.scores"} (local predictive scores), \code{"Q.scores.ort"} (global orthogonal
#'   scores), and \code{"T.scores.ort"} (local orthogonal scores). Defaults to
#'   \code{c("Q.scores", "T.scores")}.
#' @return A long data frame with one row per sample–component–score-type combination, containing
#'   the following columns:
#'   \describe{
#'     \item{\code{sample.id}}{Sample name (factor).}
#'     \item{\code{sample.id.number}}{Integer position of the sample (factor).}
#'     \item{\code{block.id}}{Block index, or \code{"Global"} / \code{"Global.ort"} for Q scores.}
#'     \item{\code{block.name}}{Block name, or \code{"Global"} / \code{"Global.ort"} for Q scores (factor).}
#'     \item{\code{dim}}{Component number.}
#'     \item{\code{scores.type}}{One of \code{"Global"}, \code{"Global.ort"}, \code{"Local"}, or \code{"Local.ort"} (factor).}
#'     \item{\code{scores.type.dim}}{Concatenation of scores type and component number, e.g. \code{"Q.scores1"} (factor).}
#'     \item{\code{value}}{Score value.}
#'   }
#' @seealso \code{\link{MakeComDimLoadingsTable}}, \code{\link{ComDim_PCA}}
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50) # 10 rows and 50 columns
#' b2 <- matrix(rnorm(800), 10, 80) # 10 rows and 80 columns
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' model <- ComDim_PCA(mb, ndim = 2)
#' tbl <- MakeComDimScoresTable(model)
#' @export
MakeComDimScoresTable <- function(model,
                                  blocks = NULL,
                                  dim = NULL,
                                  dim.ort = NULL,
                                  include = c("Q.scores", "T.scores", "Q.scores.ort", "T.scores.ort")[1:2]) {
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

  if (is.null(include)) {
    include <- c("Q.scores", "T.scores")
  } else { # Added to match words regardless upper or lowercase
    include <- toupper(include)
    include <- gsub("T.SCORES", "T.scores", include)
    include <- gsub("Q.SCORES", "Q.scores", include)
    include <- gsub("ORT", "ort", include)
  }

  if (any(!(include %in% c("Q.scores", "T.scores", "Q.scores.ort", "T.scores.ort")))) {
    stop("Some elements in 'include' were not recognized'.")
  }


  if (!(is.null(blocks))) {
    if (is.character(blocks)) {
      block.ids <- match(blocks, names(model@T.scores))
      if (any(is.na(block.ids))) {
        warning(
          "Some block names were not found in the MultiBlock: %s",
          paste0(blocks[which(is.na(block.ids))], sep = ", ")
        )
        block.ids <- block.ids[!is.na(block.ids)]
      }
    } else if (is.numeric(blocks)) {
      block.ids <- intersect(seq_along(names(model@T.scores)), blocks)
      if (length(blocks) != length(block.ids)) {
        warning(
          "Some block names were not found in the MultiBlock: %s",
          paste0(setdiff(seq_along(names(model@T.scores)), blocks), sep = ", ")
        )
      }
    }
  } else {
    block.ids <- seq(1, length(names(model@T.scores)))
  }

  if (length(block.ids) == 0) {
    stop("0 blocks were selected. The table cannot be built.")
  }

  # Check number of samples is always the same across all Q and T.
  nsamples <- nrow(model@Q.scores)
  for (i in block.ids) {
    if (nsamples != nrow(model@T.scores[[i]])) {
      stop("'model' is not an output from a ComDim analysis.")
    }
  }

  df <- NULL
  if ("Q.scores" %in% include) {
    for (v in dim) { # For each of the dimensions
      df2 <- data.frame(
        sample.id = rownames(model@Q.scores),
        sample.id.number = seq(1, nsamples),
        block.id = rep("Global", nsamples),
        block.name = rep("Global", nsamples),
        dim = rep(v, nsamples),
        scores.type = rep("Global", nsamples),
        scores.type.dim = paste0(rep("Q.scores", nsamples), rep(v, nsamples)),
        value = model@Q.scores[, v]
      )
      if (is.null(df)) {
        df <- df2
      } else {
        df <- rbind.data.frame(df, df2)
      }
    }
  }

  if (dim.ort != 0) {
    if ("Q.scores.ort" %in% include) {
      for (v in dim.ort) { # For each of the dimensions
        df2 <- data.frame(
          sample.id = rownames(model@Q.scores),
          sample.id.number = seq(1, nsamples),
          block.id = rep("Global.ort", nsamples),
          block.name = rep("Global.ort", nsamples),
          dim = rep(v, nsamples),
          scores.type = rep("Global.ort", nsamples),
          scores.type.dim = paste0(rep("Q.scores.ort", nsamples), rep(v, nsamples)),
          value = model@Orthogonal$Q.scores[, v]
        )
        if (is.null(df)) {
          df <- df2
        } else {
          df <- rbind.data.frame(df, df2)
        }
      }
    }
  }

  if ("T.scores" %in% include) {
    for (i in block.ids) { # For each of the dimensions
      for (v in dim) {
        df2 <- data.frame(
          sample.id = rownames(model@Q.scores),
          sample.id.number = seq(1, nsamples),
          block.id = rep(i, nsamples),
          block.name = rep(names(model@T.scores)[i], nsamples),
          dim = rep(v, nsamples),
          scores.type = rep("Local", nsamples),
          scores.type.dim = paste0(rep("T.scores", nsamples), rep(v, nsamples)),
          value = model@T.scores[[i]][, v]
        )
        if (is.null(df)) {
          df <- df2
        } else {
          df <- rbind.data.frame(df, df2)
        }
      }
    }
  }

  if (dim.ort != 0) {
    if ("T.scores.ort" %in% include) {
      for (i in block.ids) { # For each of the dimensions
        for (v in dim.ort) {
          df2 <- data.frame(
            sample.id = rownames(model@Q.scores),
            sample.id.number = seq(1, nsamples),
            block.id = rep(i, nsamples),
            block.name = rep(names(model@Orthogonal$T.scores)[i], nsamples),
            dim = rep(v, nsamples),
            scores.type = rep("Local.ort", nsamples),
            scores.type.dim = paste0(rep("T.scores.ort", nsamples), rep(v, nsamples)),
            value = model@Orthogonal$T.scores[[i]][, v]
          )
          if (is.null(df)) {
            df <- df2
          } else {
            df <- rbind.data.frame(df, df2)
          }
        }
      }
    }
  }
  rownames(df) <- NULL

  # Add factors to data
  df$sample.id <- factor(df$sample.id, levels = df$sample.id[seq(1, nsamples)])
  df$sample.id.number <- factor(df$sample.id.number, levels = c("Global", seq(1, nsamples)))
  df$block.name <- factor(df$block.name, levels = c(names(model@T.scores), "Global", "Global.ort"))
  df$scores.type <- factor(df$scores.type, levels = unique(df$scores.type))
  df$scores.type.dim <- factor(df$scores.type.dim, levels = unique(df$scores.type.dim))

  return(df)
}
