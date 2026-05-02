## All class definitions should go in here.

############################################################
## Class unions
setClassUnion("data.frameNULL", c("data.frame", "NULL"))
setClassUnion("matrixNULL", c("matrix", "NULL"))
setClassUnion("vectorNULL", c("vector", "NULL", "character"))
setClassUnion("vectorlist", c("vector", "list"))
setClassUnion("listNULL", c("list", "NULL"))

############################################################

## Validate MultiBlock class
#' checkMultiBlock
#'
#' Validation function used during MultiBlock building.
#' @param object A MultiBlock object.
#' @return TRUE if valid, or a character vector of error messages.
#' @keywords internal
checkMultiBlock <- function(object) {
  errors <- character()

  if (any(!(names(object@Batch)) %in% names(object@Data))) {
    msg <- "
    There is at least a Batch vector with no corresponding Data.
    Check Batch and Data names."
    errors <- c(msg, errors)
  }
  if (any(!(names(object@Metadata)) %in% names(object@Data))) {
    msg <- "
    There is at least a Metadata data.frame with no corresponding Data.
    Check Batch and Data names."
    errors <- c(msg, errors)
  }
  if (any(!(names(object@Variables)) %in% names(object@Data))) {
    msg <- "
    There is at least an element in the list of Variables with no
    corresponding Data. Check Batch and variable block names."
    errors <- c(msg, errors)
  }
  samples_is_list <- is.list(object@Samples)

  for (i in seq_along(object@Data)) {
    block_name <- names(object@Data)[i]

    # Check sample count — handle both vector Samples and per-block list Samples
    if (samples_is_list) {
      if (!(block_name %in% names(object@Samples))) {
        msg <- sprintf("
      Block '%s' has no corresponding entry in the Samples list.", block_name)
        errors <- c(msg, errors)
      } else if (nrow(object@Data[[i]]) != length(object@Samples[[block_name]])) {
        msg <- "
      The number of samples is not correct for at least one block."
        errors <- c(msg, errors)
      }
    } else {
      if (nrow(object@Data[[i]]) != length(object@Samples)) {
        msg <- "
      The number of samples is not correct for at least one block. Use
      Use MultiBlock() with 'ignore.size = TRUE' to build
      a MultiBlock with blocks of different sample number."
        errors <- c(msg, errors)
      }
    }

    if (ncol(object@Data[[i]]) != length(object@Variables[[i]])) {
      msg <- "
      The number of variables is not correct for at least one block."
      errors <- c(msg, errors)
    }

    if (names(object@Data)[i] %in% names(object@Batch)) {
      if (nrow(object@Data[[i]]) != length(object@Batch[[i]])) {
        msg <- "
        The length of Batch is not correct for at least one block."
        errors <- c(msg, errors)
      }
    }
    if (names(object@Data)[i] %in% names(object@Metadata)) {
      if (nrow(object@Data[[i]]) != nrow(object@Metadata[[i]])) {
        msg <- "
        The length of Metadata is not correct for at least one block."
        errors <- c(msg, errors)
      }
    }

    if (length(errors) != 0) {
      break
    }
  }

  if (length(errors) == 0) TRUE else errors
}


############################################################

## MultiBlock
#'
#' @title MultiBlock object
#'
#' @description Object of the type \code{MultiBlock}, to use as input for ComDim
#' analyses.
#'
#' @slot Samples vector with the sample names. If this data is not available,
#' the slot will be filled with integers.
#'
#' @slot Data A list with the data-blocks.
#'
#' @slot Variables A character vector with the variable names. If this data is
#' not available, the slot will be filled with integers.
#'
#' @slot Batch A list with the vectors with the batch information for each
#' data-block. Optional.
#'
#' @slot Metadata A list with the samples metadata.
#'
#' @name ComDim-classes
#' @importFrom methods new setClass setClassUnion setGeneric setMethod validObject is

setClass("MultiBlock",
  representation(
    Samples = "vectorlist",
    Data = "list",
    Variables = "list", # List of vectors
    Batch = "list",
    Metadata = "list"
  ),
  prototype(
    Samples = vector(),
    Data = list(),
    Variables = list(),
    Batch = list(),
    Metadata = list()
  ),
  validity = checkMultiBlock
)


############################################################
## ComDim
#' @aliases ComDim
#'
#' @title ComDim object
#'
#' @description The output of a ComDim analysis.
#'
#' @slot Method The algorithm used in the core of the ComDim analysis (ex. PCA,
#' PLS,...)
#'
#' @slot ndim The number of components.
#'
#' @slot Q.scores The Global scores.
#'
#' @slot T.scores The Local scores.
#'
#' @slot P.loadings The Loadings
#'
#' @slot Saliences The Saliences
#'
#' @slot R2X The explained variance of the MultiBlock.
#'
#' @slot R2Y For regression or discriminant models, the explained variance of
#' the Y-block.
#'
#' @slot Q2 For regression or discriminant models, the predicted variance of
#' the Y-block.
#'
#' @slot DQ2 For discriminant models, the predicted discriminant variance of the
#'  Y-block.
#'
#' @slot Singular The singular values.
#'
#' @slot Mean The mean values for each variable in the MultiBlock.
#'
#' @slot Norm The norm values for each variable in the MultiBlock.
#'
#' @slot PLS.model For ComDim analyses using PLS as the core algorithm, contains
#' the W, B, B0 and Y matrices.
#'
#' @slot cv For ComDim_KOPLS, it contains the index of the samples used during
#' the cross-validation.
#'
#' @slot Prediction A list with the predicted Y, the decision rule used for
#' sample classification, the sensitivity, the specificity, and a confusion
#' matrix.
#'
#' @slot Metadata A list with the per-sample metadata for each block.
#'
#' @slot variable.block A vector with the same length as the P.loadings,
#' indicating the block each variable belongs to.
#'
#' @slot runtime The used running time.
#'
#' @name ComDim-classes
#' @importFrom methods new

setClass(
  "ComDim",
  representation(
    Method = "character",
    ndim = "numeric",
    Q.scores = "matrix",
    T.scores = "list",
    P.loadings = "matrix", # List of P
    Saliences = "matrix",
    Orthogonal = "list",
    VIP = "vector",
    VIP.block = "list",
    R2X = "vector",
    R2Y = "vector",
    Q2 = "vector",
    DQ2 = "vector",
    Singular = "vector",
    Mean = "list",
    Norm = "list",
    PLS.model = "list", # W, U, B, B0, Y
    cv = "list", # cvTrainIndex & cvTestIndex
    Prediction = "list", #  Y.pred, trueClass, sensSpec, confusionMatrix, nclasses, decisionRule
    Metadata = "listNULL",
    variable.block = "vector",
    runtime = "numeric"
  ),
  prototype(
    Method = NULL,
    ndim = NULL,
    Q.scores = NULL,
    T.scores = list(),
    P.loadings = matrix(),
    Saliences = matrix(),
    Orthogonal = list(),
    VIP = vector(),
    VIP.block = list(),
    R2X = vector(),
    R2Y = vector(),
    Q2 = vector(),
    DQ2 = vector(),
    Singular = vector(),
    Mean = list(),
    Norm = list(),
    PLS.model = NULL, # W, U, B, B0, Y, Y.Pred,
    cv = NULL, # cvTrainIndex & cvTestIndex
    Prediction = NULL, #  Y.pred, decisionRule, trueClass, sensSpec, confusionMatrix, nclasses, decisionRule
    Metadata = NULL,
    variable.block = NULL,
    runtime = NULL
  )
) # It would be nice to add a checkComDim like checkMultiBlock

############################################################
## S4 Generics for MultiBlock accessors

#' @rdname blockNames-MultiBlock-method
#' @export
setGeneric("blockNames", function(x, ...) methods::standardGeneric("blockNames"))
#' @rdname blockNames-set-MultiBlock-method
#' @export
setGeneric("blockNames<-", function(x, value) methods::standardGeneric("blockNames<-"))

#' @rdname sampleNames-MultiBlock-method
#' @export
setGeneric("sampleNames", function(x) methods::standardGeneric("sampleNames"))
#' @rdname sampleNames-set-MultiBlock-method
#' @export
setGeneric("sampleNames<-", function(x, value) methods::standardGeneric("sampleNames<-"))

#' @rdname variableNames-MultiBlock-method
#' @export
setGeneric("variableNames", function(x, ...) methods::standardGeneric("variableNames"))
#' @rdname variableNames-set-MultiBlock-method
#' @export
setGeneric("variableNames<-", function(x, value) methods::standardGeneric("variableNames<-"))

############################################################
## S4 Methods for MultiBlock

# ---- blockNames ----

#' blockNames
#'
#' Return the block names of a MultiBlock object.
#' @param x A MultiBlock object.
#' @param ... Not used. Present for S4 generic dispatch compatibility.
#' @param slot A string: \code{"Data"} (default), \code{"Batch"}, or
#'   \code{"Metadata"}, indicating which slot to retrieve names from.
#' @return A character vector with the block names of the requested slot.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' blockNames(mb) # c("b1", "b2")
#' blockNames(mb, "Data") # same
#' @export
setMethod("blockNames", "MultiBlock", function(x, slot = "Data") {
  if (length(slot) > 1) slot <- "Data"
  switch(tolower(slot),
    data     = names(x@Data),
    batch    = names(x@Batch),
    metadata = names(x@Metadata),
    stop("'slot' must be 'Data', 'Batch', or 'Metadata'.")
  )
})

#' blockNames<-
#'
#' Set the block names of a MultiBlock object.
#' Renames the Data and Variables slots, and updates Batch and Metadata
#' names to stay consistent.
#' @param x A MultiBlock object.
#' @param value A character vector of new block names (same length as the
#'   number of blocks).
#' @return The updated MultiBlock object.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' blockNames(mb) <- c("block1", "block2")
#' blockNames(mb) # c("block1", "block2")
#' @export
setMethod("blockNames<-", "MultiBlock", function(x, value) {
  if (length(names(x@Batch)) != 0) {
    names(x@Batch) <- value[match(names(x@Batch), names(x@Data))]
  }
  if (length(names(x@Metadata)) != 0) {
    names(x@Metadata) <- value[match(names(x@Metadata), names(x@Data))]
  }
  names(x@Data) <- value
  names(x@Variables) <- value
  x
})

# ---- sampleNames ----

#' sampleNames
#'
#' Return the sample names of a MultiBlock object.
#' @param x A MultiBlock object.
#' @return A vector with the sample names.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' rownames(b1) <- paste0("s", 1:10)
#' mb <- MultiBlock(Data = list(b1 = b1))
#' sampleNames(mb) # "s1" ... "s10"
#' @export
setMethod("sampleNames", "MultiBlock", function(x) {
  x@Samples
})

#' sampleNames<-
#'
#' Set the sample names of a MultiBlock object.
#' @param x A MultiBlock object.
#' @param value A vector of new sample names. Must have the same length as
#'   the current number of samples.
#' @return The updated MultiBlock object.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' mb <- MultiBlock(Data = list(b1 = b1))
#' sampleNames(mb) <- paste0("patient_", 1:10)
#' sampleNames(mb)
#' @export
setMethod("sampleNames<-", "MultiBlock", function(x, value) {
  if (!is.vector(value)) {
    warning("The provided sample names are not a vector. Sample names were not changed.")
    return(x)
  }
  if (length(x@Samples) != length(value)) {
    warning("The length of the provided vector does not match. Sample names were not changed.")
    return(x)
  }
  x@Samples <- value
  x
})

# ---- variableNames ----

#' variableNames
#'
#' Return the variable names of a MultiBlock object.
#' @param x A MultiBlock object.
#' @param ... Not used. Present for S4 generic dispatch compatibility.
#' @param block Optional. A vector of block names or indices to retrieve.
#'   When omitted, variable names for all blocks are returned as a named list.
#'   When a single block is specified, a plain vector is returned.
#' @return A named list of variable-name vectors (all blocks), or a single
#'   vector when exactly one block is requested.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' variableNames(mb) # named list: b1 = 1:50, b2 = 1:80
#' variableNames(mb, "b1") # 1:50
#' variableNames(mb, 1:2) # same as variableNames(mb)
#' @export
setMethod("variableNames", "MultiBlock", function(x, block) {
  if (missing(block)) block <- seq_along(x@Data)
  if (is.list(block)) block <- unlist(block)

  if (length(block) == 1) {
    x@Variables[[block]]
  } else {
    thevars <- list()
    if (is.character(block)) {
      for (i in block) thevars[[i]] <- x@Variables[[i]]
    } else {
      bn <- blockNames(x)
      for (i in block) thevars[[bn[i]]] <- x@Variables[[bn[i]]]
    }
    thevars
  }
})

#' variableNames<-
#'
#' Set the variable names of a MultiBlock object.
#' \code{value} must be a named list with one entry per block, where each
#' entry is a vector of variable names whose length matches the number of
#' columns in that block.
#'
#' To update a single block, use standard list-replacement chaining:
#' \code{variableNames(mb)[["b1"]] <- newNames}
#'
#' @param x A MultiBlock object.
#' @param value A named list of variable-name vectors, one per block.
#' @return The updated MultiBlock object.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#'
#' # Replace all at once:
#' variableNames(mb) <- list(
#'   b1 = paste0("v", 1:50),
#'   b2 = paste0("v", 1:80)
#' )
#'
#' # Replace a single block using chaining:
#' variableNames(mb)[["b1"]] <- paste0("feat_", 1:50)
#' @export
setMethod("variableNames<-", "MultiBlock", function(x, value) {
  if (!is.list(value)) {
    stop("'value' must be a named list of variable-name vectors, one per block.")
  }
  if (is.null(names(value)) || any(names(value) == "")) {
    stop("All elements of 'value' must be named.")
  }
  for (nm in names(value)) {
    if (!(nm %in% names(x@Data))) {
      stop(sprintf("Block '%s' not found in the MultiBlock.", nm))
    }
    if (length(value[[nm]]) != ncol(x@Data[[nm]])) {
      stop(sprintf(
        "The length of variable names for block '%s' (%d) does not match ncol (%d).",
        nm, length(value[[nm]]), ncol(x@Data[[nm]])
      ))
    }
    x@Variables[[nm]] <- value[[nm]]
  }
  x
})

# ---- ncol / nrow ----

#' ncol for MultiBlock
#'
#' Return the number of columns (variables) in each block of a MultiBlock.
#' @param x A MultiBlock object.
#' @return A named integer vector with the number of columns per block.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' ncol(mb) # c(b1 = 50, b2 = 80)
#' @export
setMethod("ncol", "MultiBlock", function(x) {
  sapply(x@Data, ncol)
})

#' nrow for MultiBlock
#'
#' Return the number of rows (samples) in each block of a MultiBlock.
#' For a valid MultiBlock all blocks share the same number of rows, so the
#' result is a named integer vector with one (identical) value per block.
#' @param x A MultiBlock object.
#' @return A named integer vector with the number of rows per block.
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' nrow(mb) # c(b1 = 10, b2 = 10)
#' @export
setMethod("nrow", "MultiBlock", function(x) {
  sapply(x@Data, nrow)
})
