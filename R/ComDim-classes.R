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
checkMultiBlock <- function(object){

  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Validation function used during MultiBlock building.
  #' @param object A MultiBlock object.
  #' @return The building of the MultiBlock is interrupted if an error is found.
  #' @export

  errors <- character()

  if(any(!(names(object@Batch)) %in% names(object@Data))){
    msg <- "
    There is at least a Batch vector with no corresponding Data.
    Check Batch and Data names."
    errors <- c(msg, errors)
  }
  if(any(!(names(object@Metadata)) %in% names(object@Data))){
    msg <- "
    There is at least a Metadata data.frame with no corresponding Data.
    Check Batch and Data names."
    errors <- c(msg, errors)
  }
  if(any(!(names(object@Variables)) %in% names(object@Data))){
    msg <- "
    There is at least an element in the list of Variables with no
    corresponding Data. Check Batch and variable block names."
    errors <- c(msg, errors)
  }
  for(i in 1:length(object@Data)){
    if(nrow(object@Data[[i]]) != length(object@Samples)){
      msg <- "
      The number of samples is not correct for at least one block. Use
      BuildMultiBlock() with 'ignore.size = TRUE' to build a MultiBlock
      with blocks of different sample number."
      errors <- c(msg, errors)
    }
    if(ncol(object@Data[[i]]) != length(object@Variables[[i]])){
      msg <- "
      The number of variables is not correct for at least one block."
      errors <- c(msg, errors)
    }

    if(names(object@Data)[i] %in% names(object@Batch)){
      if(nrow(object@Data[[i]]) != length(object@Batch[[i]])){
        msg <- "
        The length of Batch is not correct for at least one block."
        errors <- c(msg, errors)
      }
    }
    if(names(object@Data)[i] %in% names(object@Metadata)){
      if(nrow(object@Data[[i]]) != nrow(object@Metadata[[i]])){
        msg <- "
        The length of Metadata is not correct for at least one block."
        errors <- c(msg, errors)
      }
    }

    if(length(errors) != 0) {
      break
    }
  }

  if(length(errors) == 0 ) TRUE else errors
}


############################################################

## MultiBlock
#' @aliases MultiBlock
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
#' data-block. Facultative.
#'
#' @slot Metadata A list with the samples metadata.
#'
#' @name ComDim-classes
#' @importFrom methods new

setClass("MultiBlock",
               representation(Samples = "vectorlist",
                              Data = "list",
                              Variables = "list", # List of vectors
                              Batch = "list",
                              Metadata = "list"),
               prototype(Samples = vector(),
                         Data = list(),
                         Variables = list(),
                         Batch = list(),
                         Metadata = list()),
               validity = checkMultiBlock)


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
#' @slot PLS.model For ComDim analyses using PLS as core algorithm, it contains
#' the W, B, B0 and Y).
#'
#' @slot cv For ComDim_KOPLS, it contains the index of the samples used during
#' the cross-validation.
#'
#' @slot Prediction A list with the predicted Y, the decision rule used for
#' sample classification, the sensitivity, the specificity, and a confusion
#' matrix.
#'
#' @slot Metadata vector with the sample names. If this data is not available,
#' the slot will be filled with integers.
#'
#' @slot variable.block A vector with the same length than the P.loadings
#' indicating the block each variable belongs to.
#'
#' @slot runtime The used running time.
#'
#' @name ComDim-classes
#' @importFrom methods new

setClass("ComDim",
         representation(Method = "character",
                        ndim = "numeric",
                        Q.scores = "matrix",
                        T.scores = "list",
                        P.loadings = "matrix", # List of P
                        Saliences = "matrix",
                        Orthogonal = "list",
                        R2X = "vector",
                        R2Y = "vector",
                        Q2 = "vector",
                        DQ2 = "vector",
                        Singular = "vector",
                        Mean = "list",
                        Norm = "list",
                        PLS.model = "list", # W, U, B, B0, Y
                        cv = "list", #cvTrainIndex & cvTestIndex
                        Prediction = "list", #  Y.pred, trueClass, sensSpec, confusionMatrix, nclasses, decisionRule
                        Metadata = "listNULL",
                        variable.block = "vector",
                        runtime = "numeric"),
         prototype(Method = NULL,
                   ndim = NULL,
                   Q.scores = NULL,
                   T.scores = list(),
                   P.loadings = matrix(),
                   Saliences = matrix(),
                   Orthogonal = list(),
                   R2X = vector(),
                   R2Y = vector(),
                   Q2 = vector(),
                   DQ2 = vector(),
                   Singular = vector(),
                   Mean = list(),
                   Norm = list(),
                   PLS.model = NULL, # W, U, B, B0, Y, Y.Pred,
                   cv = NULL, #cvTrainIndex & cvTestIndex
                   Prediction = NULL, #  Y.pred, decisionRule, trueClass, sensSpec, confusionMatrix, nclasses, decisionRule
                   Metadata = NULL,
                   variable.block = NULL,
                   runtime = NULL)
) # It would be nice to add a  checkComDim like I did with checkMultiBlock
