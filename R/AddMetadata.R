# Previously known as BuildMultiBlock. Now only used when there is metadata or batch to add.
# This function will change completely since it will be based on MultiBlock().
AddMetadata <- function(newBlock = newBlock, newVars = NULL, newSamples = NULL,
                        batches = NULL, metadata = NULL) {


#' AddMetadata
#'
#' Creates a MultiBlock object containing only one block. This function is used when Batch or Metadata information must be added.
#' @param newBlock The data to be incorporated into the MultiBlock structure. It can be a matrix or a data.frame.
#' @param newVars The variable names (optional). If NULL, variables names will be taken from the column names in newBlock. If none of them exist, variables will be identified with integers.
#' @param newSamples The sample names (optional). The sample replicates should have the same name. If NULL, sample names will be taken from the row names in newBlock if there are. If none of them exist, samples will be identified with integers.
#' @param batches A vector indicating the batch number for every sample (optional).
#' @param metadata A data.frame descriptive of newBlock (optional).
#' @return The MultiBlock
#' @examples
#'  b1 = matrix(rnorm(500),10,50) # 10 rows and 50 columns
#'  b2 = matrix(rnorm(800),10,80) # 10 rows and 80 columns
#'  b3 = matrix(rnorm(700),10,70) # 10 rows and 70 columns
#'  # Build AddMetadata by adding one data block at a time:
#'  mb <- AddMetadata(b1, newSamples = paste0('sample_',1:10))
#'  mb <- AddMetadata(b2, growingMB = mb)
#'  mb <- AddMetadata(b3, growingMB = mb)
#' @export

  if(!(is.matrix(newBlock))){
    if(!(is.data.frame(newBlock))){
      stop('newBlock is not a matrix nor a data frame.')
    } else {
      newBlock <- as.matrix(newBlock)
    }
  }

  # Check newVars
  if (is.null(newVars)) {
    if (is.null(colnames(newBlock))){
      newVars <- seq(1,ncol(newBlock),1)
    } else {
      newVars <- colnames(newBlock)
    }
  } else {
    if (ncol(newBlock) != length(newVars)) {
      print(sprintf('Please use a newVars object with %s items.', ncol(newBlock)))
      stop('The length of newVars is incorrect.')
    }
  }

  # Check newSamples
  if (is.null(newSamples)) {
    if (is.null(rownames(newBlock))){
      newSamples <- seq(1,nrow(newBlock),1)
    }
  } else {
    if (nrow(newBlock) != length(newSamples)) {
      print(sprintf('Please use a newSamples object with %s items.', nrow(newBlock)))
      stop('The length of newSamples is incorrect.')
    }
  }

  colnames(newBlock) <- NULL
  rownames(newBlock) <- NULL

  # Check Batches
  if (!(is.null(batches))) {
    if (!(is.numeric(batches))){
      print("'batches' should be a numeric vector.")
      stop("The 'batches' object is incorrect.")
    } else {
      if (nrow(newBlock) != length(batches)) {
        print(sprintf("Please use a 'batches' object with %s items.", nrow(newBlock)))
        stop("The length of 'batches' is incorrect.")
      }
    }
  }

  # Check Metadata
  if (!(is.null(metadata))){
    if(is.data.frame(metadata)){
      if(nrow(metadata) != nrow(newBlock)){
        warning("The 'metadata' has an incorrect size and it was subsequently ignored.")
      }
    } else {
      warning("'metadata' must be a data.frame.")
    }

  }

  # Get the block name
  block_name <- as.character(substitute(newBlock))
  if(length(block_name) > 1) {
    block_name <- "X"
  }
  Data <- list(newBlock)
  names(Data) <- block_name

  if(is.null(batches) && is.null(metadata)){
    MB <- MultiBlock(Samples = newSamples,
                     Data = Data,
                     Variables = list(newVars))
  } else if (is.null(batches) && !is.null(metadata)){
    MB <- MultiBlock(Samples = newSamples,
                     Data = Data,
                     Variables = list(newVars),
                     Metadata = list(metadata))
    names(MB@Metadata) <- block_name
  } else {
    MB <- MultiBlock(Samples = newSamples,
                     Data = Data,
                     Variables = list(newVars),
                     Batch = list(batches),
                     Metadata = list(metadata))
    names(MB@Batch) <- names(MB@Metadata) <- block_name
  }

  names(MB@Variables) <- block_name

  return(MB)
}
