ncolMultiBlock <- function(MB, block) {

  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Return the number of columns in each block of the MultiBlock.
  #' @param MB A MultiBlock.
  #' @param block The blocks to calculate the number of columns.
  #' @return A vector with the number of columns.
  #' @export
  #'

  if(any(class(MB) == "MultiBlock")){
    cols <- vector()

    if(missing(block)){
      block <- 1:length(MB@Data)
    }

    for(i in 1:length(block)){
      cols[i] <- ncol(MB@Data[[ block[i] ]])
    }

    return(cols)
  }
}
