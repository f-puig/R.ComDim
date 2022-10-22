nrowMultiBlock <- function(MB, block) {

  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Return the number of rows in each block of the MultiBlock.
  #' @param MB A MultiBlock.
  #' @param block The blocks to calculate the number of rows
  #' @return A vector with the number of rows
  #' @export
  #'

  if(any(class(MB) == "MultiBlock")){
    rows <- vector()

    if(missing(block)){
      block <- 1:length(MB@Data)
    }

    for(i in 1:length(block)){
      rows[i] <- nrow(MB@Data[[ block[i] ]])
    }

    return(rows)
  }
}
