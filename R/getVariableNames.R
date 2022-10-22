getVariableNames <- function(MB, block) {

  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Get the variable names.
  #' @param MB A MultiBlock object.
  #' @param block A vector with the list of blocks.
  #' @return The vector with the variable names.
  #' @export
  #'

  if(any(class(MB) != "MultiBlock")){
    stop("'MB' is not a MultiBlock.")
  }

  if(missing(block)){
    block <- 1:length(MB@Data)
  }

  if(is.list(block)){
    block <- unlist(block)
  }


  if(length(block) == 1){
    thevars <- MB@Variables[[ block ]]
  } else {
    thevars <- list()
    if(is.character(block)){
      for(i in block){
        thevars[[ i ]] <- MB@Variables[[ i ]]
      }
    } else {
      blocknames <- getBlockNames(MB)
      for(i in block){
        thevars[[ blocknames[i] ]] <- MB@Variables[[ blocknames[i] ]]
      }
    }

  }

  return(thevars)
}


