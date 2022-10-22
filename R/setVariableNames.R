setVariableNames <- function(MB, newVars, block) {

  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Get the variable names.
  #' @param MB A MultiBlock object.
  #' @param newVars A vector with the new variable names.
  #' @param block A vector with the list of blocks.
  #' @return The MultiBlock object
  #' @export
  #'

  if(any(class(MB) != "MultiBlock")){
    stop("'MB' is not a MultiBlock.")
  }

  if(missing(block)){
    if(length(newVars) == length(MB@Data)){
      block <- 1:length(MB@Data)
    } else {
      stop('The blocks to replace the variable names must be indicated.')
    }
  }

  if(is.list(block)) {
    block <- unlist(block) # Convert to vector
  }

  if(length(block) > 1 & is.vector(newVars)){
    stop('To replace variable names from more than 2 blocks, newVars must be a list.')
  }

  if(is.list(newVars)){
    if(length(newVars) != length(block)){
      stop('newVars and block should have the same length.')
    }
    for(i in 1:length(newVars)){
      if(ncol(MB@Data[[ block[i] ]]) == length(newVars[[i]])){
        MB@Variables[[ block[i] ]] <- newVars[[i]]
      } else {
        stop('The length of newVars is not correct for at least 1 of the blocks.')
      }
    }
  } else {
    if(ncol(MB@Data[[ block ]]) == length(newVars)){
      MB@Variables[[ block ]] <- newVars
    } else {
      stop('The length of newVars is not correct for at least 1 of the blocks.')
    }
  }

  return(MB)
}


