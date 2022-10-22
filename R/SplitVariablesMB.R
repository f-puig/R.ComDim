SplitVariablesMB <- function(data = NULL, metadata = NULL, columns = NULL, minblock = 0){

  #' Internal function used for ExpandMultiBlock
  #' Splits a data.frame or a matrix into smaller sets and returns a MultiBlock object.
  #' @param data A data.frame or a matrix. Samples are in rows and variables in columns.
  #' @param metadata The metadata associated with the variables. It has the same number of rows than columns are in data.
  #' @param columns The columns from metadata that will guide the split. Optional
  #' @param minblock The minimal number of variables that need to be contained in a block in order to be kept.
  #' @return The MultiBlock.
  #' @export


  if(is.null(data)){
    stop('No data was provided.')
  }


  if(is.null(columns)){
    if(!is.null(metadata)){
      if(is.null(dim(metadata))){
        metadata <- as.matrix(metadata) # The provided metadata might be a vector.
        columns <- 1
      } else {
        columns <- seq(1,ncol(metadata),1)
      }
    }
  } else {
    if(is.null(dim(metadata))){
      metadata <- as.matrix(metadata) # The provided metadata might be a vector.
    }
  }

  if(is.null(metadata)){
    stop('No metadata was provided.')
  }

  if(nrow(metadata) != ncol(data)){
    cat('The dimensions of metadata and data do not agree. \n')
    stop('The number of rows in metadata should match the number of columns in data.')
  } else {
    rown <- rownames(data)
    data <- t(data) # I work with the transposed matrix because when I subset,
    # if I didn't, when only one element is picked,
    # this will become a vector and produce problems.
    colnames(data) <- rown # Column names were lost after transposition.
    #data <- as.matrix(data)
  }

  list_data <- list()
  list_metadata <- list()
  list_data[[1]] <- data
  list_metadata[[1]] <- metadata
  blocks <- SplitVariables_in(list_data = list_data, list_metadata = list_metadata,
                              columns = columns, it=0, minblock = minblock)

  if(length(blocks) == 0){
    cat('0 blocks were made. Check whether the minblock threshold was set too high.')
  }
  return(blocks)
}


SplitVariables_in <- function(list_data = list_data, list_metadata = list_metadata,
                              columns = columns, it=0, minblock = minblock){
  # Update fields
  if(it != 0){
    aa_data <- list()
    aa_metadata <- list()
    namesi <- names(list_data)
    for(i in 1:length(list_data)){
      namesj <- names(list_data[[i]])
      for(j in 1:length(list_data[[i]])){
        aa_data[[paste(c(namesi[i],namesj[j]), collapse = '.')]] <- list_data[[i]][[j]]
        aa_metadata[[paste(c(namesi[i],namesj[j]), collapse = '.')]] <- list_metadata[[i]][[j]]
      }
    }
    list_metadata <- aa_metadata
    list_data <- aa_data
  }

  growingMB <- NULL
  for(k in 1:length(list_metadata)){ # Number of blocks at the given iteration
    var1 <- unique(list_metadata[[k]][,columns[1]])

    if(length(columns) == 1){ # If the last level is reached
      for(i in 1:length(var1)){
        ind <- which(list_metadata[[k]][,columns[1]] == var1[i])
        if(minblock <= length(ind)){
          xx <- list_data[[k]][ind,]
          if(length(ind) == 1){
            xx <- as.matrix(xx)
          } else {
            xx <- t(xx)
          }
          namesk <- names(list_data)[k]
          namesf <- paste(c(namesk,var1[[i]]), collapse = '.')

          if(is.null(growingMB)){
            growingMB <- BuildMultiBlock(xx)
            names(growingMB@Data) <- namesf
            names(growingMB@Variables) <- namesf
          } else {
            growingMB2 <- BuildMultiBlock(xx)
            names(growingMB2@Data) <- namesf
            names(growingMB2@Variables) <- namesf
            growingMB <- BuildMultiBlock(growingMB, growingMB2)
          }
        }
        if(i == length(var1) && k == length(list_data)){
          return(growingMB)
        }
      }
    } else {
      for(i in 1:length(var1)){
        if(i == 1){
          aa <- list()
          aa_data <- list()
        }
        ind <- which(list_metadata[[k]][,columns[1]] == var1[i])
        aa[[var1[[i]]]] <- list_metadata[[k]][ind,]
        aa_data[[var1[i]]] <- list_data[[k]][ind,]
        if(i == length(var1)){
          list_metadata[[k]] <- aa
          list_data[[k]] <- aa_data
        }
      }
      if(k == length(list_metadata)){
        return(SplitVariables_in(list_data = list_data,
                             list_metadata = list_metadata,
                             columns = columns[-1], it = it+1,
                             minblock = minblock))
      }
    }
  }
}

