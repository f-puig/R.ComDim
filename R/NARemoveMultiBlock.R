# This function requires the package imputeLCMD and zoo.
NARemoveMultiBlock <- function(MB = MB,
                               blocks = NULL,
                               minfrac = 0.5,
                               method = c('none', 'zero', 'median', 'discard', 'fixed.value',
                                          'fixed.value.all', 'fixed.noise',
                                          'random.noise', 'QRILC')[8],
                               constant = 0,
                               factor.NA = 0.5,
                               sd.noise = 0.3,
                               seed.number = NULL,
                               tune.sigma = 1,
                               showWarning = TRUE)
                               {

  #' NARemoveMultiBlock
  #'
  #' Remove NA from a MultiBlock object.
  #' @param MB The MultiBlock object.
  #' @param blocks The blocks to apply the NA imputation. It can be a vector of integers or a vector with the block names (Optional).
  #' @param minfrac minimum fraction of samples with elements different than NA and Infinite necessary. If this number is not reached for a variable, this variable will be discarded. The default value is 0.5.
  #' @param method The method to use for imputation of the NA values. 'none' will not do any transformation besides the variable discard based on the minfrac. 'zero' will convert all NA to 0 values. 'median' will replace the NA values with the median value of the column. 'fixed.value' will replace the NA by a constant number given by the user. 'fixed.value.all' will add a constant value to all the elements (NAs and not NAs) in the multi-block.  'discard' will delete the columns containing NA values. 'fixed.noise' and 'random.noise' will replace the NA values with noise. 'QRILC' will impute values using the Quantile Regression Imputation Left-Censored data method.
  #' @param constant For the 'fixed.value' method, the value by which the NAs are replaced. For the 'fixed.value.all', the value added to each element in the multi-block. The default number is 0.
  #' @param factor.NA For the 'fixed.noise' and the 'random.noise' method. For the 'fixed.noise' approach, this is the factor by which the minimal non-NA value of the column is multiplied. For the 'random.value' approach, this value will be multiplied to the minimal value of each column to determine the mean of the normal distribution. The default number is 0.5.
  #' @param sd.noise For the 'random.noise' method, this is the factor used to define the SD of the normal distribution of the random noise. The SD will be equal to sd.noise multiplied by the mean of the random noise. The default number is 0.3.
  #' @param seed.number For the 'random.noise' method, this is the seed to create the random numbers (Optional).
  #' @param tune.sigma For the 'QRILC' method. A scalar used to tune the standard deviation (if the complete data distribution is not Gaussian). The default value is tune.sigma = 1, and it corresponds to the case where the complete data distribution is Gaussian.
  #' @param showWarning If TRUE, it will return a warning in case there is any variable with only NAs in the multi-block.
  #' @return The MultiBlock
  #' @examples
  #' b1 = matrix(rnorm(500),10,50)
  #' b2 = matrix(rnorm(800),10,80)
  #' b2[c(2,3,5),c(1,2,3)] <- NA
  #' # Build multi-block by adding one data block at a time:
  #' mb <- BuildMultiBlock(b1, b2)
  #' mb <- NARemoveMultiBlock(mb, method = 'zero')
  #' @export


  if(is.null(method) || length(method) > 1){
    stop("A method should be provided. Try 'none', 'zero', 'median', 'discard', 'fixed.value', 'fixed.noise', 'random.noise' or 'QRILC'.")
  }

  if(!("MultiBlock" %in% class(MB))){
    stop("'MB' is not a MultiBlock object.")
  }

  method <- tolower(method)
  if(!(method %in% c('none', 'zero', 'median', 'discard', 'fixed.noise', 'fixed.value', 'fixed.value.all', 'random.noise', 'qrilc'))){
    stop("The provided method is not in the list of available methods. Try 'zero', 'median', 'discard', 'fixed.value', 'fixed.value.all', 'fixed.noise', 'random.noise' or 'QRILC'.")
  }

  if(!(is.null(blocks))){
    if(is.character(blocks)){
      block.ids <- match(blocks, names(MB@Data))
      if(any(is.na(block.ids))){
        warning('Some block names were not found in the MultiBlock: %s',
                paste0(blocks[which(is.na(block.ids))], sep=', '))
        block.ids <- block.ids[!is.na(block.ids)]
      }
    } else if(is.numeric(blocks)){
      block.ids <- intersect(1:length(MB@Data), blocks)
      if(length(blocks) != length(block.ids)){
        warning('Some block names were not found in the MultiBlock: %s',
                paste0(setdiff(1:length(MB@Data), blocks), sep=', '))
      }
    }
  } else {
    block.ids <- seq(1,length(MB@Data))
  }
  if(length(block.ids) == 0){
    warning('0 blocks were selected. Data was not modified.')
    return(MB)
  }

  if(is.null(minfrac)){
    if(showWarning){
      warning("'minfrac' was set to 0.") # The variable will be kept even if all elements are NAs and Infinites.
    }
    minfrac <- 0
  } else if (minfrac < 0){
    if(showWarning){
      warning("'minfrac' was set to 0.") # The variable will be kept even if all elements are NAs and Infinites.
    }
    minfrac <- 0
  } else if (minfrac > 1){
    if(showWarning){
      warning("'minfrac' was set to 1.") # If there is a NA or a Inf in a variable, this will be deleted.
    }
    minfrac <- 1
  }
  minfrac <- round(nrow(MB@Data[[1]])*minfrac)

  to_ignore <- 0
  to_delete <- vector() # Blocks that will be deleted because only contained NA.
  for(i in block.ids){

    if(any(is.complex(MB@Data[[i]]))){
      stop('The multi-block contain complex values that cannot be processed with this function.')
    }

    if(any(is.logical(MB@Data[[i]]))){
      stop('The multi-block contain logical values that cannot be processed with this function.')
    }

    numna <- apply(MB@Data[[i]], 2, function(x) sum(is.na(x)))
    colna <- which(numna >= minfrac)
    if(length(colna) > 0){
      MB@Data[[i]] <- MB@Data[[i]][,setdiff(1:length(MB@Variables[[i]]),colna)]
      MB@Variables[[i]] <- MB@Variables[[i]][setdiff(1:length(MB@Variables[[i]]),colna)]
      to_ignore <- to_ignore + length(colna)
    }
    if(length(MB@Variables[[i]]) == 0){
      to_delete <- append(to_delete, i)
    }
  }
  if(length(to_delete) != 0){
    to_delete <- rev(to_delete)
    names_to_delete <- vector()
    old_names <- names(MB@Data)
    for(i in to_delete){
      names_to_delete <- append(names_to_delete, names(MB@Data)[i])
      if(names(MB@Data)[i] %in% names(MB@Batch)){
        MB@Batch[[ names(MB@Data)[i] ]] <- NULL
      }
      if(names(MB@Data)[i] %in% names(MB@Metadata)){
        MB@Metadata[[ names(MB@Data)[i] ]] <- NULL
      }
      MB@Data[[i]] <- NULL
      MB@Variables[[i]] <- NULL
    }
    kept_names <- setdiff(old_names[block.ids],names_to_delete)
    if(length(kept_names) != 0){
      block.ids <- match(kept_names,names(MB@Data))
    } else {
      stop('The selected blocks contained only NA values. They cannot be used in ComDim.')
    }
  }

  for(i in block.ids){

    if(method == 'none'){

      next # This method can be interesting to use if we only want to apply the minfrac criteria.

    } else if(method == 'zeros'){

      pos <- which(is.na(MB@Data[[i]]))
      MB@Data[[i]][pos] <- rep(0,length(pos))

    } else if(method == 'median'){
      for(nc in 1:ncol(MB@Data[[i]])){
        if(any(is.na(MB@Data[[i]][,nc]))){
          pos <- which(is.na(MB@Data[[i]][,nc]))
          MB@Data[[i]][pos,nc] <- rep(median(MB@Data[[i]][,nc], na.rm = TRUE), length(pos))
        }
      }
    } else if(method == 'fixed.value'){
      for(nc in 1:ncol(MB@Data[[i]])){
        if(any(is.na(MB@Data[[i]][,nc]))){
          pos <- which(is.na(MB@Data[[i]][,nc]))
          MB@Data[[i]][pos,nc] <- rep(constant, length(pos))
        }
      }
    } else if(method == 'fixed.value.all'){
      MB@Data[[i]] <- MB@Data[[i]] + constant # This will add constant to all cells except where there is a NA.
      for(nc in 1:ncol(MB@Data[[i]])){
        if(any(is.na(MB@Data[[i]][,nc]))){
          pos <- which(is.na(MB@Data[[i]][,nc]))
          MB@Data[[i]][pos,nc] <- rep(constant, length(pos)) # And here we fill the positions with NA.
        }
      }
    } else if(method == 'discard'){
      pos <- vector()
      for(nc in 1:ncol(MB@Data[[i]])){
        if(any(is.na(MB@Data[[i]][,nc]))){
          pos <- append(pos, nc)
        }
      }
      if(length(pos) == ncol(MB@Data[[i]])){
        MB[[i]] <- NULL
        warning(sprint("The block %s was removed because only contained NA values.", names(MB[i])))
      } else if(length(pos) != 0){
        MB@Data[[i]] <- MB@Data[[i]][,setdiff(1:ncol(MB@Data[[i]]),pos)]
        MB@Variables[[i]] <- MB@Variables[[i]][setdiff(1:ncol(MB@Data[[i]]),pos)]
      }
    } else if(method == 'fixed.noise'){
      for(nc in 1:ncol(MB@Data[[i]])){
        if(any(is.na(MB@Data[[i]]))){
          minval <- min(MB@Data[[i]][,nc], na.rm = TRUE)
          pos <- which(is.na(MB@Data[[i]][,nc]))
          MB@Data[[i]][pos,nc] <- rep(minval*factor.NA, length(pos))
        }
      }
    } else if(method == 'random.noise'){

      if(!is.null(seed.number)){
        old <- .Random.seed
        on.exit({.Random.seed <<- old}) # Seed will be restaured on exit.
        set.seed(seed.number)
      }

      for(nc in 1:ncol(MB@Data[[i]])){
        if(any(is.na(MB@Data[[i]][,nc]))){
          minval <- min(MB@Data[[i]][,nc], na.rm = TRUE)
          pos <- which(is.na(MB@Data[[i]][,nc]))
          MB@Data[[i]][pos,nc] <- rnorm(length(pos),
                                        mean = minval*factor.NA,
                                        sd = minval*factor.NA*sd.noise)
        }
      }
    } else if(method == 'QRILC'){
      if(any(is.na(MB@Data[[i]]))){
        for(nc in 1:ncol(MB@Data[[i]])){
          if(any(is.na(MB@Data[[i]][,nc]))){
            MB@Data[[i]][,nc] <- imputeLCMD::impute.QRILC(
              as.matrix(MB@Data[[i]][,nc]), tune.sigma = tune.sigma)[[1]]
          }
        }
      }
    }
  }

  if(to_ignore != 0 && showWarning){
    warning(sprintf('%s variable(s) only contained NAs and were deleted.', to_ignore))
  }

  return(MB)
}
