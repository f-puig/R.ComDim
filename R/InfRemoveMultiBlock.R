InfRemoveMultiBlock <- function(MB = MB,
                               blocks = NULL,
                               minfrac = 0.5,
                               method = c('none', 'fixed.noise', 'random.noise')[3],
                               factor.NA = 0.5,
                               sd.noise = 0.3,
                               seed.number = NULL,
                               showWarning = TRUE)
                               {

  #' InfRemoveMultiBlock
  #'
  #' Remove NA from a MultiBlock structure.
  #' @param MB The MultiBlock structure.
  #' @param blocks The blocks to apply the NA imputation. It can be a vector of integers or a vector with the block names. Facultative.
  #' @param minfrac minimum fraction of samples with elements different than NA and Infinite necessary. If this number is not reached for a variable, this variable will be discarded. The default value is 0.5.
  #' @param method The method to use for imputation of the Inf values. 'none' will not do any transformation besides the variable discard based on the minfrac. 'fixed.noise' and 'random.noise' will replace the Inf values with noise, estimated from the minimum and maximum non-NA and non-inf values in each column. 'random.noise' differs from 'fixed.noise' in that the value to replace the Inf in 'fixed.noise' is constant while for 'random.noise' it is not. 'random.noise' obtain these values from two normal distributions (to replace the Inf and the -Inf, respectively).
  #' @param factor.NA Used in the 'fixed.noise' and the 'random.noise' method. For the 'fixed.noise' approach, -Inf will be replaced by the minimal non-NA value of the column multiplied by this number. Positive infinite values will be replaced by the max value multiplied by the inverse of this number. For the 'random.value' approach, these two numbers will be used as the mean of the normal distributions to impute the negative and the positive infinite values, respectively. The default number is 0.5.
  #' @param sd.noise For the 'random.noise' method, this is the factor used to define the SD of the normal distribution of the random noise. The SD will be equal to sd.noise multiplied by the mean of the random noise. The default number is 0.3.
  #' @param seed.number For the 'random.noise' method, this is the seed to create the random numbers. Facultative.
  #' @param showWarning If TRUE, it will return a warning in case there is any variable with only NAs in the multi-block.
  #' @return The multi-block
  #' @examples
  #' b1 = matrix(rnorm(500),10,50)
  #' b2 = matrix(rnorm(800),10,80)
  #' b2[c(2,3,5),c(1,2,3)] <- Inf
  #' # Build multi-block by adding one data block at a time:
  #' mb <- BuildMultiBlock(b1, newSamples = paste0('sample_',1:10))
  #' mb <- BuildMultiBlock(b2, growingMB = mb)
  #' mb <- InfRemoveMultiBlock(mb, method = 'zero')
  #' @export


  if(is.null(method) || length(method) > 1){
    stop("A method should be provided. Try 'none', 'fixed.noise' or 'random.noise'.")
  }

  if(!("MultiBlock" %in% class(MB))){
    stop("'MB' is not a MultiBlock object.")
  }

  method <- tolower(method)
  if(!(method %in% c('none', 'fixed.noise', 'random.noise'))){
    stop("The provided method is not in the list of available methods. Try 'none', 'fixed.noise' or 'random.noise'.")
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
  } else if (minfrac > 0){
    if(showWarning){
      warning("'minfrac' was set to 1.") # If there is a NA or a Inf in a variable, this will be deleted.
    }
    minfrac <- 1
  }
  minfrac <- round(nrow(MB@Data[[1]])*minfrac)

  to_ignore <- 0
  to_delete <- vector() # Blocks that will be deleted because only contained NA or Inf.
  for(i in block.ids){

    if(any(is.complex(MB@Data[[i]]))){
      stop('The multi-block contain complex values that cannot be processed with this function.')
    }

    if(any(is.logical(MB@Data[[i]]))){
      stop('The multi-block contain logical values that cannot be processed with this function.')
    }

    # Delete a column if only contains NAs and Infinite.
    numna <- apply(MB@Data[[i]], 2, function(x) sum(c(is.na(x),is.infinite(x))))
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
    old_names <- names(MB)
    for(i in to_delete){
      names_to_delete <- append(names_to_delete, names(MB[i]))
      MB[[i]] <- NULL
    }
    kept_names <- setdiff(old_names[block.ids],names_to_delete)
    if(length(kept_names) != 0){
      block.ids <- match(kept_names,names(MB))
    } else {
      stop('The selected blocks contained only NA or Infinite values. They cannot be used in ComDim.')
    }
  }


  to_ignore <- 0
  to_delete <- vector() # Blocks that will be deleted because only contained NA.
  for(i in block.ids){

    if(any(is.complex(MB@Data[[i]]))){
      stop('The multi-block contain complex values that cannot be processed with this function.')
    }

    if(any(is.logical(MB@Data[[i]]))){
      stop('The multi-block contain logical values that cannot be processed with this function.')
    }

    numna <- apply(MB@Data[[i]], 2, function(x) sum(c(is.na(x),is.infinite(x))))
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
      stop('The selected blocks contained only Infinite or NA values. They cannot be used in ComDim.')
    }
  }


  for(i in block.ids){

    if(method == 'none'){

      next # This method can be interesting to use if we only want to apply the minfrac criteria.

    } else if(method == 'fixed.noise'){
      for(nc in 1:ncol(MB@Data[[i]])){
        if(any(is.infinite(MB@Data[[i]]))){
          minval <- min(MB@Data[[i]][,nc], na.rm = TRUE)
          maxval <- max(MB@Data[[i]][,nc], na.rm = TRUE)
          pos <- which(is.infinite(MB@Data[[i]][,nc] & MB@Data[[i]][,nc] > 0))
          if(length(pos) > 0) {
            MB@Data[[i]][pos,nc] <- rep(minval*(1/factor.NA), length(pos))
          }
          neg <- which(is.infinite(MB@Data[[i]][,nc] & MB@Data[[i]][,nc] < 0))
          if(length(neg) > 0){
            MB@Data[[i]][neg,nc] <- rep(minval*factor.NA, length(neg))
          }
        }
      }
    } else if(method == 'random.noise'){

      if(!is.null(seed.number)){
        old <- .Random.seed
        on.exit({.Random.seed <<- old}) # Seed will be restaured on exit.
        set.seed(seed.number)
      }

      for(nc in 1:ncol(MB@Data[[i]])){
        if(any(is.infinite(MB@Data[[i]][,nc]))){
          minval <- min(MB@Data[[i]][,nc], na.rm = TRUE)
          maxval <- max(MB@Data[[i]][,nc], na.rm = TRUE)
          pos <- which(is.infinite(MB@Data[[i]][,nc] & MB@Data[[i]][,nc] > 0))
          neg <- which(is.infinite(MB@Data[[i]][,nc] & MB@Data[[i]][,nc] < 0))
          if(length(neg) > 0){
            MB@Data[[i]][neg,nc] <- rnorm(length(neg),
                                          mean = minval*factor.NA,
                                          sd = minval*factor.NA*sd.noise)
          }
          if(length(pos) > 0){
            MB@Data[[i]][pos,nc] <- rnorm(length(pos),
                                          mean = minval*(1/factor.NA),
                                          sd = minval*(1/factor.NA)*sd.noise)
          }
        }
      }
    }
  }

  if(to_ignore != 0 && showWarning){
    warning(sprintf('%s variable(s) only contained Infinite and NAs. They were deleted.', to_ignore))
  }

  return(MB)
}
