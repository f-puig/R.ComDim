# Pareto, mean-center, auto-scale, geometric mean, ranknorm

NormalizeMultiBlock <- function(MB = MB,
                               blocks = NULL,
                               method = c('none', 'auto', 'mean', 'pareto',
                                          'norm', 'geometric', 'ranknorm')[5],
                               infinite.as.NA = FALSE,
                               constant = 0,
                               offset = 3/8,
                               showWarning = TRUE)
                               {

  #' NARemoveMultiBlock
  #'
  #' Normalize all blocks from a multi-block structure. The ranknorm transform is based on that from the RNOmni package (genomedoi:10.1111/biom.13214).
  #' @param MB The multi-blocks structure.
  #' @param blocks The blocks to apply the normalization. It can be a vector of integers or a vector with the block names. Facultative.
  #' @param method The method to use for data normalization. 'none' will not do any transformation besides replacing the Inf values to NA if infinite.as.NA is set to TRUE. 'auto' stands for auto-scaling, 'mean' for mean-centering, 'pareto' for Pareto transform (mean-centered and divided by the square root of the standard deviation), 'norm' for mean-centered and divided by its norm, 'geometric' for geometric mean transform, and 'ranknorm' for the rank-based inverse normal transform.
  #' @param infinite.as.NA If TRUE, Infinite values are converted to NA. The presence of NA or Inf values will cause ComDim to stop, but they can be transformed to numerical values with the functions 'NARemoveMultiBlock' and 'InfRemoveMultiBlock', respectively.
  #' @param constant For the 'geometric' method, the value added to each element in the multi-block. The default number is 0.
  #' @param offset For the 'ranknorm' method. Defaults to (3/8), correspond to the Blom transform.
  #' @param showWarning If TRUE, it will return a warning in case there is any variable with only NAs in the multi-block.
  #' @return The multi-block
  #' @examples
  #' b1 = matrix(rnorm(500),10,50)
  #' b2 = matrix(rnorm(800),10,80)
  #' b2[c(2,3,5),c(1,2,3)] <- NA
  #' # Build multi-block by adding one data block at a time:
  #' mb <- BuildMultiBlock(b1, newSamples = paste0('sample_',1:10))
  #' mb <- BuildMultiBlock(b2, growingMB = mb)
  #' mb <- NormalizeMultiBlock(mb, method = 'auto')
  #' @export


  if(is.null(method) || length(method) > 1){
    stop("A method should be provided. Try 'none', 'auto', 'mean', 'pareto', 'norm', 'geometric', or 'ranknorm'.")
  }

  if(!("MultiBlock" %in% class(MB))){
    stop("'MB' is not a MultiBlock object.")
  }

  method <- tolower(method)
  if(!(method %in% c('none','auto', 'mean', 'pareto', 'norm', 'geometric', 'ranknorm'))){
    stop("The provided method is not in the list of available methods. Try 'none', 'auto', 'mean', 'norm', 'pareto', 'geometric or 'ranknorm'.")
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
      block.ids <- intersect(1:length(MB), blocks)
      if(length(blocks) != length(block.ids)){
        warning('Some block names were not found in the MultiBlock: %s',
                paste0(setdiff(1:length(MB), blocks), sep=', '))
      }
    }
  } else {
    block.ids <- seq(1,length(MB@Data))
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

    if(any(is.infinite(MB@Data[[i]]))){
      if(infinite.as.NA){
        MB@Data[[i]][is.infinite(MB@Data[[i]])] <- NA
      } else{
        cat('The multi-block contain infinite values that cannot be processed.\n')
        warning(sprintf("You can convert them to NA with infinite.as.NA = TRUE or using the function 'InfRemoveMultiBlock'."))
      }

    }

    # Delete rows containing only NA
    numna <- apply(MB@Data[[i]], 2, function(x) sum(is.na(x)))
    colna <- which(numna >= nrow(MB@Data[[i]]))
    if(length(colna) > 0){
      MB@Data[[i]] <- MB@Data[[i]][,setdiff(1:length(MB@Variables[[i]]),colna)]
      MB@Variables[[i]] <- MB@Variables[[i]][setdiff(1:length(MB@Variables[[i]]),colna)]
      to_ignore <- to_ignore + length(colna)
    }
    if(length(MB@Variables[[i]]) == 0){
      to_delete <- append(to_delete, i)
    }
  }
  # Delete blocks containing only NA
  if(length(to_delete) != 0){
    to_delete <- rev(to_delete)
    names_to_delete <- vector()
    old_names <- names(MB@Data)
    for(i in to_delete){
      names_to_delete <- append(names_to_delete, names(MB@Data)[i])
      MB@Data[[ i ]] <- NULL
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
      next # This method can be interesting to use if we only want to transform the Inf values to NA.
    }
    else if(method == 'auto'){

      MB@Data[[i]] <- scale(MB@Data[[i]], center = TRUE, scale = TRUE)

    } else if(method == 'mean'){

      MB@Data[[i]] <- scale(MB@Data[[i]], center = TRUE, scale = FALSE)

    } else if(method == 'pareto'){

      MB@Data[[i]] <- scale(MB@Data[[i]], center = TRUE,
                            scale = apply(MB@Data[[i]], 2, function(x) sqrt(sd(x, na.rm = TRUE)) ))

    } else if(method == 'norm'){

      Xmean <- scale(MB@Data[[i]], center = TRUE, scale = FALSE)
      Norm_X <- sqrt(sum(Xmean*Xmean))
      MB@Data[[i]] <- Xmean/Norm_X

    }else if(method == 'geometric'){

      for(nc in 1:length(MB@Variables[[i]])){
        notna <- which(!is.na(MB@Data[[i]][,nc]))
        if(length(notna) == 0) { next}
        if(any(MB@Data[[i]][notna,nc] < 0)){
          stop('The multi-block contains negative values that cannot be normalized with the geometric mean.')
        }
        MB@Data[[i]][notna,nc] <- prod(MB@Data[[i]][notna,nc])^(1/length(MB@Data[[i]][notna,nc]))
      }

    } else if(method == 'ranknorm'){

      for(nc in 1:length(MB@Variables[[i]])){
        notna <- which(!is.na(MB@Data[[i]][,nc]))
        n <- length(notna)
        r <- rank(MB@Data[[i]][notna,nc])
        MB@Data[[i]][notna,nc] <- qnorm((r - offset)/(n - 2 * offset + 1))
      }
    }
  }

  if(to_ignore != 0 && showWarning){
    warning(sprintf('%s variable(s) only contained NAs and were deleted.', to_ignore))
  }

  return(MB)
}
