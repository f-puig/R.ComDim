ProcessMultiBlock <- function(MB = MB,
                              blocks = NULL,
                              vars = NULL,
                              FUN = NULL,
                              FUN.SelectVars = NULL,
                              FUN.SelectBlocks = NULL
                              )
                              {

  #' ProcessMultiBlock
  #'
  #' Apply a custom function to Process the MultiBlocks and/or to select some blocks. If multiple functions are used, the order followed is FUN, FUN.SelectVars, and FUN.SelecBlocks.
  #' @param MB The MultiBlocks structure.
  #' @param blocks The blocks to apply the processing. It can be a vector of integers or a vector with the block names. Facultative.
  #' @param vars The variables to kept. It is a list object with the same length as 'blocks'. Each element in the list contains the position of the variables in the nth block.
  #' @param FUN The function used for data processing.
  #' @param FUN.SelectBlocks The function used to select blocks. If applied to a numeric data matrix, it should return TRUE/FALSE.
  #' @param FUN.SelectVars The function used to select variables. If applied to a numeric data matrix, it should return a vector of TRUE/FALSE with the same length as ncol.
  #' @return The MultiBlock
  #' @examples
  #' b1 = matrix(rnorm(500),10,50)
  #' b2 = matrix(rnorm(800),10,80)
  #' # Build MultiBlock by adding one data block at a time:0))
  #' mb <- BuildMultiBlock(b1, newSamples = paste0('sample_',1:1
  #' mb <- BuildMultiBlock(b2, growingMB = mb)
  #' # Prepare custom functions
  #' BY100 <- function(x) { 100*(x/max(x)) } # Normalize each block to 100%
  #' norm100000 <- function(x){norm(x)>100000} # Keep blocks with a norm higher than 100000.
  #' thr5 <- function(x) {x > max(x, na.rm = TRUE) * 0.05} # Keep variables larger than 5% in the block.
  #' mb <- ProcessMultiBlock(MB, FUN = BY100, FUN.SelectBlocks = norm100000, FUN.SelectVars = thr5)

  #' @export


  if(is.null(FUN) && is.null(FUN.SelectBlocks) &&
     is.null(FUN.SelectVars) && is.null(blocks) &&
     is.null(vars)){
     warning("'FUN','FUN.SelectBlocks' and 'FUN.SelectVars' arguments were missing. 'MB' was not processed.")
     return(MB)
  }

  if(!("MultiBlock" %in% class(MB))){
    stop("'MB' is not a MultiBlock object.")
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
  block.fixed <- setdiff(seq(1,length(MB@Data)), block.ids)

  if(!is.null(vars)){

    if(!is.list(vars)){
      stop("'variables' must be a list object.")
    } else if(length(vars) != length(block.ids)){
      stop("'variables' should have the same length as 'blocks'.")
    }

    for(i in vars){
      if(is.character(vars[[i]])){
        xx <- match(vars[[i]], MB@Variables[[i]])
        if(any(is.na(vars[[i]]))){
          warning(sprintf('Some variable names were not found in block %s: %s', i,
                  paste0(vars[[i]][which(is.na(xx))], sep=', ')))
          vars[[i]] <- xx[!is.na(xx)]
        }
      } else if(is.numeric(vars[[i]])){
        xx <- intersect(1:length(MB@Variables[[i]]), vars[[i]])
        if(length(xx) != length(MB@Variables[[i]])){
          warning(sprintf('Some variables were not found in the MultiBlock %s: %s', i,
                  paste0(setdiff(vars[[i]], 1:length(MB@Variables[[i]])), sep=', ')))
        }
        vars[[i]] <- xx[!is.na(xx)]
      }
    }
  }

  if(!is.null(FUN)){

    for(i in block.ids){

      MB@Data[[i]] <- FUN(MB@Data[[i]])

    }

  }

  if(!is.null(FUN.SelectVars)){

    for(i in block.ids){

      pos <- which(FUN.SelectVars(MB@Data[[i]]))

      if(length(pos) != 0){
        MB@Data[[i]] <- MB@Data[[i]][,pos]
        MB@Variables[[i]] <- MB@Variables[[i]][pos]
      }

    }

  }

  if(!is.null(FUN.SelectBlocks)){

    block.kept <- vector()
    for(i in block.ids){

     if(!FUN.SelectBlocks(MB@Data[[i]])) {
       block.kept <- append(block.kept, block.ids[i])
     }

    }

    block.rm <- setdiff(block.ids, block.kept)
    if(length(block.rm) != 0){
      block.rm <- rev(block.rm)
      for(i in block.rm){
        MB@Data[[i]] <- NULL
        MB@Variables[[i]] <- NULL
      }
    }
    #block.ids <- block.kept

  }

  return(MB)
}
