MultiBlock2matrix <- function(MB = MB,
                             blocks = NULL,
                             vars = NULL
                             )
                             {

  #' MultiBlock2matrix
  #'
  #' Combines the blocks of a MultiBlock into a matrix. Useful for calculating summary information of the whole MB (ex. maximum value).
  #' @param MB The MultiBlocks structure.
  #' @param blocks The blocks to combine into a matrix. It can be a vector of integers or a vector with the block names.
  #' @param vars The variables to kept. It is a list object with the same length as 'blocks'. Each element in the list contains the position of the variables in the nth block.
  #' @return The MultiBlock
  #' @examples
  #' b1 = matrix(rnorm(500),10,50)
  #' b2 = matrix(rnorm(800),10,80)
  #' # Build MultiBlock by adding one data block at a time:0))
  #' mb <- BuildMultiBlock(b1, newSamples = paste0('sample_',1:1
  #' mb <- BuildMultiBlock(b2, growingMB = mb)
  #' Convert data to a matrix.
  #' mb <- MultiBlock2matrix(MB)
  #' max(mb)

  #' @export

  if(!("MultiBlock" %in% class(MB))){
    stop("'MB' is not a MultiBlock object.")
  }

  if(!(is.null(blocks))){
    if(is.character(blocks)){
      block.ids <- match(blocks, names(MB@Data))
      if(all(is.na(block.ids))){
        stop('None of the block names were found in the MultiBlock.')
      }
      if(any(is.na(block.ids))){
        warning('Some block names were not found in the MultiBlock: %s',
                paste0(blocks[which(is.na(block.ids))], sep=', '))
        block.ids <- block.ids[!is.na(block.ids)]
      }
    } else if(is.numeric(blocks)){
      block.ids <- intersect(1:length(MB@Data), blocks)
      if(length(block.ids) == 0){
        stop('None of the block names were found in the MultiBlock.')
      }
      if(length(blocks) != length(block.ids)){
        warning('Some block names were not found in the MultiBlock: %s',
                paste0(setdiff(1:length(MB@Data), blocks), sep=', '))
      }
    }
  } else {
    block.ids <- seq(1,length(MB@Data))
  }

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
          warning(sprintf('Some variables were not found in the MultiBlock object %s: %s', i,
                  paste0(setdiff(vars[[i]], 1:length(MB@Variables[[i]])), sep=', ')))
        }
        vars[[i]] <- xx[!is.na(xx)]
      }
    }
  } else {
    vars <- list()
    for(i in block.ids){
      vars[[i]] <- 1:length(MB@Variables[[i]])
    }
  }

  cnames <- vector()
  rep_names <- vector()
  block_names <- vector()
  common <- FALSE

  for(i in block.ids){
    if(i == block.ids[1]){
      mat <- MB@Data[[i]][,vars[[i]]]
    } else {
      mat <- cbind(mat, MB@Data[[i]][,vars[[i]]])
    }

    names_new <- MB@Variables[[i]]
    if(any(names_new %in% cnames)){
      common <- TRUE
      rep_names <- append(rep_names, names_new[which(names_new %in% cnames)])
      names_new[which(names_new %in% cnames)] <-
        paste0(names_new[which(names_new %in% cnames)],
               sprintf(".%s", names(MB@Data)[i]))
    }
    block_names <- append(block_names, rep(names(MB@Data)[i], length(names_new)))
    cnames <- append(cnames, names_new)
  }

  # Now replace the name of the variables that had a replicated name.
  if(length(rep_names) != 0){
    unique_rep_names <- as.character(unique(rep_names))
    cnames <- as.character(cnames)
    for(i in unique_rep_names){
      pos <- which(cnames == i)
      cnames[pos] <- paste0(i,".",block_names[match(i, cnames)])
    }
  }

  if(common){
    warning("Some variable names were repeated across blocks.
  The repeated variable names were renamed.")
  }

  rownames(mat) <- MB@Samples
  colnames(mat) <- cnames

  return(mat)
}
