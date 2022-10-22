MakeComDimLoadingsTable <- function(model,
                            blocks = NULL,
                            dim = NULL,
                            dim.ort = NULL) {

#' MakeComDimLoadingsTable
#'
#' Creates a long dataframe with the global and local scores (Q.scores and T.scores) suitable for ggplot.
#' @param model The output from a ComDim analysis.
#' @param blocks The blocks from which the data will be extracted. It can be a vector of integers or a vector with the block names. Facultative.
#' @param dim The dimensions from which the data will be extracted. It can be a vector of integers. Facultative.
#' @param dim.ort The orthogonal dimensions from which the data will be extracted. It can be a vector of integers. If NULL, all orthogonal components are extracted. If 0, none of the orthogonal components are extracted. Facultative.
#' @return A table with the selected information
#' @examples
#' b1 = matrix(rnorm(500),10,50) # 10 rows and 50 columns
#' b2 = matrix(rnorm(800),10,80) # 10 rows and 80 columns
#' b3 = matrix(rnorm(700),10,70) # 10 rows and 70 columns
#' b4 = matrix(rnorm(500),10,50) # 10 rows and 50 columns
#' # Build multi-block by adding one data block at a time:
#' mb <- BuildMultiBlock(b1, newSamples = paste0('sample_',1:10))
#' mb <- BuildMultiBlock(b2, growingMB = mb)
#' mb2 <- BuildMultiBlock(b3, newSamples = paste0('sample_',1:10))
#' mb2 <- BuildMultiBlock(b4, growingMB = mb2)
#' # Combine the two multi-blocks:
#' mb3 <- CombineMultiBlock(newMB = mb1, growingMB = mb2)
#' @export

  if(class(model) != 'ComDim'){
    stop("'model' is  not the output of a ComDim analysis.")
  }

  if(is.null(dim)){
    dim <- 1:model@ndim
  }

  if(is.null(dim.ort) && length(model@Orthogonal) != 0){
    if(model@Orthogonal$nort != 0){
      dim.ort <- 1:model@Orthogonal$nort
    }
  } else if(is.null(dim.ort) && length(model@Orthogonal) == 0){
    dim.ort <- 0
  }

  block_names <- names(model@T.scores)

  if(!(is.null(blocks))){
    if(is.character(blocks)){
      block.ids <- match(blocks, block_names)
      if(any(is.na(block.ids))){
        warning('Some block names were not found in the MultiBlock: %s',
                paste0(blocks[which(is.na(block.ids))], sep=', '))
        block.ids <- block.ids[!is.na(block.ids)]
      }
    } else if(is.numeric(blocks)){
      block.ids <- intersect(1:length(block_names), blocks)
      if(length(blocks) != length(block.ids)){
        warning('Some block names were not found in the MultiBlock: %s',
                paste0(setdiff(1:length(block_names), blocks), sep=', '))
      }
    }
  } else {
    block.ids <- seq(1,length(block_names))
  }


  if(length(block.ids) == 0){
    stop('0 blocks were selected. The table cannot be built.')
  }


  for(i in block.ids){ # For each of the dimensions

    posvars <- which(model@variable.block == block_names[i])

    for(v in dim){
      df2 <- data.frame(
        variable.id = rownames(model@P.loadings[posvars,]),
        variable.id.number = posvars,
        block.id = rep(i, length(posvars)),
        block.name = rep(block_names[i], length(posvars)),
        dim = rep(v, length(posvars)),
        value = model@P.loadings[posvars,v]
      )
      if(v == dim[1] && i == block.ids[1]){
        df <- df2
      } else {
        df <- rbind.data.frame(df,df2)
      }
    }

    if(dim.ort != 0){
      for(v in dim.ort){
        df2 <- data.frame(
          variable.id = rownames(model@Orthogonal$P.loadings.ort[posvars,]),
          variable.id.number = posvars,
          block.id = rep(i, length(posvars)),
          block.name = rep(block_names[i], length(posvars)),
          dim = rep(v, length(posvars)),
          value = model@Orthogonal$P.loadings.ort[posvars,v]
        )
        if(v == dim[1] && i == block.ids[1]){
          df <- df2
        } else {
          df <- rbind.data.frame(df,df2)
        }
      }
    }
  }
  rownames(df) <- NULL


  # Add factors to data
  df$variable.id <- factor(df$variable.id, levels = unique(df$variable.id))
  df$variable.id.number <- factor(df$variable.id.number)
  df$block.name <- factor(df$block.name, levels = block_names)

  return(df)

}
