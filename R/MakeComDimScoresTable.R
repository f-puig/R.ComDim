MakeComDimScoresTable <- function(model,
                            blocks = NULL,
                            dim = NULL,
                            dim.ort = NULL,
                            include = c('Q.scores','T.scores', 'Q.scores.ort', 'T.scores.ort')[1:2]) {

#' MakeComDimTable
#'
#' Creates a long dataframe with the global and local scores (Q.scores and T.scores) suitable for ggplot.
#' @param model The output from a ComDim analysis.
#' @param blocks The blocks from which the data will be extracted. It can be a vector of integers or a vector with the block names. Facultative.
#' @param dim The dimensions from which the data will be extracted. It can be a vector of integers. Facultative.
#' @param dim.ort The orthogonal dimensions from which the data will be extracted. It can be a vector of integers. If NULL, all orthogonal components are extracted. If 0, none of the orthogonal components are extracted. Facultative.
#' @param include The list of ComDim information that will be added into a table. Facultative.
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

  if(is.null(include)){
    include <- c('Q.scores','T.scores')
  } else { # Added to match words regardless upper or lowercase
    include <- toupper(include)
    include <- gsub('T.SCORES','T.scores',include)
    include <- gsub('Q.SCORES','Q.scores',include)
    include <- gsub('ORT','ort',include)
  }

  if(any(!(include %in% c('Q.scores','T.scores','Q.scores.ort','T.scores.ort')))){
    stop("Some elements in 'include' were not recognized'.")
  }


  if(!(is.null(blocks))){
    if(is.character(blocks)){
      block.ids <- match(blocks, names(model@T.scores))
      if(any(is.na(block.ids))){
        warning('Some block names were not found in the MultiBlock: %s',
                paste0(blocks[which(is.na(block.ids))], sep=', '))
        block.ids <- block.ids[!is.na(block.ids)]
      }
    } else if(is.numeric(blocks)){
      block.ids <- intersect(1:length(names(model@T.scores)), blocks)
      if(length(blocks) != length(block.ids)){
        warning('Some block names were not found in the MultiBlock: %s',
                paste0(setdiff(1:length(names(model@T.scores)), blocks), sep=', '))
      }
    }
  } else {
    block.ids <- seq(1,length(names(model@T.scores)))
  }

  if(length(block.ids) == 0){
    stop('0 blocks were selected. The table cannot be built.')
  }

  # Check number of samples is always the same across all Q and T.
  nsamples <- nrow(model@Q.scores)
  for(i in block.ids){
    if(nsamples != nrow(model@T.scores[[i]])){
      stop("'model' is not an output from a ComDim analysis.")
    }
  }

  df <- NULL
  if('Q.scores' %in% include){
    for(v in dim){ # For each of the dimensions
      df2 <- data.frame(
        sample.id = rownames(model@Q.scores),
        sample.id.number = seq(1,nsamples),
        block.id = rep('Global', nsamples),
        block.name = rep('Global', nsamples),
        dim = rep(v, nsamples),
        scores.type = rep('Global', nsamples),
        scores.type.dim = paste0(rep('Q.scores', nsamples), rep(v, nsamples)),
        value = model@Q.scores[,v]
        )
      if(is.null(df)){
        df <- df2
      } else {
        df <- rbind.data.frame(df, df2)
      }
    }
  }

  if(dim.ort != 0){
    if('Q.scores.ort' %in% include){
      for(v in dim.ort){ # For each of the dimensions
        df2 <- data.frame(
          sample.id = rownames(model@Q.scores),
          sample.id.number = seq(1,nsamples),
          block.id = rep('Global.ort', nsamples),
          block.name = rep('Global.ort', nsamples),
          dim = rep(v, nsamples),
          scores.type = rep('Global.ort', nsamples),
          scores.type.dim = paste0(rep('Q.scores.ort', nsamples), rep(v, nsamples)),
          value = model@Orthogonal$Q.scores[,v]
        )
        if(is.null(df)){
          df <- df2
        } else {
          df <- rbind.data.frame(df, df2)
        }
      }
    }
  }

  if('T.scores' %in% include){
    for(i in block.ids){ # For each of the dimensions
      for(v in dim){
        df2 <- data.frame(
          sample.id = rownames(model@Q.scores),
          sample.id.number = seq(1,nsamples),
          block.id = rep(i, nsamples),
          block.name = rep(names(model@T.scores)[i], nsamples),
          dim = rep(v, nsamples),
          scores.type = rep('Local', nsamples),
          scores.type.dim = paste0(rep('T.scores', nsamples), rep(v, nsamples)),
          value = model@T.scores[[i]][,v]
          )
        if(is.null(df)){
          df <- df2
        } else {
          df <- rbind.data.frame(df, df2)
        }
      }
    }
  }

  if(dim.ort != 0){
    if('T.scores.ort' %in% include){
      for(i in block.ids){ # For each of the dimensions
        for(v in dim.ort){
          df2 <- data.frame(
            sample.id = rownames(model@Q.scores),
            sample.id.number = seq(1,nsamples),
            block.id = rep(i, nsamples),
            block.name = rep(names(model@Orthogonal$T.scores)[i], nsamples),
            dim = rep(v, nsamples),
            scores.type = rep('Local.ort', nsamples),
            scores.type.dim = paste0(rep('T.scores.ort', nsamples), rep(v, nsamples)),
            value = model@Orthogonal$T.scores[[i]][,v]
            )
          if(is.null(df)){
            df <- df2
          } else {
            df <- rbind.data.frame(df, df2)
          }
        }
      }
    }
  }
  rownames(df) <- NULL

  # Add factors to data
  df$sample.id <- factor(df$sample.id, levels = df$sample.id[seq(1,nsamples)])
  df$sample.id.number <- factor(df$sample.id.number, levels = c('Global',seq(1,nsamples)))
  df$block.name <- factor(df$block.name, levels = c(names(model@T.scores), 'Global'))
  df$scores.type <- factor(df$scores.type, levels = unique(df$scores.type))
  df$scores.type.dim <- factor(df$scores.type.dim, levels = unique(df$scores.type.dim))

  return(df)

}
