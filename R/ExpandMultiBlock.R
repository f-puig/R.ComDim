ExpandMultiBlock <- function(data = NULL, metadata = NULL, minblock = 0, loquace = TRUE){

  #' ExpandMultiBlock
  #' Splits data into several blocks, allowing variables to appear on more than one block at the same time.
  #' @param data A data.frame or a matrix. Samples are in rows and variables in columns.
  #' @param metadata A 2-column data.frame detailing how columns in data must be propagated. Entries in the first column must match the row names in data (variable  names). The second column details the new blocks.
  #' @param minblock The minimal number of variables that need to be contained in a block in order to be kept.
  #' @param loquace To display warning
  #' @return The MultiBlock.
  #' @examples
  #' data(mouse_ds)
  #' lipidsMB <- ExpandMB(data = lipids, metadata = metadata_lipids,
  #'   minblock = 0, loquace = FALSE)
  #' @export

  if(is.null(data)){
    stop('No data was provided.')
  } else {
    if(any(duplicated(rownames(data)))){
      stop('There are duplicated rownames. Data must have variables in rows and samples in columns.')
    }
  }

  if(is.null(metadata)){
    stop('No metadata was provided.')
  }

  if(any(!(rownames(data) %in% metadata[,2]))){

    aa <- rownames(data)[!(rownames(data) %in% metadata[,2])]
    if(loquace){
      warning(sprintf('%s row names are not present in metadata: %s',
                      length(aa),
                      paste0(aa, collapse = ", ")))
    }

    data <- data[rownames(data) %in% metadata[,2],]
  }

  if(any(!(metadata[,2] %in% rownames(data)))){
    metadata <- metadata[ metadata[,2] %in% rownames(data),]
  }

  if(minblock != 0){
    xx <- as.data.frame(table(metadata[,1]))
    xx <- xx[xx$Freq >= minblock,]
    metadata <- metadata[metadata[,1] %in% xx$Var1,]
    data <- data[rownames(data) %in% metadata[,2],]
  }

  expanded_data <- matrix(ncol = nrow(metadata), nrow = ncol(data))

  for(i in 1:nrow(data)){
    pos <- which(metadata[,2] == rownames(data)[i])
    expanded_data[,pos] <- kronecker(matrix(1,1,length(pos)),as.matrix(data[i,]))
    # Similar function to repmat from base R
  }

  colnames(expanded_data) <- paste(metadata[,1], metadata[,2], sep = ".")
  rownames(expanded_data) <- colnames(data)

  # Remove rows with all NAs
  pos_na <- which(apply(expanded_data, 2, function(x) all(is.na(x))))

  if(length(pos_na) != 0){
    data <- expanded_data[,-pos_na]
    metadata <- metadata[-pos_na,]
  } else {
    data <- expanded_data
  }

  splitMB <- SplitVariablesMB(data = data,
                              metadata = metadata[,1],
                              minblock = minblock)

  return(splitMB)

}
