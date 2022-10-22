MultiBlock <- function(Samples = NULL,
                       Data,
                       Variables = NULL,
                       Batch = NULL,
                       Metadata = NULL){

  #' MultiBlock
  #'
  #' Combines the blocks of a MultiBlock into a matrix. Useful for calculating summary information of the whole MB (ex. maximum value).
  #' @param Samples The vector with the sample names. Facultative.
  #' @param Data A list with the data-blocks. They can be matrices, data.frames or MultiBlocks.
  #' @param Variables The vector with the variable names. Facultative.
  #' @param Batch A list with the vectors of the batch information. Facultative.
  #' @param Metadata A list with the metadata information of each data-block.
  #' @return The MultiBlock
  #' @examples
  #' b1 = matrix(rnorm(500),10,50) # 10 rows and 50 columns
  #' b2 = as.data.frame(matrix(rnorm(800),10,80)) # 10 rows and 80 columns
  #' b2[c(2,3,4),c(5,7,8)] <- NA # Making some data missing, just for fun.
  #' rownames(b2) <- LETTERS[5:14] # Adding some sample names
  #' b3 = MultiBlock(Samples = 1:10,
  #'                 Data = list(s1 = b1, s2 = b2),
  #'                 Variables = list(s1 = 1:ncol(b1),
  #'                                  s2 = 1:ncol(b2)))

  #' @export

  if(is.null(Batch)){
    Batch <- list()
  }

  if(is.null(Metadata)){
    Metadata <- list()
  }

  samples2 <- list()
  variables2 <-list()
  if(is.list(Data)){
    for(i in 1:length(Data)){
      if(is.data.frame(Data[[i]])){
        if(is.null(rownames(Data[[i]]))){
          samples2[[i]] <- 1:nrow(Data[[i]])
        } else {
          samples2[[i]] <- rownames(Data[[i]])
        }
        if(is.null(colnames(Data[[i]]))){
          variables2[[i]] <- 1:col(Data[[i]])
        } else {
          variables2[[i]] <- colnames(Data[[i]])
        }
        Data[[i]] <- as.matrix(Data[[i]])
        if(length(dimnames(Data[[i]])) != 0){
          dimnames(Data[[i]]) <- NULL
        }
      }
    }
  }

  if(is.null(Samples)){
    samples_common <- vector()
    for(i in 1:length(samples2)){
      if(i == 1){
        samples_common <- samples2[[i]]
      } else {
        samples_common <- intersect(samples_common, samples2[[i]])
      }
    }
    if(length(samples_common) == length(samples2[[1]])){
      Samples <- samples_common
    }
  }

  if(is.null(Variables)){
    Variables <- variables2
  }

  MB <- new("MultiBlock",
            Samples = Samples,
            Data = Data,
            Variables = Variables,
            Batch = Batch,
            Metadata = Metadata)

  if(validObject(MB)){
    return(MB)
  } else {
    stop("The 'MultiBlock' could not be built.")
  }
}
