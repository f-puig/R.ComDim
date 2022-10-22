SummarizedExperiment2MultiBlock <-function(se,
                                           colData_samplenames = NULL, # Character string with the column metadata.
                                           Batch = NULL){ # Batch is the name of the column that will be considered as Batch
  #' SummarizedExperiment2MultiBlock 
  #'
  #' Converts a SummarizedExperiment into a MultiBlock.
  #' @param se The SummarizedExperiment .
  #' @param colData_samplenames Character string with the column metadata.
  #' @param Batch The name of the column that will be considered as Batch.
  #' @return The MultiBlock

  #' @export


  # Batch is the name inlistData with the replicate information.
  numassays <- length(assays(se))

  for(i in 1:numassays){
    Data <- t(as.matrix(assays(se)[[i]]))
    Variables <- list(colnames(Data))
    Data <- list(Data)
    block_name <- names(assays(se))[[i]]
    if(length(block_name) > 1) {
      block_name <- "X"
    }
    names(Variables) <- names(Data) <- block_name

    Data[[block_name]] <- Data[[block_name]][match(rownames(Data[[block_name]]), se@colData@rownames),]
    if(any(is.na(match(rownames(Data[[block_name]]), se@colData@rownames)))){
      Data[[block_name]] <- Data[[block_name]][!is.na(match(rownames(Data[[block_name]]), se@colData@rownames)),]
    }
    rownames(Data[[block_name]]) <- NULL
    colnames(Data[[block_name]]) <- NULL


    if(i == 1){ # First MB

      if(!is.null(Batch)){
        if(length(Batch) > 1) {
          stop("'Batch' must be a character string of length 1.")
        }
        Batch <- se@colData@listData[[ batch ]]
        noBatch <- setdiff(names(se@colData@listData), Batch)
        if(is.null(colData_samplenames)){
          MB <- MultiBlock(Samples = as.character(Metadata[[block_name]][[ colData_samplenames ]]),
                           Data = Data,
                           Variables = Variables,
                           Batch = Batch)
        } else {
          Metadata <- as.data.frame(se@colData@listData)
          Metadata <- Metadata[,noBatch]
          if(!(colData_samplenames %in% colnames(Metadata))){
            stop("The string in 'colData_samplenames' is not a column name of colData.")
          }
          Metadata <- Metadata[match(Metadata[[ colData_samplenames ]],
                                     se@colData@rownames),]
          if(any(is.na(match(Metadata[[ colData_samplenames ]],
                             se@colData@rownames)))){
            Metadata <- Metadata[!is.na(match(Metadata[[ colData_samplenames ]],
                                              se@colData@rownames)),]
          }
          Metadata <- list(Metadata)
          names(Metadata) <- block_name
          MB <- MultiBlock(Samples = as.character(Metadata[[block_name]][[ colData_samplenames ]]),
                           Data = Data,
                           Variables = Variables,
                           Batch = Batch,
                           Metadata = Metadata)
        }

      } else { # If Batch is NULL
        if(!is.null(colData_samplenames)){
          Metadata <- as.data.frame(se@colData@listData)
          Metadata <- Metadata[match(Metadata[[ colData_samplenames ]],
                                     se@colData@rownames),]
          if(any(is.na(match(Metadata[[ colData_samplenames ]],
                             se@colData@rownames)))){
            Metadata <- Metadata[!is.na(match(Metadata[[ colData_samplenames ]],
                                              se@colData@rownames)),]
          }
          Metadata <- list(Metadata)
          names(Metadata) <- block_name
          MB <- MultiBlock(Samples = as.character(Metadata[[block_name]][[ colData_samplenames ]]),
                           Data = Data,
                           Variables = Variables,
                           Metadata = Metadata)
        } else {
          MB <- MultiBlock(Samples = as.character(Metadata[[block_name]][[ colData_samplenames ]]),
                         Data = Data,
                         Variables = Variables)
        }
      }

    } else { # Iterations i > 1
      MB2 <- MultiBlock(Samples = rownames(Data[[block_name]]),
                       Data = Data,
                       Variables = Variables)
      MB <- BuildMultiBlock(MB, MB2)
    }
  }

  return(MB)

}
