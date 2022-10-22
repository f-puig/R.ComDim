MultiAssayExperiment2MultiBlock <-function(se,
                                           colData_samplenames = NULL, # Character string with the column metadata.
                                           Batch = NULL){ # Batch is the name of the column that will be considered as Batch

  #' MultiAssayExperiment2MultiBlock 
  #'
  #' Converts a MultiAssayExperiment2MultiBlock into a MultiBlock.
  #' @param se The SummarizedExperiment .
  #' @param colData_samplenames Character string with the column metadata. If colData_samplenames is NULL, the metadata will not be stored.
  #' @param Batch The name of the column that will be considered as Batch.
  #' @return The MultiBlock.

  #' @export

  # Here the metadata of the MultiAssayExperiment2MultiBlock is stored in the first MultiBlock.

  se <- intersectColumns(se[, , names(experiments(se))])

  simplelist <- assays(se)

  samplemap <- sampleMap(se)

  for(i in 1:length(simplelist@listData)){
    Data <- list(t(as.matrix(simplelist@listData[[i]])))
    names(Data) <- names(simplelist@listData)[i]
    Variables <- list(rownames(simplelist@listData[[i]]))
    names(Variables) <- names(Data)
    eq.names <- samplemap@listData$primary[match(rownames(Data[[1]]), samplemap@listData$colname)]
    Data[[1]] <- Data[[1]][match(se@colData@rownames, eq.names),]
    rownames(Data[[1]]) <- NULL
    colnames(Data[[1]]) <- NULL

    if (i == 1){
      Samples <- se@colData@rownames
      if(is.null(colData_samplenames)){
        Metadata <- NULL
      } else {
        Metadata <- as.data.frame(colData(se))
      }
      if(is.null(Metadata) && is.null(Batch)){
        MB <- MultiBlock(Samples = Samples,
                         Data = Data,
                         Variables = Variables)
      } else if(is.null(Metadata) && !is.null(Batch)){
        stop('There is no metadata in the MultiAssayExperiment to extract the Batch')
      } else if(!is.null(Metadata) && is.null(Batch)){
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
        Metadata <- list(x = Metadata)
        names(Metadata) <- names(Data)

        MB <- MultiBlock(Samples = Samples,
                         Data = Data,
                         Variables = Variables,
                         Metadata = Metadata)
      } else {
        if(batch %in% colnames(Metadata)){
          Batch <- Metadata[[ batch ]]
          noBatch <- setdiff(colnames(Metadata), Batch)
          Metadata <- Metadata[,noBatch]
          Metadata <- Metadata[match(Metadata[[ colData_samplenames ]],
                                     se@colData@rownames),]
          if(any(is.na(match(Metadata[[ colData_samplenames ]],
                             se@colData@rownames)))){
            Metadata <- Metadata[!is.na(match(Metadata[[ colData_samplenames ]],
                                              se@colData@rownames)),]
          }
          Metadata <- list(x = Metadata)
          names(Metadata) <- names(Data)
          MB <- MultiBlock(Samples = Samples,
                           Data = Data,
                           Variables = Variables,
                           Batch = Batch,
                           Metadata = Metadata)
        } else {
          stop("The string in batch was not found in the metadata.")
        }
      }
    } else {
      MB2 <- MultiBlock(Samples = Samples,
                       Data = Data,
                       Variables = Variables)
      MB <- BuildMultiBlock(MB, MB2)
    }
  }

  return(MB)
}
