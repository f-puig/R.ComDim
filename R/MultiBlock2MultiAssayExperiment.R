MultiBlock2MultiAssayExperiment <- function(MB,
                                            MSEmetadata = NULL){ # List of unstructured metadata describing the overall content of the object.

  #' MultiBlock2MultiAssayExperiment
  #'
  #' Converts MultiBlock into a MultiAssayExperiment.
  #' @param MB The MultiBlock object.
  #' @param MSEmetadata List of unstructured metadata describing the overall content of the object.
  #' @return The MultiAssayExperiment object.
  #'
  #' @export

  if(class(MB) != 'MultiBlock'){
    stop("'MB' is  not a MultiBlock.")
  }

  ## Create ExperimentList
  exp <- list()
  for(i in 1:length(MB@Data)){
    exp[[ names(MB@Data)[i] ]] <- t(MB@Data[[i]])
    colnames( exp[[ names(MB@Data)[i] ]] ) <- MB@Samples
    rownames( exp[[ names(MB@Data)[i] ]] ) <- MB@Variables[[ names(MB@Data)[i] ]]
  }
  exp <- ExperimentList(exp)

  ## Create colData
  colDat <- matrix(MB@Samples, nrow = length(MB@Samples), ncol = 1)
  for(i in 1:length(MB@Data)){
    if(names(MB@Data)[i] %in% names(MB@Metadata)){
      colDat <- cbind.data.frame(colDat, MB@Metadata[[ names(MB@Metadata)[i] ]])
    }
    if(names(MB@Data)[i] %in% names(MB@Batch)){
      colDat <- cbind.data.frame(colDat, MB@Batch[[ names(MB@Metadata)[i] ]])
    }
  }
  rownames(colDat) <- MB@Samples

  # SampleMap
  for(i in 1:length(MB@Data)){
    if(i == 1){
      SEsampleMap <- data.frame(assay = rep(names(MB@Data)[i], length(MB@Samples)),
                                primary = MB@Samples,
                                colname = MB@Samples)
    } else {
      SEsampleMap <- rbind.data.frame(SEsampleMap,
                                      data.frame(assay = rep(names(MB@Data)[i],
                                                             length(MB@Samples)),
                                                 primary = MB@Samples,
                                                 colname = MB@Samples))
    }
  }


  if(is.null(MSEmetadata)){
    se <- MultiAssayExperiment(experiments=exp,
                               colData=colDat,
                               sampleMap=SEsampleMap)
  } else {
    se <- MultiAssayExperiment(experiments=exp,
                               colData=colDat,
                               sampleMap=SEsampleMap,
                               metadata=MSEmetadata)
  }

  return(se)
}
