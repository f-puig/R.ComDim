FilterSamplesMultiBlock <- function(MB = MB, samples =  getSampleNames(MB)) {


  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Predict_MB applies an existent ComDim model to a MultiBlock object.
  #' @param MB A MultiBlock object for which prediction is desired.
  #' @param samples A vector containing the names of the samples to keep.
  #' @return A ComDim object.
  #' @export

  if(class(MB) != 'MultiBlock'){
    stop("'MB' is  not a MultiBlock.")
  }

  all_samples <- getSampleNames(MB)
  samples <- intersect(samples, all_samples)
  block_names <- getBlockNames(MB)
  if(length(samples) == 0){
    stop('None of the provided sample names was found in MB.')
  } else {
    pos <- match(samples, all_samples)
    ntable <- length(block_names)
    MB@Samples <- MB@Samples[pos]
    for(i in 1:ntable){
      MB@Data[[ i ]] <- MB@Data[[ i ]][pos,]
      if(block_names[i] %in% names(MB@Batch)){
        MB@Batch[[i]] <- MB@Batch[[i]][pos]
      }
      if(block_names[i] %in% names(MB@Metadata)){
        MB@Metadata[[i]] <- MB@Metadata[[i]][pos,]
      }
    }
  }

  return(MB)
}
