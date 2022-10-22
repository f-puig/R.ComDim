setSampleNames <- function(MB, sample.names) {

  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Get the sample names.
  #' @param MB A MultiBlock object.
  #' @param sample.names A vector with the new sample names.
  #' @return The MultiBlock object.
  #' @export
  #'

  if(any(class(MB) != "MultiBlock")){
    stop("'MB' is not a MultiBlock.")
  }

  if(!is.vector(sample.names)) {
    warning('The provided sample names are not in a vector format. Sample names were not changed.')
    return(MB)
  }

  if(length(MB@Samples) != length(sample.names)) {
    warning('The length of the provided vector is not correct. Sample names were not changed.')
    return(MB)
  }

  MB@Samples <- sample.names

  return(MB)
}

