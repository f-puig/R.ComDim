getSampleNames <- function(MB) {

  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Get the sample names.
  #' @param MB A MultiBlock object.
  #' @return A vector with the Sample names.
  #' @export
  #'

  if(any(class(MB) != "MultiBlock")){
    stop("'MB' is not a MultiBlock.")
  }

  theNames <- MB@Samples
  return(theNames)
}

