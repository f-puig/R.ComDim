getBlockNames <- function(MB, DataOrBatchOrMetadata) {

  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Return the block names.
  #' @param MB A MultiBlock object.
  #' @param DataOrBatchOrMetadata A string with either "Data", "Batch" or
  #' "Metadata". Facultative.
  #' @return A vector with the block names.
  #' @export
  #'

  if(any(class(MB) != "MultiBlock")){
    stop("'MB' is not a MultiBlock.")
  }

  if(!missing(DataOrBatchOrMetadata) && length(DataOrBatchOrMetadata) > 1) {
    DataOrBatchOrMetadata <- "Data"
  } else if(missing(DataOrBatchOrMetadata)) {
    DataOrBatchOrMetadata <- "Data"
  }

  if(tolower(DataOrBatchOrMetadata) == "data") {
    theNames <- names(MB@Data)
  }

  if(tolower(DataOrBatchOrMetadata) == "batch") {
    theNames <- names(MB@Batch)
  }

  if(tolower(DataOrBatchOrMetadata) == "metadata") {
    theNames <- names(MB@Metadata)
  }

  return(theNames)
}
