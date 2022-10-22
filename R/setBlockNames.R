setBlockNames <- function(MB, block.names) {

  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Set the block names.
  #' @param MB A MultiBlock object.
  #' @param block.names A vector with the block names.
  #' @return The MultiBlock object.
  #' @export
  #'

  if(any(class(MB) != "MultiBlock")){
    stop("'MB' is not a MultiBlock.")
  }

  if(length(names(MB@Batch)) != 0) {
    names(MB@Batch) <- block.names[
      match(names(MB@Batch), names(MB@Data))]
  }
  if(length(names(MB@Metadata)) != 0) {
    names(MB@Metadata) <- block.names[
      match(names(MB@Metadata), names(MB@Data))]
  }
  names(MB@Data) <- block.names
  names(MB@Variables) <- block.names

  return(MB)
}
