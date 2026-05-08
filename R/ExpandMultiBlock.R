#' ExpandMultiBlock
#'
#' Splits data into several blocks, allowing variables to appear in more than one block
#' simultaneously. Each variable is duplicated into every block to which it is assigned according
#' to the \code{metadata} mapping table.
#' @param data A data.frame or matrix with samples in rows and variables in columns.
#' @param metadata A 2-column data.frame describing how variables are assigned to blocks. The
#'   first column gives the block name; the second column gives the variable name, and must match
#'   the column names of \code{data}. A variable may appear in multiple rows (and therefore in
#'   multiple blocks).
#' @param minblock Integer. Blocks with fewer than \code{minblock} variables are discarded. Use
#'   \code{0} (default) to keep all blocks regardless of size.
#' @param loquace Logical. If \code{TRUE} (default), a warning is printed listing column names in
#'   \code{data} that have no matching entry in the second column of \code{metadata}.
#' @details
#'   For each row in \code{metadata} that matches a variable in \code{data}, the variable's
#'   values are copied into the corresponding block column. Column names in the resulting expanded
#'   matrix are formed as \code{<block>.<variable>}. Variables with all-\code{NA} values after
#'   expansion are removed. If no matches exist between \code{data} column names and
#'   \code{metadata}, \code{NULL} is returned with a warning. If \code{minblock} filtering
#'   removes all blocks, \code{NULL} is returned.
#' @return A \code{MultiBlock} object whose blocks are defined by the first column of
#'   \code{metadata}, or \code{NULL} if no valid blocks could be constructed.
#' @seealso \code{\link{MultiBlock}}, \code{\link{ProcessMultiBlock}}
#' @examples
#' data(mouse_ds)
#' lipidsMB <- ExpandMultiBlock(data = lipids, metadata = metadata_lipids,
#'   minblock = 0, loquace = FALSE)
#' @export
ExpandMultiBlock <- function(data = NULL, metadata = NULL, minblock = 0, loquace = TRUE) {
  if (is.null(data)) {
    stop("No data was provided.")
  }

  if (is.null(metadata)) {
    stop("No metadata was provided.")
  }

  if (is.numeric(metadata[, 2])) {
    metadata[, 2] <- as.character(metadata[, 2])
  }

  meta_vars <- as.character(metadata[, 2])

  if (any(duplicated(rownames(data)))) {
    stop("There are duplicated rownames. Data must have samples in rows and variables in columns.")
  }

  if (!any(colnames(data) %in% meta_vars)) {
    warning(paste(
      "0 coincidences were found between the column names of 'data' and the",
      "second column of 'metadata'. Either the input data or the metadata",
      "data.frame may be in an incorrect format."
    ))
    return(NULL)
  }

  if (any(!(colnames(data) %in% meta_vars))) {
    aa <- colnames(data)[!(colnames(data) %in% meta_vars)]
    if (loquace) {
      warning(sprintf(
        "%d column names are not present in metadata: %s",
        length(aa),
        paste0(aa, collapse = ", ")
      ))
    }
    data <- data[, colnames(data) %in% meta_vars, drop = FALSE]
  }

  if (any(!(meta_vars %in% colnames(data)))) {
    metadata <- metadata[meta_vars %in% colnames(data), ]
    meta_vars <- as.character(metadata[, 2])
  }

  if (minblock != 0) {
    xx <- as.data.frame(table(metadata[, 1]))
    xx <- xx[xx$Freq >= minblock, ]
    metadata  <- metadata[metadata[, 1] %in% xx$Var1, ]
    meta_vars <- as.character(metadata[, 2])
    data <- data[, colnames(data) %in% meta_vars, drop = FALSE]
  }

  expanded_data <- matrix(ncol = nrow(metadata), nrow = nrow(data))

  for (i in seq_len(ncol(data))) {
    pos <- which(meta_vars == colnames(data)[i])
    expanded_data[, pos] <- kronecker(matrix(1, 1, length(pos)), as.matrix(data[, i]))
  }

  colnames(expanded_data) <- paste(metadata[, 1], meta_vars, sep = ".")
  rownames(expanded_data) <- rownames(data)

  # Remove columns with all NAs
  pos_na <- which(apply(expanded_data, 2, function(x) all(is.na(x))))
  if (length(pos_na) != 0) {
    expanded_data <- expanded_data[, -pos_na]
    metadata <- metadata[-pos_na, ]
    meta_vars <- as.character(metadata[, 2])
  }

  block_assignments <- as.character(metadata[, 1])
  unique_blocks <- unique(block_assignments)

  data_list <- list()
  for (bn in unique_blocks) {
    idx <- which(block_assignments == bn)
    mat <- expanded_data[, idx, drop = FALSE]
    colnames(mat) <- meta_vars[idx]
    data_list[[bn]] <- mat
  }

  if (length(data_list) == 0) {
    warning("0 blocks were made. Check whether the minblock threshold was set too high.")
    return(NULL)
  }

  MB <- MultiBlock(Data = data_list)
  return(MB)
}
