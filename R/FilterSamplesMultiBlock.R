#' FilterSamplesMultiBlock
#'
#' Retain a subset of samples in a MultiBlock object.
#' @param MB A MultiBlock object.
#' @param samples A vector of sample names to keep. Names not found in the
#'   MultiBlock are silently ignored; the order of the retained samples
#'   follows the order given here.
#' @return The MultiBlock object restricted to the requested samples. The \code{Batch} and
#'   \code{Metadata} slots are also subsetted and reordered to match.
#' @seealso \code{\link{MultiBlock}}, \code{\link{AddMetadata}}, \code{\link{ProcessMultiBlock}}
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' rownames(b1) <- rownames(b2) <- paste0("sample_", 1:10)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' mb <- FilterSamplesMultiBlock(mb, samples = paste0("sample_", 1:5))
#' @export
FilterSamplesMultiBlock <- function(MB, samples = sampleNames(MB)) {
  if (!("MultiBlock" %in% class(MB))) {
    stop("'MB' is not a MultiBlock object.")
  }

  all_samples <- sampleNames(MB)
  samples <- intersect(samples, all_samples)
  block_names <- blockNames(MB)

  if (length(samples) == 0) {
    stop("None of the provided sample names was found in MB.")
  }

  pos <- match(samples, all_samples)
  ntable <- length(block_names)

  MB@Samples <- MB@Samples[pos]

  for (i in seq_len(ntable)) {
    MB@Data[[i]] <- MB@Data[[i]][pos, , drop = FALSE]
    if (block_names[i] %in% names(MB@Batch)) {
      MB@Batch[[block_names[i]]] <- MB@Batch[[block_names[i]]][pos]
    }
    if (block_names[i] %in% names(MB@Metadata)) {
      MB@Metadata[[block_names[i]]] <- MB@Metadata[[block_names[i]]][pos, , drop = FALSE]
    }
  }

  return(MB)
}
