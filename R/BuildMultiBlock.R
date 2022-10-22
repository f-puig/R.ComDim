BuildMultiBlock <- function(..., ignore.names = FALSE, ignore.size = FALSE) {


#' BuildMultiBlock
#' Creates a MultiBlock object
#' @param ... Data.frames, matrices, or MultiBlocks to be combined.
#' @param ignore.names The MultiBlock is created without checking the sample names. Names are replaced by integers. Useful when the nomenclature to name samples is not consistent across blocks but there exist sample correspondence across blocks.
#' @param ignore.size The MultiBlock is created without checking the sample size. This is set as TRUE when dealing with replicate data and we want to consider each replicate batch as a stand alone block. See SplitRW().
#' @return The MultiBlock.
#' @examples
#'  b1 = matrix(rnorm(500),10,50) # 10 rows and 50 columns
#'  b2 = as.data.frame(matrix(rnorm(800),10,80)) # 10 rows and 80 columns
#'  b3 = new("MultiBlock", Samples = 1:10, Data = list(s1 = b1, s2 = b1),
#'  Variables = list(s1 = 1:ncol(b1), s2 = letters[1:length(b1)]))
#'  # Build MultiBlock by adding one data block at a time:
#'  mb <- BuildMultiBlock(b1, b2, b3)
#' @export

  if(length(list(...)) != 0){
    add_param <- list(...) # The different blocks to combine

    names_param <- sapply(substitute(list(...))[-1], deparse)

    if(any(duplicated(names_param))){
      stop('The different blocks to combine must have different names.')
    }
  } else {
    stop("No blocks were provided.")
  }

  # First step. Check if the provided elements are data.frames, matrices, Blocks, or Multiblocks.
  # If any is none of the above, then stop.
  for(i in 1:length(add_param)){
    if(!(any(class(add_param[[i]]) %in% c("data.frame", "matrix", "MultiBlock")))){
      stop("At least one of the provided blocks has not a suitable class.
           It must be 'data.frame', 'matrix', or 'MultiBlock'")
    }
  }

  # Now prepare one block at a time.
  # Capture sample names and block length
  samples_list <- list()
  variables_list <- list()
  block_length <- vector()
  #samples_numeric <- TRUE
  k <- 1 # Number of blocks

  for(i in 1:length(add_param)) {

    if(any(class(add_param[[i]]) == "matrix" || class(add_param[[i]]) == "data.frame")) {
      if(is.data.frame(add_param[[i]])){
        if(any('character' %in% lapply(add_param[[i]], class))){
          stop('There is at least one data.frame with character data.')
        }
        add_param[[i]] <- as.matrix(add_param[[i]])
      }

      if(!is.null(rownames(add_param[[i]]))){
        samples_list[[k]] <- rownames(add_param[[i]])
      } else {
        samples_list[[k]] <- 1:nrow(add_param[[i]])
      }

      if(!is.null(colnames(add_param[[i]]))){
        newVars <- colnames(add_param[[i]])
      } else {
        newVars <- 1:ncol(add_param[[i]])
      }
      if(length(dimnames(add_param[[i]])) != 0){
        dimnames(add_param[[i]]) <- NULL
      }


      if(i == 1){
        rownames(add_param[[i]]) <- NULL
        colnames(add_param[[i]]) <- NULL

        Data <- list(x = add_param[[i]])

        growingMB <- MultiBlock(Samples = samples_list[[k]],
                                Data = Data,
                                Variables = list(x = newVars))
        names(growingMB@Data)[k] <- names_param[i]
        names(growingMB@Variables)[k] <- names_param[i]

      } else {
        rownames(add_param[[i]]) <- NULL
        colnames(add_param[[i]]) <- NULL
        growingMB@Data[[ names_param[i] ]] <- add_param[[i]]
        if(ncol(add_param[[i]]) == 1){
          growingMB@Data[[ names_param[i] ]] <- as.matrix(growingMB@Data[[ names_param[i] ]])
        }
        growingMB@Variables[[ names_param[i] ]] <- newVars
      }

      block_length <- append(block_length, nrow(add_param[[i]]))
      k <- k + 1

    } else { # If it is a multiblock

      if(k > 1) {
        if(any(names(add_param[[i]]@Data) %in% names(growingMB@Data))){
          dupe <- which(names(add_param[[i]]@Data) %in% names(growingMB@Data))
          names(add_param[[i]]@Data)[dupe] <- paste0(names(add_param[[i]]@Data)[dupe], ".new")
          warning("There were at least 2 blocks with the same name. The suffix '.new' was added to one of them.")
        }
      }

      for(j in 1:length(add_param[[i]]@Data)){

        samples_list[[k]] <- add_param[[i]]@Samples

        xx <- names(add_param[[i]]@Data)[j]

        if(is.data.frame(add_param[[i]]@Data[[j]])){
          add_param[[i]]@Data[[j]] <- as.matrix(add_param[[i]]@Data[[j]])
        }

        rownames(add_param[[i]]@Data[[j]]) <- NULL
        colnames(add_param[[i]]@Data[[j]]) <- NULL

        if(j == 1 && i == 1){

          Data <- list(x = add_param[[i]]@Data[[j]])

          growingMB <- MultiBlock(Samples = 1:nrow(add_param[[i]]@Data[[j]]),
                                  Data = Data,
                                  Variables = list(x = add_param[[i]]@Variables[[j]]))
          names(growingMB@Data)[k] <- xx
          names(growingMB@Variables)[k] <- xx

        } else {
          growingMB@Data[[ xx ]] <- as.matrix(add_param[[i]]@Data[[j]])
          growingMB@Variables[[ xx ]] <- add_param[[i]]@Variables[[j]]
        }

        if(xx %in% names(add_param[[i]]@Batch)){
          growingMB@Batch[[ xx ]] <- add_param[[i]]@Batch[[ xx ]]
        }
        if(xx %in% names(add_param[[i]]@Metadata)){
          growingMB@Metadata[[ xx ]] <- add_param[[i]]@Metadata[[ xx ]]
        }

        block_length <- append(block_length, nrow(add_param[[i]]@Data[[j]]))
        k <- k + 1

      }
    }

    #if(is.character(samples_list[[k-1]])){
    #  samples_numeric <- FALSE
    #}
  }

  #if(ignore.names){
  #  samples_numeric <- TRUE
  #}

  if(ignore.names){
    if(length(unique(block_length)) != 1 && !ignore.size){
      stop("Blocks do not have the same length. If common samples across blocks
           must be selected, then provide character names for all samples
           across the blocks.")
    }

    if(length(unique(block_length)) != 1){
      growingMB@Samples <- list()
      for(i in 1:length(block_length)){
        xx <- names_param[i]
        growingMB@Samples[[ xx ]] <- samples_list[[i]]
      }
      warning("
Current object is not a suitable 'MultiBlock' object for ComDim since sample
numbers are not common across blocks. If data contains replicates, replicates
can be split into separate blocks with SplitRW(). Otherwise, provide character
names for the sample names and use BuildMultiBlock() with ignore.names = FALSE
to remove uncommon samples across blocks.")
    } else {
      growingMB@Samples <- 1:unique(block_length)
    }

    for(i in 1:length(samples_list)){
      samples_list[[i]] <- 1:nrow(growingMB@Data[[i]])
    }
  }

  if(!ignore.names){
    ignore.size <- FALSE # ignore.names=F overrides ignore.size
  }

  if(!ignore.size) {
    # Check order for samples
    common_samples <- vector()
    for(i in 1:length(samples_list)){

      if(any(duplicated(samples_list[[i]]) && ignore.names)){
        stop('There are duplicate sample names in at least one block. Set ignore.names = TRUE and/or ignore.size = TRUE.')
      }

      if(i == 1){
        common_samples <- samples_list[[i]]
      } else {
        common_samples <- intersect(common_samples, samples_list[[i]])
      }
    }
    if(length(common_samples) == 0) {
      stop('There are no common sample names across blocks.')
    }

    if(length(common_samples) != max(block_length)){
      warning('Some samples are not common across blocks. They will be discarded.')
    }

    growingMB@Samples <- common_samples
    for(i in 1:length(growingMB@Data)){
      growingMB@Data[[i]] <- as.matrix(growingMB@Data[[i]][match(common_samples,
                                                       samples_list[[i]]),])
      if(names(growingMB@Data)[i] %in% names(growingMB@Batch)){
        growingMB@Batch[[ names(growingMB@Data)[i] ]] <-
          growingMB@Batch[[ names(growingMB@Data)[i] ]][match(common_samples,
                                                              samples_list[[i]])]
      }

      if(names(growingMB@Data)[i] %in% names(growingMB@Metadata)){
        growingMB@Metadata[[ names(growingMB@Data)[i] ]] <-
          growingMB@Metadata[[ names(growingMB@Data)[i] ]][match(common_samples,
                                                                 samples_list[[i]]),]
      }
    }
  }

  if(is.list(growingMB@Samples) && !identical(names(growingMB@Samples), names(growingMB@Data))){
    mis1 <- setdiff(names(growingMB@Samples), names(growingMB@Data))
    mis2 <- setdiff(names(growingMB@Data), names(growingMB@Samples))
    stop(sprintf("Multi-Block '%s' was not properly built. '%s' fields should contain the '%s' name instead of '%s.'",
                 mis1, mis1, mis1, mis2))
  }

  return(growingMB)
}
