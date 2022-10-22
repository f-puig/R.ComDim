SplitRW <- function( MB = MB, checkSampleCorrespondence = FALSE,
                        batchNormalisation = TRUE, showSampleCorrespondence = TRUE) {

#' SplitRW - Building the RW list compatible with ComDim_PCA().
#'
#' To split a multi-block into a list containing smaller blocks, comprising each data from one batch.
#' @param  MB The multi-block structure built with BuildMultiblock().
#' @param  checkSampleCorrespondence If FALSE, the same number of samples and the sample order are assumed for all the batches. If TRUE, only the samples found in all replicate blocks will be included in the final structure.
#' @param  batchNormalisation If TRUE, each RW block is divided by its norm and the square root of the number of replicates.
#' @param  showSampleCorrespondence If TRUE, print in the console a data frame with the sample correspondence across blocks.
#' @return The list containing the multi-block composed of the replicate-blocks.
#' @examples
#' b1 = matrix(rnorm(500),10,50)
#' batch_b1 = rep(1,10)
#' b2 = matrix(rnorm(800),30,80)
#' batch_b2 = c(rep(1,10),rep(2,10),rep(3,10))
#' # Generate the multi-block (mb)
#' mb <- BuildMultiBlock(b1, batches = batch_b1)
#' mb <- BuildMultiBlock(b2, growingMB = mb, batches = batch_b2, equalSampleNumber = FALSE)
#' rw <- SplitRW(mb)
#' @export

  newMB <- list()

  if(class(MB) != 'MultiBlock'){
    stop("'MB' is  not a MultiBlock.")
  }

  if(!is.list(MB@Samples)){
    samples <- MB@Samples
    MB@Samples <- list()
    xx <- names(MB@Data)
    for(i in 1:length(xx)){
      MB@Samples[[ xx[i] ]] <- samples
    }
    rm(samples)
  }

  # Check that the batch information is correctly provided for all the blocks:
  give_error <- 0

  batch_names <- list()
  samples_per_batch <- list()

  if(length(MB@Batch) != 0){
    for (i in names(MB@Data)){ # Check the batch information is consistent across RWs
      if(i %in% names(MB@Batch)){
        batch_names[[i]] <- sort(unique(MB@Batch[[i]]))
        sampled<-as.vector(NULL)
        for (j in 1:length(batch_names[[i]])){
          sampled[j] <- length(which(MB@Batch[[i]] == batch_names[[i]][j]))
        }
        samples_per_batch[[i]] <- sampled
        if(length(unique(sampled)) > 1){
          warning(sprintf('The replicate blocks in %s have a different number of samples.', names(MB@Data)[i]))
        }
      } else {
        batch_names[[i]] <- "Batch 1"
        samples_per_batch[[i]] <- MB@Samples[[i]]
      }
    }
    if(length(unique(unlist(samples_per_batch))) > 1 && checkSampleCorrespondence == FALSE){
      print('Information is missing regarding the splitting.')
      print('Using checkSampleCorrespondence as TRUE is recommended.')
      stop('The data cannot be split into replicate blocks.')
    }

    for (i in names(MB@Batch)){
      if(checkSampleCorrespondence == TRUE){
        for(j in 1:length(batch_names[[i]])){
          batch_position <- which(MB@Batch[[i]] == batch_names[[i]][j])
          if (i == names(MB@Batch)[1] && j == 1) {
            replicate_names <- MB@Samples[[i]][batch_position]
          } else {
            replicate_names <- intersect(replicate_names, MB@Samples[[i]][batch_position])
          }
        }
      } else {
        replicate_names <- 1:min(lengths(MB@Samples)) # Used just to count the number of samples in the smallest block.
      }
    }

    if (checkSampleCorrespondence == TRUE && length(replicate_names)== 0) {
      print('There are 0 samples in common across the replicate blocks')
      give_error <- 1
    }

    if (give_error) {
      stop('The data cannot be split into replicate blocks.')
    }
  } else {
    warning("The MultiBlock does not contain 'Batch' information. It cannot be split into replicate blocks.")
    return(MB)
  }


  ## PROCEED WITH THE RW SPLITTING.
  # df_SampleNames is a table printed in the console (go to the end of the script for more info)
  if(showSampleCorrespondence){
    df_SampleNames <- matrix(, nrow = length(replicate_names), ncol = 0)
  }

  k <-1
  for (i in names(MB@Batch)){

    for (j in 1:length(batch_names[[i]])){

      batch_position <- which(MB@Batch[[i]] == batch_names[[i]][j])
      replicate_position <- as.vector(NULL)
      if (checkSampleCorrespondence == TRUE) {
        for(namei in replicate_names){
          replicate_position <- c(replicate_position,
                                  intersect(which(MB@Samples[[i]] == namei), batch_position))
        }
      } else {
        replicate_position <- batch_position
      }
      # Keep only the common samples across blocks
      if (length(replicate_position) != length(replicate_names)){
        print(sprintf('There are sample duplicates in batch %s from block %s.',as.character(batch_names[[i]][j]),as.character(i)))
        print('Duplicate samples should be removed.')
        stop('Existence of duplicate samples within one or more batches.')
      }

      sorted <- order(MB@Samples[[i]][replicate_position])

      df_SampleNames <- cbind(df_SampleNames, MB@Samples[[i]][replicate_position[sorted]])

      growingMB <- MultiBlock(Samples = replicate_names,
                              Data = list(s1 = MB@Data[[i]][replicate_position[sorted],]),
                              Variables = list(s1 = MB@Variables[[i]]),
                              Batch = list(s1 = rep(batch_names[[i]][j],length(sorted)))
                              )

      if(length(batch_names[[i]]) == 1){
        growingMB <- setBlockNames(growingMB, i)
      } else {
        growingMB <- setBlockNames(growingMB,
                                   paste(i, as.character(batch_names[[i]][j]), sep='_')
                                   )
      }

      if(batchNormalisation){ # Divide each RW block by its norm and the sqrt(number_of_replicates) to avoid overrepresentation of the data in the multi-set.
        normed <- growingMB@Data[[1]] / norm(growingMB@Data[[1]], type = "F")
        growingMB@Data[[1]] <- normed / sqrt(length(batch_names[[i]]))
      }

      #if(length(names(MB[[i]])) > 5){ # In case there are additional fields
      if(length(MB@Metadata[[i]]) != 0){ # In case there is metadata

        growingMB@Metadata[[getBlockNames(growingMB)]] <-  MB@Metadata[[i]][replicate_position[sorted],]

      }
      k <- k + 1

      if(j == 1){
        mbb <- growingMB
      } else {
        mbb <- BuildMultiBlock(mbb, growingMB)
      }
    }
    if (i == names(MB@Batch)[1]){
      newMB <- mbb
    } else {
      newMB <- BuildMultiBlock(newMB, mbb)
    }
  }


  ## SHOW A DATA-FRAME WITH ALL THE SAMPLE NAMES FOR ALL THE (REPLICATE) BLOCKS.
  ## It is wise to use it when checkSampleCorrespondence is set to FALSE.
  if(showSampleCorrespondence){
    colnames(df_SampleNames) <- names(newMB@Data)
    print('The sample names are:')
    print(df_SampleNames)
  }

  return(newMB)
}
