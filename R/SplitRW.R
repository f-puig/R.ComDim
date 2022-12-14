SplitRW <- function( MB = MB, checkSampleCorrespondence = FALSE, batchNormalisation = TRUE) {

#' SplitRW - Building the RW list compatible with ComDim_PCA().
#'
#' To split a multi-block into a list containing smaller blocks, comprising each data from one batch.
#' @param  MB The multi-block structure built with BuildMultiblock().
#' @param  checkSampleCorrespondence If FALSE, the same number of samples and the sample order are assumed for all the batches. If TRUE, only the samples found in all replicate blocks will be included in the final structure.
#' @param  batchNormalisation If TRUE, each RW block is divided by its norm and the square root of the number of replicates.
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

  # Check MB
  if(is.list(MB)){
    for(i in 1:length(MB)){
      if(is.list(MB[[i]])){
        namesMB <- names(MB[[i]])
        if (is.null(MB[[i]]$Batch)){
          print(sprintf('The Batch information was not included for Block %s.', i))
          stop('This information should be included using the batches parameter.')
        }
        if(any(!(namesMB[1:3] %in% c('Data','Variables','Samples')))){
          stop(sprintf('Information is missing in block %s.', names(MB)[i]))
        } else {
          if(nrow(MB[[i]]$Data) != length(MB[[i]]$Samples)){
            stop(sprintf('The length of the Samples vector in block %s is incorrect.',i))
          }
          if(ncol(MB[[i]]$Data) != length(MB[[i]]$Variables)){
            stop(sprintf('The length of the Variables vector in block %s is incorrect.',i))
          }
          if(nrow(MB[[i]]$Data) != length(MB[[i]]$Batch)){
            stop(sprintf('The length of the Batches vector in block %s is incorrect.',i))
          }
        }
      } else {
        print(sprintf('The block %s is not a list.', i))
        stop(sprintf('The block %s must be a list.', i))
      }
    }
  } else {
    print('The provided multi-block is not a list.')
    stop('The multi-block must be a list.')
  }

  # Check that the batch information is correctly provided for all the blocks:
  give_error <- 0

  batch_names <- list()
  samples_per_batch <- list()

  for (i in 1:length(MB)){ # Check the batch information is consistent across RWs
    batch_names[[i]] <- sort(unique(MB[[i]]$Batch))
    sampled<-as.vector(NULL)
    for (j in 1:length(batch_names[[i]])){
      sampled[j] <- length(which(MB[[i]]$Batch == batch_names[[i]][j]))
    }
    samples_per_batch[[i]] <- sampled
    if(length(unique(sampled)) > 1){
      warning(sprintf('The replicate blocks in %s have a different number of samples.', names(MB)[i]))
    }
  }
  if(length(unique(unlist(samples_per_batch))) > 1 && checkSampleCorrespondence == FALSE){
    print('Information is missing regarding the splitting.')
    print('Using checkSampleCorrespondence as TRUE is recommended.')
    stop('The data cannot be split into replicate blocks.')
  }

  for (i in 1:length(MB)){
    if(checkSampleCorrespondence == TRUE){
      for(j in 1:length(batch_names[[i]])){
        batch_position <- which(MB[[i]]$Batch == batch_names[[i]][j])
        if (i == 1 && j == 1) {
          replicate_names <- MB[[i]]$Samples[batch_position]
        } else {
          replicate_names <- intersect(replicate_names, MB[[i]]$Samples[batch_position])
        }
      }
    }
  }

  if (checkSampleCorrespondence == TRUE && length(replicate_names)== 0) {
    print('There are 0 samples common across the replicate blocks')
    give_error <- 1
  }

  if (give_error) {
    stop('The data cannot be split into replicate blocks.')
  }

  ## PROCEED WITH THE RW SPLITTING.
  k <-1
  for (i in 1:length(MB)){
    for (j in 1:length(batch_names[[i]])){

      growingMB <- list() # The Replicate block

      batch_position <- which(MB[[i]]$Batch == batch_names[[i]][j])
      if (checkSampleCorrespondence == TRUE) {
        replicate_position <- as.vector(NULL)
        for(pos in 1:length(replicate_names)){
          replicate_position <- c(replicate_position,grep(replicate_names[pos],MB[[i]]$Samples))
        } # Keep only the common samples across blocks
        batch_position <- intersect(batch_position,replicate_position)
        if (length(batch_position) != length(replicate_names)){
          print(sprintf('There are sample duplicates in batch %s from block %s.',as.character(batch_names[[i]][j]),as.character(i)))
          print('Duplicate samples should be removed.')
          stop('Existence of duplicate samples within one or more batches.')
        }
      }

      sorted <- order(MB[[i]]$Samples[batch_position])

      #growingMB$Name <- paste(MB[[i]]$Name, as.character(batch_names[[i]][j]), sep='_R')
      if(batchNormalisation){ # Divide each RW block by its norm and the sqrt(number_of_replicates) to avoid overrepresentation of the data in the multi-set.
        normed <- MB[[i]]$Data[batch_position[sorted],] / norm(MB[[i]]$Data[batch_position[sorted],], type = "F")
        growingMB$Data <- normed / sqrt(length(batch_names[[i]]))

      } else {
        growingMB$Data <- MB[[i]]$Data[batch_position[sorted],]
      }

      growingMB$Variables <- MB[[i]]$Variables
      growingMB$Samples <- MB[[i]]$Samples[batch_position[sorted]]
      growingMB$Batch <- rep(batch_names[[i]][j],length(batch_position))

      #if(length(names(MB[[i]])) > 5){ # In case there are additional fields
      if(length(names(MB[[i]])) > 4){ # In case there are additional fields
        newFields <- setdiff(names(MB[[i]]),names(growingMB))

        for(fields in 1:length(newFields)) {
          if(length(MB[[i]][[newFields[fields]]]) == length(MB[[i]]$Samples)) {
            growingMB[[newFields[fields]]] <- MB[[i]][[newFields[fields]]][batch_position[sorted]]
          } else {
            warning(sprintf('The list %s in block %s has an incorrect size.', newFields[fields], as.character(i)))
            print(sprintf('The list %s from block %s was ignored.', newFields[fields],as.character(i)))
          }
        }
      }

      #newMB[[k]] <- growingMB
      newMB[[ paste(names(MB)[i], as.character(batch_names[[i]][j]), sep='_R') ]]<- growingMB

      k <- k + 1

    }
  }

  ## SHOW A DATA-FRAME WITH ALL THE SAMPLE NAMES FOR ALL THE (REPLICATE) BLOCKS.
  df_SampleNames <- matrix(, nrow = length(newMB[[k-1]]$Samples), ncol = length(newMB))
  #labelSampleNames <- as.vector(NULL)

  for (k in 1:length(newMB)){
    df_SampleNames[,k] <- newMB[[k]]$Samples
    #labelSampleNames[k] <- names(newMB[[k]]
  }

  colnames(df_SampleNames) <- names(newMB)
  print('The sample names are:')
  print(df_SampleNames)
  return(newMB)
}
