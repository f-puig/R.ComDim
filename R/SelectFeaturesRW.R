SelectFeaturesRW <- function(RW = RW, results = results, ndim = NULL, blocks = NULL, threshold_cor = 1, threshold_cov = 1, mean.RW = T, plots = "NO") {

#' SelectFeaturesRW
#'
#' Finds the important variables presenting a coordinated response in all the replicate-blocks
#' @param RW The object used as input in the ComDim analysis.
#' @param results The output object obtained in the ComDim analysis.
#' @param ndim The number of the component for which the most important variables are to be calculated.
#' @param blocks A vector with the indices or the names for the replicate blocks of the same data type.
#' @param threshold_cor The "times" parameter used to calculate the threshold in the following formula: cor(variable) > times * sd(cor(variables)). Minimal value that can be assigned to threshold_cor is 1.
#' @param threshold_cov The "times" parameter used to calculate the threshold in the following formula: cov(variable) > times * sd(cov(variables)). Minimal value that can be assigned to threshold_cor is 1.
#' @param mean.RW Logical value to indicate whether the RW data must be mean-centered (TRUE) or not (FALSE).
#' @param plots Parameter to indicate whether a plot must be produced. Possible values are "NO" for no plots, "separated" for plotting each plot individually, and "together" to plot all the plots in the same grid.
#' @return An object with 2 lists. The first list contains the important variables presenting a positive relationship with the scores, while the second list contains the most important variables presenting a negative relationship.
#' @examples
#'  b1 = matrix(rnorm(500),10,50)
#'  batch_b1 = rep(1,10)
#'  b2 = matrix(rnorm(800),30,80)
#'  batch_b2 = c(rep(1,10),rep(2,10),rep(3,10))
#'  mb <- BuildMultiBlock(b1, batches = batch_b1)
#'  mb <- BuildMultiBlock(b2, growingMB = mb, batches = batch_b2, equalSampleNumber = FALSE)
#'  rw <- SplitRW(mb)
#'  results<-ComDim_PCA(rw, 2) # In this analysis, we used 2 components.
#'  features <- SelectFeaturesRW(RW = rw, results = results, ndim = 1, blocks = c(2,3,4))
#' @export

  # Create Covariance and Correlation variables to avoid warnings during package checking.
  Covariance <- NULL
  Correlation <- NULL

  # Check data
  if(class(RW) != 'MultiBlock'){
    stop("'RW' is  not a MultiBlock.")
  }

  if (is.list(RW@Samples)){
    stop('The replicate blocks to evaluate have a different sample size.')
  }
  sample_number <- length(getSampleNames(RW))

  if (length(unique(lengths(getVariableNames(RW, blocks)))) != 1){
    stop('The replicate blocks to evaluate have a different variable size.')
  }

  if (!("ComDim" %in% class(results))){
    stop("'results' is not the output from a ComDim analysis.")
  }

  if (ndim > results@ndim){
    stop("The 'ndim' value chosen for this analysis is too high.")
  }

  if(nrow(results@Q.scores) != length(getSampleNames(RW))){
    stop("The size of the local scores do not match the size of 'RW'.")
  }

  # Check if result contains the blocks of local loadings
  if(length(results@variable.block) != sum(lengths(RW@Variables))){
    stop("The total variable number of the ComDim model does not match the total variable number of 'RW'.")
  }

  if (!(is.logical(mean.RW))){
    stop("The parameter 'mean.RW' must be logical.")
  }

  if(!(is.character(plots))){
    if(is.logical(plots)){
      if (plots == FALSE){
        plots <- "NO"
      } else {
        plots <- "separated"
      }
    }
  }

  if(!(requireNamespace("pracma", quietly = TRUE))){
    warning("The 'pracma' package is needed but not installed." )
  }

  if(plots != "NO" && !(requireNamespace("ggplot2", quietly = TRUE))){
    warning("The 'ggplot2' package is needed but not installed." )
    plots <- "NO"
  }

  if(plots != "NO" && !(requireNamespace("gridExtra", quietly = TRUE))){
    warning("The 'gridExtra' package is needed but not installed." )
    plots <- "separated"
  }

  mean_blocks <- list(NULL)
  for (i in blocks){
    # Mean-center
    if (isTRUE(mean.RW) && requireNamespace("pracma", quietly = TRUE)){
      mean_blocks[[i]] <- RW@Data[[ i ]] - pracma::repmat(colMeans(RW@Data[[ i ]]),unique(sample_number),1)
    } else {
      mean_blocks[[i]] <- RW@Data[[ i ]]
    }
  }

  if(length(ndim) != 1){
    stop("The length of 'ndim' must be of 1.")
  }

  if(!(is.numeric(ndim))){
    stop("'ndim' must be numeric.")
  }

  s2 <- s1 <- matrix(, nrow = length(blocks), ncol=length(getVariableNames(RW, blocks[1])))

  # Scale Q in case it wasnt (if the data comes from ComDim, it was).
  if (requireNamespace("pracma", quietly = TRUE)){
    qs <- results@Q.scores[,ndim]%*%pracma::pinv(t(results@Q.scores[,ndim])%*%results@Q.scores[,ndim])
  } else {
    qs <- results@Q.scores[,ndim]
    warning("If 'results' was not scaled, the output from this function may be incorrect.")
  }

  # Calculate the corr and cov values.
  for (i in 1:length(blocks)){
    # Calculate the s-plot values.
    #s1(i,:) = (T' * X_Cent) ./(N - 1); % cov
    #s2(i,:) = s1(i,:) ./ (std(T) .* std(X_Cent)); %cor
    s1[i,] <- (t(qs) %*% mean_blocks[[ blocks[i] ]])/(unique(sample_number)-1)
    s2[i,] <- s1[i,]/(stats::sd(qs)*apply(mean_blocks[[ blocks[i] ]], 2, stats::sd))
  }

  s1[is.nan(s1)] <- 0
  s2[is.nan(s2)] <- 0
  s1_imp <- list(NULL)
  s2_imp <- list(NULL)
  for (i in 1:length(blocks)){
    if (is.numeric(threshold_cov) && length(threshold_cov) == 1){
      if (threshold_cov < 1) {
        threshold_cov <- 1
        warning("'threshold_cov value has been modified to 1.")
      }
      std_cov <- stats::sd(s1[i,])*threshold_cov
      s1_imp[[i]] <- c(which(s1[i,] > std_cov), which(s1[i,] < -std_cov)) # Covariate and anticovariate
    } else {
      stop("'Threshold_cov' must be a numeric object of length 1.")
    }

    if (is.numeric(threshold_cor) && length(threshold_cor) == 1){
      if (threshold_cor < 1) {
        threshold_cor <- 1
        warning("'threshold_cor value has been modified to 1.")
      }
      std_cor <- stats::sd(s2[i,])*threshold_cor
      s2_imp[[i]] <- c(which(s2[i,] > std_cor), which(s2[i,] < -std_cor))
    } else {
      stop("'Threshold_cor' must be a numeric object of length 1.")
    }
  }

  if(is.numeric(blocks)){
    blocks <- getBlockNames(RW)[blocks]
  }

  names(s1_imp) <- blocks
  names(s2_imp) <- blocks
  imp_pos <- list()
  imp_neg <- list()

  for(i in blocks){
    imp <- intersect(s1_imp[[i]],s2_imp[[i]])
    # I only use the loadings to check their sign.
    imp_pos[[i]] <- intersect(imp, which(results@P.loadings[which(results@variable.block == i),ndim] > 0))
    imp_neg[[i]] <- intersect(imp, which(results@P.loadings[which(results@variable.block == i),ndim] < 0))
  }

  common_RW_pos <- list()
  common_RW_neg <- list()
  for(i in 1:length(blocks)){
    if (i == 1){
      common_RW_pos <- imp_pos[[i]]
      common_RW_neg <- imp_neg[[i]]
    } else {
      common_RW_pos <- intersect(common_RW_pos, imp_pos[[i]])
      common_RW_neg <- intersect(common_RW_neg, imp_neg[[i]])
    }
  }

  names(common_RW_pos) <- RW@Variables[[blocks[1]]][common_RW_pos]
  names(common_RW_neg) <- RW@Variables[[blocks[1]]][common_RW_neg]
  ivs <- list (positive = common_RW_pos, negative = common_RW_neg)

  if (grepl('NO', plots, ignore.case = TRUE)) {
    # Do nothing.
  }else if (grepl('separated', plots, ignore.case = TRUE) && (requireNamespace("ggplot2", quietly = TRUE))){
    for(i in 1:length(blocks)){
      xx <- data.frame(Covariance = s1[i,], Correlation = s2[i,])
      color<-rep('#000000', length(getVariableNames(RW, blocks[i])))
      color[c(imp_pos[[i]],imp_neg[[i]])] <- rep('#FF0000', length(c(imp_pos[[i]],imp_neg[[i]])))
      labels <- rep("",length(getVariableNames(RW, blocks[i])))
      labels[c(imp_pos[[i]],imp_neg[[i]])] <- RW@Variables[[blocks[i]]][c(imp_pos[[i]],imp_neg[[i]])]
      useless<-readline(prompt=sprintf("S-plot for block %d shown. Press [enter] to continue", blocks[i]))
      print(ggplot2::ggplot(xx, ggplot2::aes(x= Covariance, y= Correlation)) +
        ggplot2::geom_point(color = color) +
        ggplot2::geom_text(label=labels) +
        ggplot2::ggtitle(sprintf("S-plot for Block %s and component %i", blocks[i], ndim)))
    }
  }else if (grepl('together', plots, ignore.case = TRUE)  && (requireNamespace("ggplot2", quietly = TRUE))){
    pi <- list(NULL)
    colori <- list(NULL)
    labelsi <- list(NULL)
    xx <- list(NULL)
    for(i in 1:length(blocks)){
      colori[[i]] <-rep('#000000', length(getVariableNames(RW, blocks[i])))
      colori[[i]][c(imp_pos[[i]],imp_neg[[i]])] <- rep('#FF0000', length(c(imp_pos[[i]],imp_neg[[i]])))
      xx[[i]] <- data.frame(Covariance = s1[i,], Correlation = s2[i,])
      labelsi[[i]] <- rep("",length(getVariableNames(RW, blocks[i])))
      labelsi[[i]][c(imp_pos[[i]],imp_neg[[i]])] <- RW@Variables[[blocks[i]]][c(imp_pos[[i]],imp_neg[[i]])]
      pi[[i]] <- ggplot2::ggplot(xx[[i]], ggplot2::aes(x= Covariance, y= Correlation)) +
        ggplot2::geom_point(color = colori[[i]]) +
        ggplot2::geom_text(label=labelsi[[i]]) +
        ggplot2::ggtitle(sprintf("S-plot for Block %s and component %i", blocks[i], ndim))
    }
    to_eval <- "gridExtra::grid.arrange("
    for(i in 1:length(blocks)){
      to_eval <- paste0(to_eval, sprintf("pi[[%d]], ",i))
    }
    to_eval <-paste0(to_eval, " nrow = 1)")
    eval(parse(text = to_eval))
  }else{
    warning("The 'plots' parameter was not set correctly. Plots will not be displayed.")
  }

  return(ivs) # ivs = important variables
}
