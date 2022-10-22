# Do the same for regression and discriminant.

ComDim_Exploratory_MB <- function(MB = MB, ndim = NULL,
                                  FUN = FUN,
                                  normalise = FALSE, threshold = 1e-10,
                                  loquace = FALSE,
                                  method = 'FUN',
                                  ...) {

#' ComDim - Finding common dimensions in multi-block datasets.
#'
#' Finding common dimensions in multi-block datasets. Code translated from the following MATLAB function: comdim_PCA_2020.m
#' @param MB A MultiBlock oject.
#' @param ndim Number of Common Dimensions
#' @param FUN The function used as core of the ComDim analysis.
#' @param normalise To apply normalisation. FALSE == no, TRUE == yes (default).
#' @param threshold The threshold limit to stop the iterations. If the "difference of fit" < threshold (1e-10 as default).
#' @param loquace To display the calculation times. TRUE == yes,  FALSE == no (default).
#' @param ... Additional arguments needed for the internal core function of the ComDim analysis.
#' @return A list containing the results from a ComDim analysis.
#' The list fields:
#' Q Global scores (nrow x ndim).
#' P Scaled Global Loadings calculated from Q and the blocks (nvars x ndim).
#' P_Loc Loadings - List containing the Loadings for every block (local vars x ndim).
#' T_Loc Local scores - List containing the (scaled) Local scores for every block (nrow x ndim).
#' saliences Weight of the original blocks in each dimension (ntable x ndim).
#' Sum_saliences_Tab Sum of Saliences for each block in a Dimension.
#' Sum_saliences_Dim Sum of Saliences for each Dimension for a block.
#' b Regression coefficients between Local and Global Scores.
#' explained Percentage explanation given by each dimension (1 x ndim).
#' runtime Period of time spent to execute the analysis (in seconds).
#' NormMB norms of the blocks (1 x 1).
#' MeanMB means of the blocks (1 x nvars).
#' SingVal Vector with singular values (1 x ndim).
#' @examples
#' b1 = matrix(rnorm(500),10,50)
#' batch_b1 = rep(1,10)
#' b2 = matrix(rnorm(800),30,80)
#' batch_b2 = c(rep(1,10),rep(2,10),rep(3,10))
#' # Generate the multi-block (mb)
#' mb <- BuildMultiBlock(b1, batches = batch_b1)
#' mb <- BuildMultiBlock(b2, growingMB = mb, batches = batch_b2, equalSampleNumber = FALSE)
#' rw <- SplitRW(mb)
#' # Do ComDim
#' results <- ComDim_PCA(rw, 2) # In this analysis, we used 2 components.
#' @export

# INITIALISATION

  progress_bar = utils::txtProgressBar(min=0, max=80, style = 3, char = "=")

  start_ini_time <- Sys.time()  # To start counting calculation time for the initialization


  if(class(MB) != 'MultiBlock'){
    stop("'MB' is  not a MultiBlock.")
  }

  ntable = length(getBlockNames(MB)) # number of blocks
  nrowMB = length(getSampleNames(MB)) # number of samples
  variable_number = ncolMultiBlock(MB) # Number of variables per block

  give_error <- 0

  # Remove all variables with 0 variance. If there are no remaining variables, give_error
  numvar <- 0
  numblock0 <- vector()
  for(i in 1:ntable){

    var0 <- apply(MB@Data[[i]], 2, function(x){var(x)})
    var0 <- which(var0 == 0)

    if(length(var0) != 0){
      numvar <- numvar + length(var0)
      if(length(var0) == length(MB@Variables[[i]])){
        numblock0 <- append(numblock0, i)
      }
      MB@Data[[i]] <- MB@Data[[i]][,setdiff(1:variable_number[i],var0)]
      MB@Variables[[i]] <- MB@Variables[[i]][setdiff(1:variable_number[i],var0)]
      variable_number[i] <- length(MB@Variables[[i]])
    }

  }
  if(numvar != 0){
    warning(sprintf('Number of variables excluded from the analysis because of their zero variance: %d', numvar))
  }
  if(length(numblock0) != 0){
    numblock0 <- rev(numblock0) # To sort in decreasing order.
    for(i in 1:length(numblock0)){
      MB@Data[[numblock0[i]]] <- NULL
      MB@Variables[[numblock0[i]]] <- NULL
      if(numblock[[i]] %in% MB@Batch){
        MB@Batch[[numblock0[i]]] <- NULL
      }
      if(numblock[[i]] %in% MB@Metadata){
        MB@Metadata[[numblock0[i]]] <- NULL
      }
    }
    warning(sprintf('Number of blocks excluded from the analysis because of their zero variance: %d', length(numblock0)))
  }

  ntable = length(getBlockNames(MB)) # number of blocks
  variable_number = ncolMultiBlock(MB) # In case the number of blocks has changed.

  if(length(MB@Data) == 0){ # In case all blocks had 0 variance...
    warning('All blocks had 0 variance')
    give_error <- 1
  }

  # Check for infinite or missing values (they cannot be handled with svd)
  for(i in 1:ntable){
    if(any(is.na(MB@Data[[i]]))){
      stop('The MB contains NA values. They must be removed first for ComDim analysis.')
    }
    if(any(is.infinite(MB@Data[[i]]))){
      stop('The MB contains infinite values. They must be removed first for ComDim analysis.')
    }
  }

  if (give_error) {
    stop('The data is not ready for ComDim.')
  } else {
    print('The data can be used for ComDim.')
  }


  if(is.null(ndim)){
    ndim <- ntable # If the number of components is not defined,
                   # the number of components to extract is equal to the number of blocks.
  }

  pieceBar <- 4+2*ntable+ndim # Number of updates in the progress bar.
  pieceBar <- 80/pieceBar
  total_progress <- pieceBar

  DimLabels <- paste0('CC',1:ndim)   # One label per component.
  TableLabels<- getBlockNames(MB) # One label per block.

  end_ini_time <- Sys.time() # To end the count of the calculation time.

  if (loquace) {
    print(sprintf("Initialisation finished after : %s millisecs", (end_ini_time - start_ini_time)*1000))
  }

  utils::setTxtProgressBar(progress_bar, value = total_progress)

# NORMALISATION

  X_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))
  Xnorm_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))

  res_calib <- list()
  temp_tabCalib <- list()
  s_r <- list()
  res_calib$SingVal <-as.vector(NULL)
  res_calib$NormMB <-as.vector(NULL)
  res_calib$MeanMB <-list()

  for (i in 1:ntable) {

    res_calib$MeanMB[[TableLabels[i]]] <- colMeans(MB@Data[[i]])
    names(res_calib$MeanMB[[TableLabels[i]]]) <- MB@Variables[[i]]

    if (normalise){
      # Normalise original blocks

      X_mean <- MB@Data[[i]] - matrix(data = rep(1,nrowMB), ncol = 1, nrow = nrowMB) %*% res_calib$MeanMB[[i]]
      XX <- X_mean*X_mean
      Norm_X <- sqrt(sum(XX))
      X_Normed <- X_mean/Norm_X
      res_calib$NormMB[i] <- Norm_X
      temp_tabCalib[[i]] <- X_Normed
      s_r [[i]] <- X_Normed

    } else {
      res_calib$NormMB[[TableLabels[i]]] <- rep(1,length(MB@Variables[[i]]))
      names(res_calib$NormMB[[TableLabels[i]]]) <- MB@Variables[[i]]
      temp_tabCalib[[i]] <- MB@Data[[i]]
      s_r [[i]] <- MB@Data[[i]]
    }

    if (i==1){
      X_mat[,1:variable_number[1]] = MB@Data[[i]]
      Xnorm_mat[,1:variable_number[1]] = temp_tabCalib[[i]]
    } else {
      beg <- sum(variable_number[1:(i-1)])+1
      ending <- sum(variable_number[1:i])
      X_mat[,beg:ending] = MB@Data[[i]]
      Xnorm_mat[,beg:ending] = temp_tabCalib[[i]]
    }

    norm_comdim <- Sys.time()
    if(loquace){

      print(sprintf("Normalization of block %s finished after : %s millisecs", i, (norm_comdim - start_ini_time)*1000))

    }
    total_progress <- total_progress + pieceBar
    utils::setTxtProgressBar(progress_bar, value = total_progress)

  }

  names(res_calib$NormMB) <- TableLabels

  nR <- nrow(Xnorm_mat)
  nC <- ncol(Xnorm_mat)

  total_progress <- total_progress + pieceBar*ntable
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  ## DO ComDim

  saliences <- matrix(, ncol = ndim, nrow  = ntable)
  Q <- matrix(, ncol = ndim, nrow  = nrowMB)
  unexplained <- ntable # Warning!: true only if the tables were set to norm=1
  varexp <- as.vector(NULL)

  for(dims in 1:ndim){
    previousfit <- 10000
    deltafit <- 1000000
    lambda <- matrix(rep(1,ntable),nrow=ntable,ncol=1)

    qini <- as.vector(s_r[[1]][,1])
    qini <- qini/sqrt(as.vector(qini %*% qini)) #random initialization of t + set to unit length
    qold <- 100

    iters <- 0
    while(norm(qold-qini, "2") > threshold && iters < 100){
      iters <- iters + 1
      qold <- qini
      q <- 0

      W <- matrix(, nrow = nrowMB, ncol = 0)
      for(j in 1:ntable){
        W <- cbind(W,sqrt(lambda[j])*s_r[[j]])
      }

      # HERE IS THE CORE COMDIM
      scores <- FUN(W = W, ndim = ndim,...)

      sv <- sqrt(t(scores[,dims]) %*% scores[,dims]) # For R2X
      qtmp <- scores[,dims] / as.vector(sv)

      for(j in 1:ntable){
        #takes each table into account for lambda after PCA
        lambda[j] <- t(qtmp)%*%(s_r[[j]]%*%t(s_r[[j]]))%*%qtmp
        q <- q + lambda[j]*qtmp
      }

      q <- q/sqrt(as.vector(t(q)%*%q)) #standardizes t
      if (abs(min(q)) > abs(max(q))){
        q <- -q
      }
      qini <- q
    } #deltafit>threshold

    saliences[,dims] <- lambda
    Q[,dims] <- q


    res_calib$SingVal[dims] <- sv^2 # Calculated from the scores (new for Exploratory)
    varexp[dims] <- res_calib$SingVal[dims]^2

    # Deflate blocks
    aux <- diag(nrowMB)-as.matrix(q)%*%t(as.matrix(q))
    for(j in 1:ntable){
      s_r[[j]] <- aux%*%s_r[[j]]
    }

    iter_comdim <- Sys.time()
    if(loquace){

      print(sprintf("Component %s determined after : %s millisecs", dims, (iter_comdim - start_ini_time)*1000))

    }
    total_progress <- total_progress + pieceBar
    utils::setTxtProgressBar(progress_bar, value = total_progress)

    }

  #res_calib$SingVal$i <- 'Singular value'
  names(res_calib$SingVal) <- DimLabels #Dimensions

  #Q$i <- data[[1]]$Samples
  #Q$v <- DimLabels
  colnames(Q) <- DimLabels
  rownames(Q) <- MB@Samples
  res_calib$Q <- Q
  rm(Q)

  ## Adding metadata. Metadata extracted from the first block
  res_calib$metadata <- list()
  if(length(MB@Metadata) != 0){
    res_calib$metadata[[1]] <- MB@Metadata
  }


  end_comdim <- Sys.time()
  if (loquace) {
    print(sprintf("Scores finished after : %s millisecs", (end_comdim - start_ini_time)*1000))
  }
  total_progress <- total_progress + pieceBar
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  meanW <- as.vector(NULL)
  for (j in 1:ntable){
    meanW <- c(meanW,mean(s_r[[j ]]))
  }
  #res_calib$s_n <- s_n
  #res_calib$s_r <- s_r
  rm(s_r)

  ##
  #explained <- list()
  explained <- varexp/sum(varexp)*100
  names(explained) <- DimLabels
  #explained$i <-'% explained'
  #explained$v <- DimLabels #Dimensions
  res_calib$explained <- explained
  rm(varexp)
  rm(explained)

  ## Calculate Sums of saliences for each Dim
  #Sum_saliences_Dim <- list()
  Sum_saliences_Dim <- colSums(saliences)
  names(Sum_saliences_Dim) <- DimLabels
  #Sum_saliences_Dim$i <-'Sum Dim Saliences'
  #Sum_saliences_Dim$v <- DimLabels #Dimensions

  res_calib$Sum_saliences_Dim <- Sum_saliences_Dim
  rm(Sum_saliences_Dim)

  #Calculate Sums of saliences for each Table
  #Sum_saliences_Tab <- list()
  Sum_saliences_Tab <- rowSums(saliences)
  #Sum_saliences_Tab$i <- 'Sum Tab Saliences'
  names(Sum_saliences_Tab) <- TableLabels
  #Sum_saliences_Tab$v <- TableLabels #tables

  res_calib$Sum_saliences_Tab <- Sum_saliences_Tab
  rm(Sum_saliences_Tab)

  #saliences$i <- TableLabels #tables
  #saliences$v <- DimLabels #Dimensions
  res_calib$saliences <- saliences
  colnames(res_calib$saliences) <- DimLabels
  rownames(res_calib$saliences) <- TableLabels

  ##Calculate Normalised concatenated Xs ('Calib') from col
  # Calculate concatenated CD loadings
  # Reorganise Loadings - 1 matrix / LV

  ## If Output = NULL, nothing else is calculated
  #if(!(is.null(output))){
      L_CD_Vec <- NULL
      L_X_Vec <- NULL
      Calib <- NULL

      nCalib <-  nrowMB
  #}

  ## Already defined in during ComDim initiallization
  #isT <- grepl(output, 'T')
  #isL <- grepl(output, 'L')
  #isP <- grepl(output, 'P')
  #isLP <- c(isP,isL)
  #isTLP <- c(isT,isP,isL)

  # tic
  # Calculate concatenated CD loadings
  # Reorganise Loadings - 1 matrix / LV
  #L_CD <- list()
  #L_X <- list()
  T_Loc <- list()
  for(i in 1:ntable){ # Prepare lists
    #L_CD[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = ncol(temp_tabCalib[[i]]))
    #L_X[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = ncol(temp_tabCalib[[i]]))
    T_Loc[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = nrowMB)
  }
  #b <- matrix(,ncol = ndim, nrow = ntable)

  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop('The package pracma is needed.')
  }

  for(j in 1:ndim){
    #T_mat <- matrix(, nrow = nrowMB, ncol = 0)

    for(i in 1:ntable){
    # Q.d are orthonormal in ICA & PCA

      temp <- t(temp_tabCalib[[i]])%*%res_calib$Q[,j] #Scaled CD 'local' Loadings

      #if  (!(is.null(isP))){
        #L_CD[[TableLabels[i]]][,j] <- temp
      #}

      #if  (!(is.null(isL))){
      #  #Q.d are orthonormal in ICA & PCA
      #  L_X[[TableLabels[i]]][,j] <- t(MB@Data[[i]])%*%Q$Data[,j] #Unscaled X 'local' Loadings;
      #}

      #T_mat <- cbind(T_mat,temp_tabCalib[[i]]%*%(temp/as.vector(t(temp)%*%temp))) #local Scores
      #T_mat <- cbind(T_mat,temp_tabCalib[[i]]%*%(temp %*% pracma::pinv(t(temp)%*%temp))) #local Scores

      #if(!(is.null(isT))){
        #T_Loc[[TableLabels[i]]][,j] <- temp_tabCalib[[i]]%*%(temp/as.vector(t(temp)%*%temp)) # local Scores
        T_Loc[[TableLabels[i]]][,j] <- temp_tabCalib[[i]]%*%(temp %*% pracma::pinv(t(temp)%*%temp)) # local Scores
      #}


      # Deflate each temp_tabCalib
      temp_tabCalib[[i]] <- temp_tabCalib[[i]]-res_calib$Q[,j]%*%t(temp)
    }


    # For each CC
    # MLR b-coefficients between Local and Global Scores
    #[b0,b[,j]] <- mlr_DB(T_mat,Q$d[,j],0)
    # b=pinv(X'*X)*X'*Y;
    #b[,j] <- pracma::pinv(t(T_mat)%*%T_mat)%*%t(T_mat)%*%res_calib$Q[,j]
    #b0 <- 0
  }

  for(i in 1:ntable){
    rownames(T_Loc[[ TableLabels[i] ]]) <- rownames(res_calib$Q)
    colnames(T_Loc[[ TableLabels[i] ]]) <- colnames(res_calib$Q)
  }

  # If Output==[], nothing else is calculated
  #if(!(is.null(output))){
    # Calculate Global Loadings
    #if(!(is.null(isP))){
      L_CD_Vec <- t(Xnorm_mat)%*%res_calib$Q # Scaled CD 'global' Loadings
    #}

    #if(!(is.null(isL))){
    #  L_X_Vec <- t(X_mat)%*%Q$Data # Unscaled X 'global' Loadings
    #}
  #}

  #### If Output==[], nothing else is calculated
  load_comdim <- Sys.time()
  if(loquace){
    print(sprintf("Loadings finished after : %s millisecs", (load_comdim - start_ini_time)*1000))
  }

  total_progress <- total_progress + pieceBar
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  rm(X_mat)
  rm(temp_tabCalib)
  rm(temp)

  ## Output
  end_function <- Sys.time()

  #res_calib$b <- b
  #colnames(res_calib$b) <- DimLabels
  #rownames(res_calib$b) <- TableLabels
  #res_calib$b$i <- TableLabels # tables
  #res_calib$b$v <- DimLabels  # Dimensions

  #if(!(is.null(isT))){
    res_calib$T_Loc <- T_Loc
    #res_calib$T_Loc$i <- res_calib$Q$i # samples
    #res_calib$T_Loc$v <- DimLabels # Dimensions
  #}

  # T_Loc is no longer needed
  rm(T_Loc)



  #if(!(is.null(isP))){

    varnames <- vector()
    for(i in 1:ntable){
      varnames <- append(varnames, MB@Variables[[i]])
    }
    #if(length(which(duplicated(varnames))) != 0){
    #  warning('There are at least two variables with the same name. Their name in the rownames of P has been altered.\n')
    #}

    res_calib$P <- L_CD_Vec
    #res_calib$P$i <- varnames # numbers of all variables;
    #res_calib$P$v <- DimLabels # Dimensions
    colnames(res_calib$P) <- DimLabels
    rownames(res_calib$P) <- varnames

    #res_calib$P_Loc <- L_CD
    #res_calib$P_Loc$i <- TableLabels # Tables
    #res_calib$P_Loc$v <- DimLabels # Dimensions

    # for(i in 1:ntable){
    #   rownames(res_calib$P_Loc[[TableLabels[i]]]) <- MB@Variables[[i]]
    #   colnames(res_calib$P_Loc[[TableLabels[i]]]) <- DimLabels
    # }
  #}

  rm(L_CD_Vec)
  #rm(L_CD)

  # if(!(is.null(isL))){
  #   res_calib$Lx$Data <- L_X_Vec
  #   res_calib$Lx$i <- t(c(1:nC)) # numbers of all variables;
  #   res_calib$Lx$v <- DimLabels # Dimensions
  #   colnames(res_calib$Lx$Data) <- DimLabels
  #   rownames(res_calib$Lx$Data) <- t(c(1:nC))
  #   #res_calib$Lx_Loc$i <- TableLabels # Tables
  #   #res_calib$Lx_Loc$v <- DimLabels # Dimensions
  #   res_calib$Lx_Loc$Data <- L_X
  #   for (i in 1:ntable){
  #     colnames(res_calib$Lx_Loc$Data[[TableLabels[i]]]) <- DimLabels
  #     rownames(res_calib$Lx_Loc$Data[[TableLabels[i]]]) <- MB@Variables[[i]]
  #   }
  # }
  # rm(L_X_Vec)
  # rm(L_X)

  ##
  rm(saliences)

  # Define block length
  belong_block <- rep(0, nrow(res_calib$P))
  k <- 1
  for(i in 1:ntable){
    belong_block[k:(k+variable_number[i]-1)] <- TableLabels[i]
    k <- k+variable_number[i]
  }

  #### If Output==[], nothing else is calculated
  #if(normalise == 1){
  #  res_calib$Xnorm_mat$i <- res_calib$Q$i # samples
  #  res_calib$Xnorm_mat$v <- t(c(1:nC)) # numbers of all variables;
  #  res_calib$Xnorm_mat$Data <- Xnorm_mat
  #}
  #rm(Xnorm_mat)
  #### If Output==[], nothing else is calculated

  end_output <- Sys.time()
  running_time <- (end_output - start_ini_time)
  if(loquace){
    print(sprintf("Analysis finished after : %s seconds", running_time))
  }
  res_calib$runtime <- running_time  # Total time of analysis

  progress_bar = utils::txtProgressBar(min=0, max=80, style = 3, char = "=")
  utils::setTxtProgressBar(progress_bar, value = 80)

  close(progress_bar)

  # Save data in a ComDim structure
  res_calib <- new("ComDim",
                   Method = method,
                   ndim = ndim,
                   Q.scores = res_calib$Q,
                   T.scores = res_calib$T_Loc,
                   P.loadings = res_calib$P,
                   Saliences = res_calib$saliences,
                   Orthogonal = list(),
                   R2X = res_calib$explained,
                   R2Y = vector(),
                   Q2 = vector(),
                   DQ2 = vector(),
                   Singular = res_calib$SingVal,
                   Mean = list(MeanMB = res_calib$MeanMB, MeanY = NULL),
                   Norm = list(NormMB = res_calib$NormMB),
                   PLS.model = list(),
                   cv = list(),
                   Prediction = list(),
                   Metadata = res_calib$metadata,
                   variable.block = belong_block,
                   runtime = as.numeric(res_calib$runtime))
  return(res_calib)
}
