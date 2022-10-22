# Do the same for regression and discriminant.

ComDim_y_MB <- function(MB = MB, y = y,
                        ndim = NULL, FUN = FUN,
                        type = c('regression','discriminant')[1],
                        orthogonalization = FALSE,
                        decisionRule = c('fixed','max')[2],
                        normalise = FALSE, threshold = 1e-10,
                        loquace = FALSE,
                        method = 'FUN',
                        ...) {

#' ComDim - Finding common dimensions in multi-block datasets.
#'
#' Finding common dimensions in multi-block datasets. Code translated from the following MATLAB function: comdim_PCA_2020.m
#' @param MB A MultiBlock oject.
#' @param y  The y vector.
#' @param ndim Number of Common Dimensions
#' @param FUN The function used as core of the ComDim analysis.
#' @param type To indicate wether ComDim is used for regression or as a discriminant tool.
#' @param orthogonalization If TRUE, the orthogonal components will be captured.
#' @param decisionRule Only used if method is set to PLS-DA. If 'fixed', samples are assigned to the class...
#' @param normalise To apply normalisation. FALSE == no (default), TRUE == yes.
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

  if(!(type %in% c('regression','discriminant'))){
    stop("'type' must be either 'regression' or 'discriminant'.")
  }

  if(type == 'discriminant' && !(decisionRule %in% c('fixed','max'))){
    stop("'decisionRule' must be either 'fixed' or 'max'.")
  }

  if(orthogonalization && ndim != 1){
    stop("'ndim' must be set to 1 when applying orthogonalization.")
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


  ## Check y matrix (convert to dummy matrix if it is not already.)
  if(type == 'discriminant'){

    tmp <- unique(as.vector(y))
    if(all(tmp %in% c(0,1))){ # Is it dummy?

      if(!is.matrix(y)){
        y <- as.matrix(y)
      }

      classVect <- rep(NA, nrow(y))

      if(ncol(y) == 1){
        classVect <- as.vector(y)
      } else {
        if(any(rowSums(y) != 1)){
          if(any(rowSums(y) == 0)){
            stop('At least one sample was not assigned to a class.')
          } else if (any(colSums(y) > 1)){
            stop('At least one sample was assigned to more than one class.')
          }
        }
      }
      if(is.null(colnames(y))){
        colnames(y) <- paste0("Class",1:ncol(y))
      }
      for(i in 1:ncol(y)){
        classVect[which(y[,i] == 1)] <- colnames(y)[i]
      }
    } else { # If not dummy
      if((is.matrix(y) || is.data.frame(y)) && ncol(y) > 1){
        stop("Predictive analysis can only be performed if 'y' is a dummy matrix or a class vector.")
      }
      classVect <- as.vector(y)
      classVect_sorted <- sort(unique(classVect))

      # Now generate the dummy matrix from y
      y <- matrix(rep(0, length(classVect_sorted)*length(classVect)),
                  ncol = length(classVect_sorted), nrow =length(classVect))
      for(i in 1:length(classVect_sorted)){
        y[which(classVect == classVect_sorted[i]),i] <- 1
      }
      colnames(y) <- classVect_sorted
      rm(classVect_sorted)
    }
  }

  if(!is.matrix(y)){ # If it's a vector and the method is PLS-R
    y <- as.matrix(y)
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

    # PLS addition
    meanX <- colMeans(Xnorm_mat)
    meanY <- colMeans(y)
    # End PLS addition

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

  # PLS addition
  ICs <- ndim
  # End of PLS addition

  ## DO ComDim

  saliences <- matrix(, ncol = ndim, nrow  = ntable)
  Q <- matrix(, ncol = ndim, nrow  = nrowMB)
  if(orthogonalization){
    saliences.ort <- saliences # The same dimensions
    Q.ort <- Q # The same dimensions
  }
  unexplained <- ntable # Warning!: true only if the tables were set to norm=1
  varexp <- as.vector(NULL)
  varexp.y <- as.vector(NULL)
  PLS2_Scores <- matrix(, ncol = ndim, nrow  = nrowMB)
  PLS2_P <- matrix(, ncol = ndim, nrow  = sum(variable_number))
  PLS2_W <- matrix(, ncol = ndim, nrow  = sum(variable_number))
  PLS2_Q <- matrix(, ncol = ndim, nrow  = ncol(as.matrix(y)))
  PLS2_U <- matrix(, ncol = ndim, nrow  = nrowMB)

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
      output <- FUN(W = W, y = y, ndim = ndim,...)
      PLS2_Scores[,dims] <- output$scores
      PLS2_P[,dims] <- output$P
      PLS2_W[,dims] <- output$W
      PLS2_Q[,dims] <- output$Q
      PLS2_U[,dims] <- output$U

      sv <- sqrt(t(PLS2_Scores[,dims]) %*% PLS2_Scores[,dims]) # For R2X
      varexp.y[dims] <- (t(PLS2_U[,dims]) %*% PLS2_U[,dims])^2 # For R2Y (see varexp)
      qtmp <- PLS2_Scores[,dims] / as.vector(sv)

      if(orthogonalization && output$ort != 0){
        varexp.ort <- sv.ort <- rep(NA, output$ort)
        q.ort <- qtmp.ort <- matrix(, nrow = length(qtmp), ncol = output$ort)
        lambda.ort <- matrix(, nrow = ntable, ncol = output$ort)
        for(j in 1:output$ort){
          sv.ort[j] <- sqrt(t(output$orthoscores[,j]) %*% output$orthoscores[,j]) # For R2X
          varexp.ort[j] <- sv.ort[j]^4
          qtmp.ort[,j] <- output$orthoscores[,j] / as.vector(sv.ort[j])
          q.ort[,j] <- rep(0,output$ort) # To fill afterwards
        }
      }
      #print(sv.ort) vector of 3
      #print(qtmp.ort) matrix of 12x3

      for(j in 1:ntable){
        #takes each table into account for lambda after PCA
        lambda[j] <- t(qtmp)%*%(s_r[[j]]%*%t(s_r[[j]]))%*%qtmp
        q <- q + lambda[j]*qtmp

        if(orthogonalization && output$ort != 0){
          for(jj in 1:output$ort){
            lambda.ort[j,jj] <- t(qtmp.ort[,jj])%*%(s_r[[j]]%*%t(s_r[[j]]))%*%qtmp.ort[,jj]
            q.ort[,jj] <- q.ort[,jj] + lambda.ort[j,jj]*qtmp.ort[,jj]
          }
        }
      }

      q <- q/sqrt(as.vector(t(q)%*%q)) #standardizes t
      if (abs(min(q)) > abs(max(q))){
        q <- -q
      }
      if(orthogonalization && output$ort != 0){
        for(jj in 1:output$ort){
          q.ort[,jj] <- q.ort[,jj]/sqrt(as.vector(t(q.ort[,jj])%*%q.ort[,jj])) #standardizes t
          if (abs(min(q.ort[,jj])) > abs(max(q.ort[,jj]))){
            q.ort[,jj] <- -q.ort[,jj]
          }
        }
      }
      qini <- q
    } #deltafit>threshold

    saliences[,dims] <- lambda
    Q[,dims] <- q
    saliences.ort <- Q.ort <- NULL
    if(orthogonalization && output$ort != 0){
      saliences.ort <- lambda.ort
      Q.ort <- q.ort
      colnames(Q.ort) <- paste0('CC',1:output$ort)
      rownames(Q.ort) <- MB@Samples
    }


    res_calib$SingVal[dims] <- sv^2 # Calculated from the scores (new for Exploratory)
    varexp[dims] <- res_calib$SingVal[dims]^2

    # Deflate blocks
    aux <- diag(nrowMB)-as.matrix(q)%*%t(as.matrix(q))
    if(orthogonalization && output$ort != 0){
      aux <- aux-as.matrix(q.ort)%*%t(as.matrix(q.ort))
    }
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
  rownames(output$y) <- NULL
  colnames(output$y) <- NULL

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
    meanW <- c(meanW,colMeans(s_r[[j ]]))
  }
  #res_calib$s_n <- s_n
  #res_calib$s_r <- s_r
  rm(s_r)

  ##
  #explained <- list()
  explained <- varexp/sum(varexp)*100
  explained.ort <- NULL
  if(orthogonalization && output$ort != 0){
    total.varexp <- varexp + sum(varexp.ort)
    explained <- 100*varexp/total.varexp
    explained.ort <- 100*varexp.ort/total.varexp
    names(explained.ort) <- colnames(Q.ort)
  }
  r2y <- 100*varexp.y/sum(varexp.y)
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
  L_X <- list()
  T_Loc_o <- T_Loc <- list()
  for(i in 1:ntable){ # Prepare lists
    #L_CD[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = ncol(temp_tabCalib[[i]]))
    #L_X[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = ncol(temp_tabCalib[[i]]))
    T_Loc[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = nrowMB)
    if(orthogonalization && output$ort != 0){
      T_Loc_o[[TableLabels[i]]] <- matrix(,ncol = output$ort, nrow = nrowMB)
    }
  }
  b <- matrix(,ncol = ndim, nrow = ntable)

  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop('The package pracma is needed.')
  }
  orthoP <- NULL
  if(orthogonalization && output$ort != 0){
    for(j in 1:output$ort){
      for(i in 1:ntable){
        tempo <- t(temp_tabCalib[[i]])%*% Q.ort[,j] #Scaled CD 'local' Loadings
        T_Loc_o[[TableLabels[i]]][,j] <- temp_tabCalib[[i]]%*%(tempo %*% pracma::pinv(t(tempo)%*%tempo)) # local Scores
        # Deflate each temp_tabCalib
        temp_tabCalib[[i]] <- temp_tabCalib[[i]]- Q.ort[,j]%*%t(tempo)
        if(i == 1){
          tempo2 <- tempo
        } else {
          tempo2 <- append(tempo2,tempo)
        }
      }
      if(j == 1){
        orthoP <- tempo2
      } else {
        orthoP <- cbind(orthoP, tempo2)
      }
    }
    orthoP <- as.matrix(orthoP)
    colnames(orthoP) <- paste0('CC',1:output$ort)
    for(i in 1:ntable){
     if(i == 1){
       Xnorm_mat <- temp_tabCalib[[i]]
     } else {
       Xnorm_mat <- cbind(Xnorm_mat, temp_tabCalib[[i]])
     }
    }
  }




  for(j in 1:ndim){
    T_mat <- matrix(, nrow = nrowMB, ncol = 0)

    for(i in 1:ntable){

      temp <- t(temp_tabCalib[[i]])%*%res_calib$Q[,j] #Scaled CD 'local' Loadings

      #if  (!(is.null(isP))){
      #L_CD[[TableLabels[i]]][,j] <- temp
      #}

      #if  (!(is.null(isL))){
      #  L_X[[TableLabels[i]]][,j] <- t(data[[i]]$Data)%*%res_calib$Q[,j] #Unscaled X 'local' Loadings;
      #}

      #T_mat <- cbind(T_mat,temp_tabCalib[[i]]%*%(temp/as.vector(t(temp)%*%temp))) #local Scores
      T_mat <- cbind(T_mat,temp_tabCalib[[i]]%*%(temp %*% pracma::pinv(t(temp)%*%temp))) #local Scores

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
    b[,j] <- pracma::pinv(t(T_mat)%*%T_mat)%*%t(T_mat)%*%res_calib$Q[,j]

    #b0 <- 0
  }

  for(i in 1:ntable){
    rownames(T_Loc[[ TableLabels[i] ]]) <- rownames(res_calib$Q)
    colnames(T_Loc[[ TableLabels[i] ]]) <- colnames(res_calib$Q)
    if(orthogonalization && output$ort != 0){
      rownames(T_Loc_o[[TableLabels[i]]]) <- rownames(Q.ort)
      colnames(T_Loc_o[[TableLabels[i]]]) <- colnames(Q.ort)
    }
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

  ## Addition for PLS
  if(identical(unlist(output$y[,1]),
               unlist((y[,1]-mean(y[,1]))/sd(y[,1])))){ # Check if data in FUN is mean-centered.
    ComDimPLS2_B <- PLS2_W %*% pracma::pinv(t(PLS2_P) %*% PLS2_W) %*% t(PLS2_Q*sd(y)+mean(y))
  } else {
    ComDimPLS2_B <- PLS2_W %*% pracma::pinv(t(PLS2_P) %*% PLS2_W) %*% t(PLS2_Q)
  }
  ComDimPLS2_B0 <- meanY - meanW %*% ComDimPLS2_B

  # End of addition for PLS

  res_calib$b <- b
  colnames(res_calib$b) <- DimLabels
  rownames(res_calib$b) <- TableLabels
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


  ## 07/08/2022 Calculate Y predictive and rescale.
  Ypred <- Xnorm_mat %*% ComDimPLS2_B

  for(i in 1:ncol(y)){
    Ypred[,i] <- unlist(Ypred[,i]) + ComDimPLS2_B0[i]
  }


  ## 07/08/2022 Predict classes
  if(type == 'discriminant'){

    predClass <- rep(NA, nrow(y))

    if(decisionRule == 'max'){

      for(i in 1:nrow(Ypred)){
        xx <- which(Ypred[i,] == max(Ypred[i,], na.rm = TRUE))
        if(length(xx) == 1) {
          predClass[i] <- colnames(y)[xx] # Y stores the Class names
        } else {
          predClass[i] <- NaN # There is a tie.
          warning(sprintf('Class for sample %s could not be predicted.', i))
        }
      }

    } else if(decisionRule == 'fixed'){
      for(i in 1:nrow(Ypred)){
        xx <- which(Ypred[,i] > 1/ncol(y))
        if(length(xx) == 1) {
          predClass[i] <- colnames(y)[xx] # Y stores the Class names
        } else {
          predClass[i] <- NaN # There is a tie.
          warning(sprintf('Class for sample %s could not be predicted.', i))
        }
      }
    }

    if(is.matrix(y)){
      meanY <-colMeans(y)
    } else {
      meanY <- mean(y)
    }

    classy <- colnames(y)
    Q2 <- DQ2 <- specvec <- sensvec <- rep(NA, length(classy))
    confusionMatrix <- list()

    for (i in 1:length(classy)){

      pos <- which(classVect == classy[i])
      neg <- which(classVect != classy[i])
      tp <- length(which(classVect == classy[i] & predClass == classy[i]))
      tn <- length(which(classVect != classy[i] & predClass != classy[i]))
      fp <- length(which(classVect != classy[i] & predClass == classy[i]))
      fn <- length(which(classVect == classy[i] & predClass != classy[i]))

      if(length(tp) == 0){
        tp <- 0
      }

      if(length(tn) == 0){
        tn <- 0
      }

      if(length(fp) == 0){
        fp <- 0
      }

      if(length(fn) == 0){
        fn <- 0
      }

      sensvec[i] <- tp/length(pos)
      specvec[i] <- tn/length(neg)

      cm <- matrix(c(tp,fp,fn,tn), ncol = 2, nrow = 2)
      rownames(cm) <- c("predClass1", "predClass0")
      colnames(cm) <- c("trueClass1", "trueClass0")
      confusionMatrix[[classy[i]]] <- cm

      ## DQ2
      E0 <- Ypred[neg,i]-y[neg,i]   # Calculate Residuals of Class0 samples
      E1 <- Ypred[pos,i]-y[pos,i]   # Calculate Residuals of Class1 samples

      E0count <- which(E0 > 0)           # Find predictions for Class0 samples larger than 0
      E1count <- which(E1 < 0)           # Find predictions for Class1 samples larger than 1

      SSE0_all <- t(E0) %*% E0 # Calculate SSE for those samples of Class0
      SSE1_all <- t(E1) %*% E1 # Calculate SSE for those samples of Class1

      SSE0 <- t(E0[E0count]) %*% E0[E0count] # Calculate SSE for those samples of Class0
      SSE1 <- t(E1[E1count]) %*% E1[E1count] # Calculate SSE for those samples of Class1

      PRESSD <- SSE0 + SSE1           # Calculate total SSE for all samples = PRESSD
      PRESS <- SSE0_all + SSE1_all

      Ym <- y[,i] - mean(y[,i])               # Calculate Total sum of squares (over all samples)
      TSS <- t(Ym) %*% Ym

      DQ2[i] <- 1 - (PRESSD/TSS)         # Calculate DQ2
      Q2[i] <- 1 - PRESS/TSS

    }
    names(Q2) <- classy
    names(DQ2) <- classy
  } else if(type == 'regression'){
    PRESS <- t(Ypred-y) %*% (Ypred-y)
    Ym <- y - mean(y)
    TSS <- t(Ym) %*% Ym
    Q2 <- 1 - PRESS/TSS
  }

  #Ypred <-Ypred*meanY # Re-scale Ypred. I think this is not needed.




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
  # Orthogonal = list(
  #   nort = res_calib$cv$OrthoLVsOptimalNum,
  #   Q.scores.ort = Q.scores.ort,
  #   T.scores.ort = T_Loc_o,
  #   P.loadings.ort = loadings.ort,
  #   Saliences.ort = Saliences.ort,
  #   #Y.ort = res_calib$ComDimOPLSModel$Up,
  #   R2X = if(res_calib$cv$OrthoLVsOptimalNum == 0){
  #     NULL
  #   } else {
  #     res_calib$ComDimOPLSModel$R2XC[2:(res_calib$cv$OrthoLVsOptimalNum+1)] #Explained variation for predictive model components after addition of Y-orthogonal model components.
  #   }
  # ),
  #

  res_calib <- new("ComDim",
                   Method = method,
                   ndim = ndim,
                   Q.scores = res_calib$Q,
                   T.scores = res_calib$T_Loc,
                   P.loadings = res_calib$P,
                   Saliences = res_calib$saliences,
                   Orthogonal = list(nort = output$ort,
                                     Q.scores.ort = Q.ort,
                                     T.scores.ort = T_Loc_o,
                                     P.loadings.ort = orthoP,
                                     Saliences.ort = saliences.ort,
                                     R2X = explained.ort),
                   R2X = res_calib$explained,
                   R2Y = r2y,
                   Q2 = Q2,
                   DQ2 = if(type == 'discriminant'){
                      DQ2
                   } else {
                      vector()
                   },
                   Singular = res_calib$SingVal,
                   Mean = list(MeanMB = res_calib$MeanMB, MeanY = meanY),
                   Norm = list(NormMB = res_calib$NormMB),
                   PLS.model = list(W = PLS2_W,
                                    B = ComDimPLS2_B,
                                    B0 = ComDimPLS2_B0,
                                    Y = y), # W, U, B, B0, Y
                   cv = list(),
                   Prediction = if(type == 'discriminant'){
                     list(Y.pred = Ypred,
                          decisionRule = decisionRule,
                          trueClass = classVect,
                          predClass = predClass,
                          Sensitivity = sensvec,
                          Specificity = specvec,
                          confusionMatrix = confusionMatrix)
                   }else {
                     list(Y.pred = Ypred)
                   },
                   Metadata = res_calib$metadata,
                   variable.block = belong_block,
                   runtime = as.numeric(res_calib$runtime))

  if(!orthogonalization){
    res_calib@Orthogonal <- list()
  }
  return(res_calib)
}
