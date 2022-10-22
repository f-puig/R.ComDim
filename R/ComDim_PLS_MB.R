# For the DA situation I should include the formula of the VIPs.
# It requires pracma::inv, which is invoked in the package.
# Need to change all matrix divisions by pracma:inv * ...
# Need to check that the result is the same than in Matlab.
# Check how the saliences are computed.
# See if the B values make sense (on June DNR was skeptical on how to deal with this situation in the multi-block scenareo due to the different scaling...)


ComDim_PLS_MB <- function(MB = MB, y = y, ndim = NULL,
                          method = c('PLS-DA','PLS-R'),
                          decisionRule = c('fixed','max')[2],
                          normalise = FALSE, threshold = 1e-10,
                          loquace = FALSE,
                          CompMethod = 'Normal', Partitions = 1) {

#' ComDim - Finding common dimensions in multi-block datasets
#'
#' Finding common dimensions in multi-block datasets. Code translated from the following MATLAB function: comdim_PCA_2020.m
#' @param MB A MultiBlock object.
#' @param y The Y-block to use in the PLS model as dependent data. A class vector or a dummy matrix.
#' @param method PLS-DA or PLS-R
#' @param decisionRule Only used if method is set to PLS-DA. If 'fixed', samples are assigned to the class...
#' @param ndim Number of Common Dimensions
#' @param normalise To apply normalisation. FALSE == no, TRUE == yes (default).
#' @param threshold The threshold limit to stop the iterations. If the "difference of fit" < threshold (1e-10 as default).
#' @param loquace To display the calculation times. TRUE == yes,  FALSE == no (default).
#' @param CompMethod To speed up the analysis for really big multi-blocks. 'Normal' (default), 'Kernel', 'PCT', 'Tall' or 'Wide'.
#' @param Partitions To speed up the analysis for really big multi-blocks. This parameter is used if CompMethod is 'Tall' or 'Wide'.
#' @return A list containing the results from a ComDim analysis.
#' The list fields:
#' Q Global scores (nrow x ndim).
#' P Scaled Global Loadings calculated from Q and the blocks (nvars x ndim).
#' P_Loc Loadings - List containing the Loadings for every block (local vars x ndim).
#' T_Loc Local scores - List containing the (scaled) Local scores for every block (nrow x ndim).
#' Lx Unscaled Global loadings for data.
#' Lx_Loc Unscaled Local loadings for data.
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
#' # Generate the Y-block
#' y <- scale(1:10, center = TRUE)
#' # Do ComDim
#' results <- ComDim_PLS(rw, y, 2) # In this analysis, we used 2 components.
#' @export

# INITIALISATION

  if(!(method == 'PLS-DA' || method == 'PLS-R')){
    stop("'method' must be 'PLS-DA' or 'PLS-R'.")
  }

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

  ## Check y matrix (convert to dummy matrix if it is not already.)
  if(method == 'PLS-DA'){

    tmp <- unique(as.vector(y))
    if(all(tmp %in% c(0,1))){ # Is it dummy?

      if(!is.matrix(y)){ # If it's a vector and the method is PLS-R
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
        stop('ComDim-PLS can only be applied if Y is a dummy matrix or a class vector')
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

  #yraw <- y # The y without normalization

  pieceBar <- 4+2*ntable+ndim # Number of updates in the progress bar.
  pieceBar <- 80/pieceBar
  total_progress <- pieceBar

  DimLabels <- paste0('CC',1:ndim)   # One label per component.
  TableLabels <- getBlockNames(MB) # One label per block.

  #isT <- grepl(output, 'T')
  #isL <- grepl(output, 'L')
  #isP <- grepl(output, 'P')

  if(CompMethod == 'Tall' && any(variable_number > 1000)){
    print("The method 'Tall' is not recommeded for large blocks (> 1000 columns).")
    print("Do you still want to continue with the analysis?")
    print("You can change now the CompMethod if you prefer: (yes/no)")
    ans1 <- tolower(readline("Type your choice: (y/n)"))
    if (ans1 == 'y' || ans1 == 'yes') {
      print("Which method do you want to use instead (Normal/Kernel/PCT/Wide)?")
      ans2 <- tolower(readline("Type your choice: (n/k/p/w)"))
      if (ans2 == 'n' || ans2 == 'normal'){
        CompMethod = 'Normal'
      } else if(ans2 == 'k' || ans2 == 'kernel'){
        CompMethod = 'Kernel'
      } else if(ans2 == 'p' || ans2 == 'pct'){
        CompMethod = 'PCT'
      } else if(ans2 == 'w' || ans2 == 'wide'){
        CompMethod = 'Wide'
      } else {
        stop("Wrong answer typed.")
      }
    } else if (ans1 == 'n' || ans1 == 'no'){
      print("The method 'Tall' will be used")
    } else {
      stop("Wrong answer typed.")
    }
  }

  end_ini_time <- Sys.time() # To end the count of the calculation time.

  if (loquace) {
    print(sprintf("Initialisation finished after : %s millisecs", (end_ini_time - start_ini_time)*1000))
  }

  utils::setTxtProgressBar(progress_bar, value = total_progress)

  # NORMALISATION

  #start_norm_time <- Sys.time()  # To start counting calculation time for the initialization

  X_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))
  Xnorm_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))

  res_calib <- list()
  temp_tabCalib <- list()
  s_n <- list()
  res_calib$SingVal <-as.vector(NULL)
  res_calib$NormMB <-as.vector(NULL)
  res_calib$MeanMB <-list()

  for (i in 1:ntable) {

    res_calib$MeanMB[[TableLabels[i]]] <- colMeans(MB@Data[[i]])
    names(res_calib$MeanMB[[TableLabels[i]]]) <- MB@Variables[[i]]
    #res_calib$MeanMB$i <- TableLabels

    if (normalise){
      # Normalise original blocks

      X_mean <- MB@Data[[i]] - matrix(data = rep(1,nrowMB), ncol = 1, nrow = nrowMB) %*% res_calib$MeanMB[[i]]
      # In a previous version (that worked), I had used as.matrix instead of matrix.
      XX <- X_mean*X_mean

      Norm_X <- sqrt(sum(XX))
      X_Normed <- X_mean/Norm_X

      res_calib$NormMB[i] <- Norm_X

      temp_tabCalib[[i]] <- X_Normed
      s_n [[i]] <- X_Normed

    } else {
      res_calib$NormMB[[TableLabels[i]]] <- rep(1,length(MB@Variables[[i]]))
      names(res_calib$NormMB[[TableLabels[i]]]) <- MB@Variables[[i]]
      #res_calib$MeanMB$Variables <- MB@Variables[[i]]
      #res_calib$MeanMB$i <- TableLabels

      temp_tabCalib[[i]] <- MB@Data[[i]]
      s_n [[i]] <- MB@Data[[i]]
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

  #res_calib$NormMB$i <- 'Norm'
  #names(res_calib$NormMB$Data) <- TableLabels

  nR <- nrow(Xnorm_mat)
  nC <- ncol(Xnorm_mat)

# COMPRESSION

  s_r <- Compress_Data_2020(s_n, CompMethod, Partitions)

  end_comp_time <- Sys.time()

  if(loquace){

    print(sprintf("Compression finished after : %s millisecs", (end_comp_time - start_ini_time)*1000))

  }

  total_progress <- total_progress + pieceBar*ntable
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  # PLS addition
  ICs <- ndim
  # End of PLS addition

  ## DO ComDim
  #ini_comdim <- Sys.time()

  #saliences <- list()
  saliences <- matrix(, ncol = ndim, nrow  = ntable)
  #Q <- list()
  Q <- matrix(, ncol = ndim, nrow  = nrowMB)
  unexplained <- ntable; # Warning!: true only if the tables were set to norm=1
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

      # usv <-svd(W) This was for ComDim_PCA

      # Addition for PLS
      PLS2_Res <- PLS2_DNR_2018(X = W,Y = y, lvs = ICs)

      PLS2_Scores[,dims] <- PLS2_Res$Scores[,1]
      PLS2_P[,dims] <- PLS2_Res$P[,1]
      PLS2_W[,dims] <- PLS2_Res$W[,1]
      PLS2_Q[,dims] <- PLS2_Res$Q[,1]
      PLS2_U[,dims] <- PLS2_Res$U[,1]

      sv <- sqrt(t(PLS2_Scores[,dims]) %*% PLS2_Scores[,dims]) # For R2X
      varexp.y[dims] <- (t(PLS2_U[,dims]) %*% PLS2_U[,dims])^2 # For R2Y (see varexp)
      qtmp <- PLS2_Scores[,dims] / as.vector(sv)
      # End of addition for PLS

      # qtmp <- usv$u[,1]

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

    ICs <- ICs - 1

    saliences[,dims] <- lambda
    Q[,dims] <- q


    # SingVal is square of the SVD singular value
    # because here it is based on X, not X.X'
    # Does NOT work for Kernel & Wide & Tall
    if(CompMethod == 'Normal' || CompMethod == 'PCT'){
      res_calib$SingVal[dims] <- sv^2 # Calculated from the scores (new for PLS)
      varexp[dims] <- res_calib$SingVal[dims]^2
    } else {
      varexp[dims] <- 0 #This will be overwritten at the end
      res_calib$SingVal[dims] <- 0 #This will be overwritten at the end
    }


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
    meanW <- c(meanW,colMeans(s_r[[j]]))
  }
  #res_calib$s_n <- s_n
  #res_calib$s_r <- s_r
  rm(s_n)
  rm(s_r)

  ##
  #explained <- list()
  explained <- 100*varexp/sum(varexp)
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
  T_Loc <- list()
  for(i in 1:ntable){ # Prepare lists
    #L_CD[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = ncol(temp_tabCalib[[i]]))
    #L_X[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = ncol(temp_tabCalib[[i]]))
    T_Loc[[TableLabels[i]]] <- matrix(,ncol = ndim, nrow = nrowMB)
  }
  b <- matrix(,ncol = ndim, nrow = ntable) # Used to calculate R2X in Kernel and Tall

  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop('The package pracma is needed.')
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

  # If Output==[], nothing else is calculated
  #if(!(is.null(output))){
    # Calculate Global Loadings
    #if(!(is.null(isP))){
      L_CD_Vec <- t(Xnorm_mat)%*%res_calib$Q # Scaled CD 'global' Loadings
    #}

    #if(!(is.null(isL))){
    #  L_X_Vec <- t(X_mat)%*%res_calib$Q # Unscaled X 'global' Loadings
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
  ComDimPLS2_B <- PLS2_W %*% pracma::pinv(t(PLS2_P) %*% PLS2_W) %*% t(PLS2_Q)
  PLS2_Res$ComDimPLS2_B <- ComDimPLS2_B
  ComDimPLS2_B0 <- meanY - meanW %*% ComDimPLS2_B
  PLS2_Res$ComDimPLS2_B0 <- ComDimPLS2_B0
  res_calib$PLS2_Res <- PLS2_Res
  # End of addition for PLS

  res_calib$b <- b
  colnames(res_calib$b) <- DimLabels
  rownames(res_calib$b) <- TableLabels
  #res_calib$b$i <- TableLabels # tables
  #res_calib$b$v <- DimLabels  # Dimensions

  #if(!(is.null(isT))){
    #res_calib$T_Loc$i <- res_calib$Q$i # samples
    #res_calib$T_Loc$v <- DimLabels # Dimensions
    res_calib$T_Loc <- T_Loc

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
    #for(i in 1:ntable){
    #  rownames(res_calib$P_Loc[[TableLabels[i]]]) <- MB@Variables[[i]]
    #  colnames(res_calib$P_Loc[[TableLabels[i]]]) <- DimLabels
    #}
  #}

  rm(L_CD_Vec)
  #rm(L_CD)

  # if(!(is.null(isL))){
  #   res_calib$Lx$i <- t(c(1:nC)) # numbers of all variables;
  #   res_calib$Lx$v <- DimLabels # Dimensions
  #   res_calib$Lx$Data <- L_X_Vec
  #   colnames(res_calib$Lx$Data) <- DimLabels
  #   rownames(res_calib$Lx$Data) <- t(c(1:nC))
  #   #res_calib$Lx_Loc$i <- TableLabels # Tables
  #   res_calib$Lx_Loc$v <- DimLabels # Dimensions
  #   res_calib$Lx_Loc$Data <- L_X
  #   for (i in 1:ntable){
  #     colnames(res_calib$Lx_Loc$Data[[TableLabels[i]]]) <- DimLabels
  #     rownames(res_calib$Lx_Loc$Data[[TableLabels[i]]]) <- data[[i]]$Variables
  #   }
  # }
  # rm(L_X_Vec)
  # rm(L_X)

  ##
  # Singular value calculated from b and saliences
  # Since :
  # b_T_Q(:,j)=saliences.d(:,j).*saliences.d(:,j)/SingVal.d(1,j);
  if(CompMethod == 'Kernel' || CompMethod == 'Tall' || CompMethod == 'Wide'){
    for(dims in 1:ndim){
      res_calib$SingVal[dims] <- t(b[,dims])%*%((saliences[,dims]^2)/as.vector(t(b[,dims])%*%b[,dims]))
    }
    names(res_calib$SingVal) <- DimLabels
    varexp <- res_calib$SingVal^2
    res_calib$explained <- varexp/sum(varexp)*100
    names(res_calib$explained) <- DimLabels
  }

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

  rm(Q)

  ## 07/08/2022 Calculate Y predictive and rescale.
  Ypred <- Xnorm_mat %*% PLS2_Res$ComDimPLS2_B

  for(i in 1:ncol(y)){
    Ypred[,i] <- unlist(Ypred[,i]) + PLS2_Res$ComDimPLS2_B0[i]
  }

  ## 07/08/2022 Predict classes
  if(method == 'PLS-DA'){

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
  } else if(method == 'PLS-R'){
    PRESS <- t(Ypred-y) %*% (Ypred-y)
    Ym <- y - mean(y)
    TSS <- t(Ym) %*% Ym
    Q2 <- 1 - PRESS/TSS
  }

  #Ypred <-Ypred*meanY # Re-scale Ypred. I think this is not needed.

  # The way R2X and R2Y are calculated must be verified.
  #TSS <- sum(diag(Xnorm_mat))
  #R2X <- matrix(rep(NA, ndim),
  #              ncol = 1,
  #              nrow = ndim)
  # R2Y <- vector()
  #
  # #for (i in 1:ndim){
  #   #rss <- sum(diag(Xnorm_mat - (T_mat[,i] %*% t(T_mat[,i]))))
  #   #R2X[i,1] <- 1 - rss/TSS
  #   for(j in 1:ncol(y)){
  #     PRESS <- sum((y[,j] - Ypred[,j])^2)
  #     TSSY <- sum((y[,j] - mean(y[,j]))^2)
  #     R2Y[j] <- 1 - PRESS/TSSY
  #   }
  #
  # #}
  # #rownames(R2X) <- paste0("CC",1:ndim)
  # names(R2Y) <- paste0("Y",1:ncol(y))


  # ssqX = sum(sum(Xnorm_mat*Xnorm_mat))
  # ssqY = sum(sum(y*y))
  #
  #
  # R2X <- R2Y <- vector()
  # for (i in 1:ndim){
  #
  #   for(j in 1:ncol(y)){
  #      # PRESS for Q2
  #     R2X[i] = ((t(T_mat[,i]) %*% T_mat[,i] %*% t(PLS2_P[,i]) %*% PLS2_P[,i])/ssqX)*100
  #     R2Y[i] = ((t(T_mat[,i]) %*% T_mat[,i] %*% t(PLS2_Q[i,1])%*% PLS2_Q[i,1])/ssqY)*100
  #
  #     if(ncol(y) == 2){
  #
  #     }
  #   }
  #
  # }
  #
  # PLS2_Res$R2X <- R2X
  # PLS2_Res$R2Y <- R2Y
  rm(T_mat)

  end_output <- Sys.time()
  running_time <- (end_output - start_ini_time) # Total time of analysis
  if(loquace){
    print(sprintf("Analysis finished after : %s millisecs", running_time*1000))
  }

  progress_bar = utils::txtProgressBar(min=0, max=80, style = 3, char = "=")
  utils::setTxtProgressBar(progress_bar, value = 80)

  close(progress_bar)

  # Save data in a ComDim structure

  res_calib <- new("ComDim",
                   Method = method, # PLS-DA or PLS-R
                   ndim = ndim,
                   Q.scores = res_calib$Q,
                   T.scores = res_calib$T_Loc,
                   P.loadings = res_calib$P,
                   Saliences = res_calib$saliences,
                   Orthogonal = list(),
                   R2X = res_calib$explained,
                   R2Y = r2y,
                   Q2 = Q2,
                   DQ2 = if(method == 'PLS-DA'){
                              DQ2
                          } else {
                              vector()
                          },
                   Singular = res_calib$SingVal,
                   Mean = list(MeanMB = res_calib$MeanMB,
                               MeanY = meanY),
                   Norm = list(NormMB = res_calib$NormMB),
                   PLS.model = list(W = PLS2_W,
                                    B = PLS2_Res$ComDimPLS2_B,
                                    B0 = PLS2_Res$ComDimPLS2_B0,
                                    Y = y), # W, U, B, B0, Y
                   cv = list(),
                   Prediction = if(method == 'PLS-DA'){ # To do?
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
                   runtime = as.numeric(running_time))

  return(res_calib)
}


Compress_Data_2020 <- function(s_n = s_n, CompMethod = CompMethod, Partitions = Partitions){

#' Compress large multi-block objects.
#'
#' Internal function of ComDim_PCA().
#' @param s_n The multi-block object.
#' @param Partitions The number of partitions.
#' @param CompMethod It can be 'Normal' (default), 'Kernel', 'PCT', 'Tall' or 'Wide'.
#' @param Partitions The number of partitions.
#' @examples
#' b1 <- matrix(rnorm(5000),10,500)
#' b2 <- matrix(rnorm(2000),10,200)
#' blist <- list(b1 = b1, b2 = b2)
#' blist <- Compress_Data_2020(blist, 'Normal', 1)
#' @return The compressed multi-block.
#' @export

  s_r <- list()
  ntable <- length(s_n)
  nR <- nrow(s_n[[1]])

  if(!is.null(CompMethod)){
    if (CompMethod == 'Tall'){
      for(i in 1:ntable) {
        MaxRank <- min(c(nrow(s_n[[i]]),ncol(s_n[[i]])))
        # Tall segmented PCT
        s_r[[i]] <- PCA_Tall_PCT_DNR(Xn = s_n[[i]], PCs = MaxRank, Partitions = Partitions)
      }
    } else if (CompMethod == 'Wide'){
      for(i in 1:ntable){
        Xi <- ColumnsPartition(Xn = s_n[[i]], Partitions = Partitions)
        T_all <-   matrix(, nrow = nrow(s_n[[i]]), ncol = 0)
        for(p in 1:Partitions){
          usv <- svd(Xi[[p]])
          T_in <- usv$u%*%diag(usv$d)
          T_all <- cbind(T_all, T_in)
        }
        s_r[[i]] <- T_all

      }
    } else if (CompMethod == 'Kernel'){
      for(i in 1:ntable){
        Kern <- (s_n[[i]]%*%t(s_n[[i]]))/nR
        usv <- svd(Kern)
        s_r[[i]] <- usv$u%*%diag(usv$d)
      }
    } else if (CompMethod == 'PCT'){
      for(i in 1:ntable){
        usv <- svd(s_n[[i]])
        s_r[[i]] <- usv$u%*%diag(usv$d)
      }
    } else if (CompMethod == 'Normal'){
      for(i in 1:ntable){
        s_r[[i]] <- s_n[[i]]
      }
    }
  }
  return(s_r)
}


PCA_Tall_PCT_DNR <-function(Xn = Xn, PCs = PCs, Partitions = Partitions){

#' Compress the multi-block using the PCT method.
#'
#' Internal function of ComDim_PCA().
#' @param Xn The block to compress.
#' @param PCs Number of PCs to keep.
#' @param Partitions The number of partitions.
#' @examples
#' b1 <- matrix(rnorm(5000),10,500)
#' b1 <- PCA_Tall_PCT_DNR(b1,10,2)
#' @return The compressed multi-block.
#' @export

  W <- RowsPartition(Xn, Partitions)
  cols <- ncol(Xn)

  D <- t(W[[1]])%*%W[[1]]

  if(Partitions > 1) {
    for(par in 2:Partitions){
      D <- D + (t(W[[par]])%*%W[[par]])
    }
  }

  vs <- eigen(D)
  vs_flipped <- vs$vectors[seq(length(vs$values),1,-1),]

  Tm <- matrix(, nrow=nrow(Xn), ncol=0)
  Taux <- matrix(, nrow=0, ncol = PCs)

  for(par in 1:Partitions){
    to <- W[[par]]%*%vs_flipped

    Taux <- rbind(Taux, to[,1:PCs])
  }

  Tm <- cbind(Tm, Taux)
  return(Tm)
}

ColumnsPartition <- function(Xn = Xn, Partitions = Partitions){

#' Calculate vertical partitions in the multi-block.
#'
#' Internal function of ComDim_PCA().
#' @param Xn A block.
#' @param Partitions The number of partitions.
#' @examples
#' b1 <- matrix(rnorm(500),10,50)
#' b1 <- ColumnsPartition(b1,2)
#' @return The partitioned block.
#' @export

  cols <- ncol(Xn)
  stride <- floor(cols/Partitions)
  remaining <- cols %% Partitions

  count <- 1
  W <-list()

  for(i in 1:Partitions){
    step <- count + stride - 1
    if(i == Partitions){
      step <- step + remaining
    }
    W[[i]] <- Xn[,count:step]
    count <- count + stride
  }
  return(W)
}

RowsPartition <- function(Xn = Xn, Partitions = Partitions){

#' Calculate horizontal partitions in the multi-block.
#'
#' Internal function of ComDim_PCA().
#' @param Xn A block.
#' @param Partitions The number of partitions.
#' @examples
#' b1 <- matrix(rnorm(500),10,50)
#' b1 <- RowsPartition(b1,2)
#' @return The partitioned block.
#' @export

  rows <- nrow(Xn)
  stride <- floor(rows/Partitions)
  remaining <- rows %% Partitions

  count <- 1
  W <- list()

  for(i in 1:Partitions){
    step <- count + stride - 1
    if(i == Partitions){
      step <- step + remaining
    }
    W[[i]] <- Xn[count:step,]
    count <- count + stride
  }
  return(W)
}

PLS2_DNR_2018 <- function(X = X, Y = Y, lvs = lvs){

  # Store original X
  X_old <- X
  px <- ncol(X)
  rows <- nrow(X)
  # if(is.vector(Y)){
  #   u_old <- as.matrix(Y)
  #   py <- 1
  # } else {
    u_old <- Y[,1] # Any column of Y can be used to initialize u.
    py <- ncol(Y)
  # }

  # Initialize matrices
  W <- matrix(data = rep(0,px*lvs), nrow = px, ncol = lvs)
  P <- W
  Q <- matrix(data = rep(0,py*lvs), nrow = py, ncol = lvs)
  Te <- matrix(data = rep(0,rows*lvs), nrow = rows, ncol = lvs)
  U <- Te

  if(is.matrix(X)){
    meanX <-colMeans(X)
  } else {
    meanX <- mean(X)
  }

  if(is.matrix(Y)){
    meanY <-colMeans(Y)
  } else {
    meanY <- mean(Y)
  }

  # Necessary for PLS_ICA where nLVs decreases to 1
  B0_mat <- B_mat <- list(NULL) # Each item in the list is a matrix with px nrow and py ncol.

  for(i in 1:lvs){
    B_mat[[i]] <- matrix(rep(0,px*py), nrow = px, ncol = py)
    B0_mat[[i]] <- matrix(rep(0,lvs*py), nrow = lvs, ncol = py)
  }

  for (i in 1:lvs){
    limite <- 0
    while (TRUE){ # Repeat until breaking the loop
      w <- t(u_old) %*% X
      w <- w / norm(w, type = '2')
      ti <- X %*% t(w)
      ti <- ti / as.vector(w %*% t(w))
      q <- t(ti) %*% Y
      q <- as.vector(q) / as.vector(t(ti) %*% ti)
      u <- (Y %*% q) / as.vector(t(q) %*% q)
      if ( norm(u - u_old, type = '2') < 1e-6) {
        break
      }
      limite <- limite+1
      if ( limite > 1000) {
        break
      }
      u_old <- u
    }
    p <- t(ti) %*% X
    p <- p / as.vector(t(ti) %*% ti)

    W[,i] <- w
    P[,i] <- p
    Q[,i] <- q
    Te[,i] <- ti
    U[,i] <- u
    X <- X - ti %*% p
    Y <- Y - ti %*% q

    i_nums <-c(1:i)
    #print(dim(pracma::pinv(t(P[,i_nums]) %*% W[,i_nums])))
    if(ncol(Y) == 1){
      B_mat[[i]] <- (W[,i_nums] %*% pracma::pinv(t(P[,i_nums]) %*% W[,i_nums])) %*% as.matrix(Q[,i_nums])
    } else {
      B_mat[[i]] <- (W[,i_nums] %*% pracma::pinv(t(P[,i_nums]) %*% W[,i_nums])) %*% t(Q[,i_nums])
    }

    B0_mat[[i]] <- meanY - meanX %*% B_mat[[i]]
  }

  B <- (W  %*% pracma::pinv(t(P) %*% W)) %*% t(Q)
  B0 <- meanY - meanX %*% B

  PLS2_Res <- list()
  PLS2_Res$B <- B
  PLS2_Res$B0 <- B0
  PLS2_Res$Scores <- Te
  PLS2_Res$P <- P
  PLS2_Res$W <- W
  PLS2_Res$U <- U
  PLS2_Res$Q <- Q

  PLS2_Res$B_mat <- B_mat
  PLS2_Res$B0_mat <- B0_mat
  PLS2_Res$Yhat <- X_old %*% B + rep(1,rows) %*% B0

  return(PLS2_Res)
}
