# It is based to ComDim OPLS from Julien Boccard.

# It requires pracma::inv, which is invoked in the package.
# Need to change all matrix divisions by pracma:inv * ...
# Need to check that the result is the same than in Matlab.
# Check how the saliences are computed.
# See if the B values make sense (on June DNR was skeptical on how to deal with this situation in the multi-block scenareo due to the different scaling...)

library(utils)
# library(ggplot2)
# library(ggpubr)
require(pracma)

ComDim_KOPLS_MB <- function(MB = MB,  # collection = The kernel matrix (un-centered) from ComDim ; see 'koplsKernel()' for details.
                        y = y,  # Y = The response matrix (un-centered/scaled). Could be binary (for discriminant analysis) or real-valued.
                        max.ort =  2, # northolvs = Max number of Y-orthogonal components (integer).
                        normalise = TRUE,
                        nrcv = NULL, # nrcv = Number of cross-validation rounds (integer).
                        # cvType = # cvType = Type of cross-validation. See also 'koplsCrossValSet()' for details.
                        #         c('nfold', # Either 'nfold' for n-fold
                        #            'mccv', # Monte Carlo CV
                        #            'mccvb')[1],  # Monte Carlo class-balanced CV.
                        decisionRule = c('fixed','max')[2],
                        preProcK = c('mc','no')[2], # preProcK = Pre-processing settings for the kernel matrix.
                          # Either 'mc' for mean-centering or 'no' for no pre-processing.
                        cvFrac = 0.7, # cvFrac = Fraction of observations in the training set during cross-validation. Only applicable for 'mccv' or 'mccvb'.
                        method = c('k-OPLS-DA','k-OPLS-R')[1], # discriminant analysis, regression.
                        #   If 'da', sensitivity and specificity will be calculated.
                        loquace = FALSE  # verbose. If zero, no output will be displayed, otherwise some
                                      # output will be displayed regarding the cross-validation progress (default).
                        ) {

    # OUTPUT
    # res_calib = Object with 'A' predictive components
    # and 'northolvs' Y-orthogonal components.
    #
    # Contains the following entries from 'ComDimOPLSCV' :
    #  cv = Cross-validation results:
    #  Q2Yhat = Total Q-square result for all Y-orthogonal components.
    #	 Q2YhatVars = Q-square result per Y-variable for all Y-orthogonal
    #       components.
    #	 Yhat = All predicted Y values as a concatenated matrix.
    #	 Tcv = Predictive score vector T for all cross-validation rounds.
    #  cvTrainIndex = Indices for the training set observations during
    #       the cross-validation rounds.
    #	 cvTestIndex = Indices for the test set observations during the
    #       cross-validation rounds.
    #  da = Cross-validation results specifically for discriminant
    #       analysis (DA) cases:
    #    predClass = Predicted class list per class and Y-orthogonal
    #       components (integer values).
    #	 trueClass = True class list per class and Y-orthogonal
    #       components (integer values).
    #	 sensSpec = Sensitivity and specificity values per class and
    #       Y-orthogonal components (integer values).
    #	 confusionMatrix = Confusion matrix during cross-validation
    #       rounds.
    #	 nclasses = Number of classes in model.
    #	 decisionRule = Decision rule used: 'max' or 'fixed'.
    #
    # Plus :
    # res_calib.cv.OrthoLVsOptimalNum = optimal number of Y-orthogonal components
    # res_calib.cv.Yhat = All predicted Y values as a concatenated matrix.
    # res_calib.ComDimOPLSModel.loadings = Loadings
    #
    # if 'da'
    # res_calib.cv.DQ2Yhat = Total Q-square result for Y-orthog components
    #
    # Calls :
    # ComDimOPLSCV
    # DQ2

  ndim = 1
  cvType = 'nfold'

  start_ini_time <- Sys.time()  # To start counting calculation time for the initialization

  if(!(method == 'k-OPLS-DA' || method == 'k-OPLS-R')){
    stop("'method' must be 'k-OPLS-DA' or 'k-OPLS-R'.")
  }

  # INITIALISATION
  if (loquace) {
    progress_bar <- utils::txtProgressBar(min = 0, max = 100, style = 3)
    pieceBar <- 100/(nrcv+1)
    total_progress <- pieceBar
    cat('\n')
    cat('Start of analysis\n')
    utils::setTxtProgressBar(progress_bar, value = 0)
    cat('\n')
  }

  if(class(MB) != 'MultiBlock'){
    stop("'MB' is  not a MultiBlock.")
  }

  ntable = length(getBlockNames(MB)) # number of blocks
  nrowMB = length(getSampleNames(MB)) # number of samples

  if(is.null(nrcv)){
    if(nrowMB >= 10){
      nrcv <- 5
    } else {
      nrcv <- round(nrowMB/5)
    }
  }

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
    stop('All blocks had 0 variance')
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

  if(is.null(ndim)){
    ndim <- ntable # If the number of components is not defined,
    # the number of components to extract is equal to the number of blocks.
  }

  DimLabels <- paste0('CC',1:ndim)   # One label per component.
  TableLabels<- getBlockNames(MB) # One label per block.


  end_ini_time <- Sys.time() # To end the count of the calculation time.

  if (loquace) {
    print(sprintf("Initialisation finished after : %s millisecs", (end_ini_time - start_ini_time)*1000))
    utils::setTxtProgressBar(progress_bar, value = pieceBar)
  }


  if(!(is.vector(y))){
    N <- nrow(y)
    m <- ncol(y)
  } else {
    N <- length(y)
    m <- 1
  }

  ## Check y matrix (convert to dummy matrix if it is not already.)
  if(!is.matrix(y)){ # If it's a vector and the method is OPLS-R
    y <- as.matrix(y)
  }

  if(method == 'k-OPLS-DA'){

    tmp <- unique(as.vector(y))
    if(all(tmp %in% c(0,1))){ # Is it dummy?

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
    nclasses <- length(unique(classVect))
  }
  MeanY <- colMeans(y)

  # if(method == 'k-OPLS-DA'){
  #
  #   decisionRule <- 'max' #move to arg... %this is a parameter for DA decision rule
  #   tmp <- unique(y)
  #   print(tmp)
  #   if(all(tmp == c(0,1))){
  #     if(m == 1){
  #       #---koplsDummy <- function(class, numClasses = length(unique(class))){
  #       labels <- unique(y) #the set of classlabels in class
  #       numClasses = length(unique(y))
  #       labels_sorted <- sort(labels) #sort labels in ascending order
  #       len_class <- length(y)#number of samples
  #       dummy <- matrix(rep(0,len_class*numClasses), nrow = len_class, ncol = numClasses)
  #       #dummy matrix initialized as a zero matrix
  #       for (i in 1:numClasses){ #for each class label
  #         ind2 <- which(y==labels_sorted[i])
  #         #find the rows (samples) that belongs to the current class, labels_sorted(i
  #         dummy[ind2,i] <- 1 #!!!write ones in the positions where we have the current class....
  #       }
  #       rm(labels, numClasses, labels_sorted, len_class, ind2)
  #     }
  #     #classVect <- koplsReDummy(y)
  #     if(!is.matrix(y)){
  #       y <- as.matrix(y)
  #     }
  #     classVect <- rep(NaN, nrow(y))
  #     for(i in 1:ncol(y)){
  #       ind2 <- which(y[,i]==1)
  #       classVect[ind2] <- i
  #     }
  #     rm(ind2)
  #
  #   } else if(all(y %% 1 == 0) && m == 1){
  #     classVect <- y
  #     labels <- unique(y+1) #the set of classlabels in class
  #     numClasses = length(unique(y))
  #     labels_sorted <- sort(labels) #sort labels in ascending order
  #     len_class <- length(y)#number of samples
  #     dummy <- matrix(rep(0,len_class*numClasses), nrow = len_class, ncol = numClasses)
  #     #dummy matrix initialized as a zero matrix
  #     for (i in 1:numClasses){ #for each class label
  #       ind2 <- which(y==labels_sorted[i])
  #       #find the rows (samples) that belongs to the current class, labels_sorted(i
  #       dummy[ind2,i] <- 1 #!!!write ones in the positions where we have the current class....
  #     }
  #     rm(labels, numClasses, labels_sorted, len_class, ind2)
  #   } else {
  #     stop("'method' is k-OPLS-DA, but y appears to be neither dummy (1 0) matrix nor a vector of (integer) class labels")
  #   }
  #   nclasses <- length(unique(classVect))
  # }

  if(cvType == 'mccvb' && method != 'k-OPLS-DA'){
    stop("Class balanced monte-carlo cross validation only applicable to 'k-OPLS-DA' method.")
  }

  if(!(any(c(cvType == 'mccvb', cvType == 'mccv', cvType == 'nfold')))){
    stop(sprintf('%s - unknown Cross-validation type', cvType))
  }

  # convert Y-scaling to more explicit format
  #YcenterType <- 'no'
  #YscaleType <- 'no'
  #if (preProcY != 'no'){
  #  YcenterType <- 'mc'
  #  if (preProcY != 'mc'){
  #    YscaleType <- preProcY
  #  }
  #}

  Yhat <- rep(NaN, m)
  YhatDaSave <- list()
  pressyVars <- list()
  pressyVarsTot <- list()
  cvTestIndex  <- vector()
  cvTrainIndex <- vector()

  # Normalisation step
  res_calib <- list()
  res_calib$SingVal$Data <-as.vector(NULL)
  res_calib$NormMB$Data <-as.vector(NULL)
  res_calib$MeanMB$Data <-list()

  for (i in 1:ntable){

    res_calib$MeanMB$Data[[TableLabels[i]]] <- colMeans(MB@Data[[i]])
    res_calib$MeanMB$Variables[[TableLabels[i]]] <- MB@Variables[[i]]

    if(normalise){
      MB@Data[[i]]<- MB@Data[[i]] - t(t(rep(1,nrowMB))) %*% t(as.matrix(res_calib$MeanMB$Data[[TableLabels[i]]]))
      res_calib$NormMB$Data[i] <- norm(MB@Data[[i]],'2')
      MB@Data[[i]] <- MB@Data[[i]]/norm(MB@Data[[i]],'2')
    }
    temp <- MB@Data[[i]] %*% t(MB@Data[[i]])

    # Generate the Kernel matrix
    if(i == 1){
      K <- (1/ntable)*temp
    } else {
      K <- K+(1/ntable)*temp
    }
  }

  AllYhat <- vector()

  # Cross-validation
  pressyTot <- pressy <- matrix(, nrow = max.ort+ndim, ncol = 1)
  for(icv in 1:nrcv){

    if (loquace){
      cat(sprintf('Please wait... cv round: %d of %d\n', icv, nrcv))
      total_progress <- total_progress + pieceBar
      utils::setTxtProgressBar(progress_bar, value = total_progress)
      cat('\n')
    }

    #set up CV -----------------

    # ComDimOPLSCrossValSet function was absorbed by ComDim_OPLS
    cvSet <- list()

    if(is.vector(y)){
      y <- as.matrix(y)
    }

    modInd <- vector()
    predInd <- vector()

    # if(cvType == 'mccvb'){ #'Monte-carlo cv - class Balanced'
    #   #check if Y is dummy or labels...
    #   tmp <- unique(c(y))
    #   if(all(tmp %in% c(0,1))){ # If y is a dummy vector
    #     # classVect <- koplsReDummy(y)
    #     ny <- nrow(y)
    #     my <- ncol(y)
    #     classVect <- rep(NA,ny)
    #     for (i in 1:my){
    #       ind <- which(y[,i] == 1)
    #       classVect[ind] <- i
    #     }
    #   } else {
    #     classVect <- y
    #   }
    #
    #   minset <- unique(classVect) #find all classlabels
    #   for (i in 1:length(minset)){
    #     #for each class
    #     currentClass <- minset[i] #current class label
    #     ind <- which(classVect==currentClass) #find all samples from current class
    #
    #     #randomize
    #     ran <- rnorm(length(ind)) #randomize
    #     sorted <- sort.int(ran, index.return = TRUE) #sort randomize number to get randomized index string
    #     ind <- ind[sorted$ix] #apply randomization on the real index vector
    #     #-end randomize
    #
    #     modelLim <- ceiling(length(ind)*cvFrac)
    #     #number of elements that will get assigned as model
    #     modInd <- c(modInd,ind[1:modelLim])
    #     predInd <- c(predInd,ind[(modelLim+1):length(ind)])
    #   }
    # }
    #
    # if(cvType == 'mccv'){ #'Monte-carlo cv'
    #
    # #randomize
    #   ran <- rnorm(length(K[,1]),1) #randomize
    #   sorted <- sort.int(ran, index.return = TRUE) #sort randomize number to get randomized index string
    #   ind <- seq(1,length(ran),1)
    #   ind <- ind[sorted$ix] #apply randomization on the real index vector
    #   modelLim <- ceiling(length(ind)*cvFrac) #number of elements that will get assigned as model
    #   modInd <- ind[1:modelLim]
    #   predInd <- ind[(modelLim+1):length(ind)]
    # }

    if(cvType == 'nfold') { #'N-Fold cross validation'
      predInd <- seq(icv, length(y[,1]), by = nrcv)
      modInd <- setdiff(1:length(y[,1]),predInd)
    }

    cvSet$cvType <- cvType
    cvSet$nfold <- nrcv
    cvSet$nfoldRound <- icv

    cvSet$KTrTrcol <- list()
    KTrTrcol <- list()
    for (ta in 1:ntable){
      KTrTrcol$Data <- MB@Data[[ta]][modInd,]
      cvSet$KTrTrcol[[ta]] <- KTrTrcol
    }

    cvSet$KTrTr <- K[modInd,modInd]
    cvSet$KTeTr <- K[predInd,modInd]
    cvSet$KTeTe <- K[predInd,predInd]

    cvSet$yTraining <- y[modInd,]
    cvSet$yTest <- y[predInd,]

    cvSet$trainingIndex <- modInd
    cvSet$testIndex <- predInd
    # --- end ComDimOPLSCrossValSet
    cvTestIndex <- append(cvTestIndex,cvSet$testIndex)
    cvTrainIndex <- append(cvTrainIndex,cvSet$trainingIndex)

    #get Kernel matrices ------- change so that this is done in the K
    #matrix only once and selected by indeces.
    KtrTrcol <- cvSet$KTrTrcol
    KtrTr <- cvSet$KTrTr
    KteTr <- cvSet$KTeTr
    KteTe <- cvSet$KTeTe

    #--- koplsScale <- function(X,YcenterType, YscaleType){
    if(!(is.matrix(cvSet$yTraining))){
      cvSet$yTraining <- as.matrix(cvSet$yTraining)
    }

    YScaleObj <- list()
    YScaleObj$meanV <- colMeans(cvSet$yTraining)
    YScaleObj$stdV <- apply(cvSet$yTraining, 2, sd)
    YScaleObj$X <- cvSet$yTraining

    YScaleObj$preProcY <- 'no' # No pre-processing for the y-block
    #---
    #YScaleObjTest <- koplsScaleApply(cvSet$yTest,YScaleObj)
    YScaleObjTest <- list()

    #YScaleObjTest$centerType <- YScaleObj$centerType
    #YScaleObjTest$scaleType <- YScaleObj$scaleType
    YScaleObjTest$meanV <- YScaleObj$meanV
    YScaleObjTest$stdV <- YScaleObj$stdV
    ncolX <- ncol(cvSet$yTest)
    nrowX <- nrow(cvSet$yTest)

    YScaleObjTest$X <- cvSet$yTest
    #---

    if (preProcK == 'mc'){
      #---KteTe <- koplsCenterKTeTe(KteTe,KteTr,KtrTr)
      nrowKer <- nrow(KteTr)
      ncolKer <- ncol(KteTr)
      Itrain <- diag(nrowKer)
      I_nTrain <- as.matrix(rep(1,ncolKer), nrow = ncolKer, ncol = 1)
      nTrain <- nrowKer
      I <- diag(ncolKer)
      I_n <- as.matrix(rep(1,nrowKer), nrow = nrowKer, ncol = 1)
      n <- ncolKer
      D_te <-  (1/nTrain) * I_n %*% t(I_nTrain)
      KteTe <- KteTe - D_te %*% t(KteTr) - KteTr %*% t(D_te) + D_te %*% KtrTr %*% t(D_te)
      rm(nrowKer, ncolKer, Itrain, I_nTrain, I, I_n, n, D_te)
      #---
      #---KteTr <- koplsCenterKTeTr(KteTr,KtrTr)
      nrowKrr <- nrow(KtrTr)
      ncolKrr <- ncol(KtrTr)
      nrowKer <- nrow(KteTr)
      ncolKer <- ncol(KteTr)
      Itrain <- diag(nrowKrr)
      I_nTrain <- as.matrix(rep(1,nrowKrr), nrow = nrowKrr, ncol = 1)
      nTrain <- nrowKrr
      I <- diag(nrowKer)
      I_n <- as.matrix(rep(1,nrowKer), nrow = nrowKer, ncol = 1)
      n <- nrowKer
      KteTr <- (KteTr-(1/nTrain) * I_n %*% t(I_nTrain) %*% KtrTr) %*% (Itrain-(1/nTrain)*I_nTrain%*%t(I_nTrain))
      rm(nrowKrr, ncolKrr, nrowKer, ncolKer, I_nTrain, nTrain, I, I_n, n)
      #---
      #---KtrTr <- koplsCenterKTrTr(KtrTr)
      nrowK <- nrow(KtrTr)
      I <- diag(nrowK)
      I_n <- matrix(rep(1,nrowK), ncol = 1, nrow = nrowK)
      KtrTr <- (I- (1/nrowK)* I_n%*%t(I_n)) %*% KtrTr %*%(I-(1/nrowK)*I_n%*%t(I_n))
      rm(I, I_n)
      #---
    }

    #estimate K-OPLS model------- Changed for ComDim
    model <- ComDimOPLSModel(KtrTrcol,YScaleObj,ndim,max.ort)

    #set up model stats----------
    if(is.vector(YScaleObjTest$X)){ # For PLS-R
      ssy <- sum(sum((YScaleObjTest$X)^2))
      ssyVars <- sum((YScaleObjTest$X)^2)
      ssx <- sum(diag(KteTe))
    } else {# For PLS-DA
      ssy <- sum(colSums((YScaleObjTest$X)^2))
      ssyVars <- colSums((YScaleObjTest$X)^2)
      ssx <- sum(diag(KteTe))
    }


    if(icv==1) {
      ssyTot <- ssy
      ssyVarsTot <- ssyVars
      ssxTot <- ssx
    } else{
      ssyTot <- ssyTot+ssy
      ssyVarsTot <- ssyVarsTot+ssyVars
      ssxTot <- ssxTot+ssx
    }

    #for each combination of Y-osc components
    for( ioax in 1:(max.ort+1)){

      for( ioay in 1:1){ #oay+1)

        # ComDimOPLS predict yYhat
        modelPredy <- ComDimOPLSPredict(KteTr,KteTe,KtrTr, model,ioax-1,0)

        #---tmp <- koplsRescale(YScaleObj,modelPredy$Yhat)
        tmp <- list()
        tmp$X <- modelPredy$Yhat
        tmp$preProcY <- 'no'
        #---

        if(ioax == 1){
          AllYhatind <- tmp$X
        } else {
          AllYhatind <- cbind(AllYhatind, tmp$X)
        }

        pressy[ioax,ioay] <- sum(colSums((YScaleObjTest$X-modelPredy$Yhat)^2))
        if(length(pressyVars) < ioax){
          pressyVars[[ioax]] <- list()
        }

        pressyVars[[ioax]][[ioay]] <- colSums((YScaleObjTest$X-modelPredy$Yhat)^2)
        if(length(pressyVarsTot) < ioax){
          pressyVarsTot[[ioax]] <- list()
          if(ioay == 1){
            pressyVarsTot[[ioax]][[ioay]] <- list()
          }
        }

        if((icv==1)) {
          pressyTot[ioax,ioay] <- pressy[ioax,ioay]
          pressyVarsTot[[ioax]][[ioay]] <- pressyVars[[ioax]][[ioay]]
        } else {
          pressyTot[ioax,ioay] <- pressyTot[ioax,ioay]+pressy[ioax,ioay]
          pressyVarsTot[[ioax]][[ioay]] <- pressyVarsTot[[ioax]][[ioay]] + pressyVars[[ioax]][[ioay]]
        }


        #if 'da' save Yhat for all rounds
        if(method == 'k-OPLS-DA'){

          if(length(YhatDaSave) < ioax){
            YhatDaSave[[ioax]] <- list()
          }

          if(icv == 1){
            YhatDaSave[[ioax]][[ioay]] <-  vector()
          }

          #+mean on Yhat
          #---tmp <- koplsRescale(YScaleObj,modelPredy$Yhat)
          tmp <- list()
          tmp$X <- modelPredy$Yhat
          tmp$preProcY <- 'no'
          #---
          YhatDaSave[[ioax]][[ioay]] <- rbind(YhatDaSave[[ioax]][[ioay]],tmp$X)
        }

        #if highest number of oscs - save Yhat and Xhat
        if(ioax == max.ort+1){
          if(icv == 1){ # && ioay==oay+1)
            Yhat <- vector()
          }

          #---tmp <- koplsRescale(YScaleObj,modelPredy$Yhat)
          tmp <- list()
          tmp$X <- modelPredy$Yhat
          tmp$preProcY <- 'no'
          #---
          Yhat <- rbind(Yhat,tmp$X)
        }
      }
    }
    if(icv == 1){
      AllYhat <- AllYhatind
    } else {
      AllYhat <- rbind(AllYhat, AllYhatind)
    }
    ## This is the end of the crossvalidation.
  }

  if (loquace) {
    cat('Finishing up...\n')
    total_progress <- total_progress + pieceBar
    utils::setTxtProgressBar(progress_bar, value = total_progress)
    cat('\n')
  }

  res_calib$ComDimOPLSModel <- ComDimOPLSModel(MB,y,ndim,max.ort)

  res_calib$cv$Yhat <- Yhat
  res_calib$cv$AllYhat <- AllYhat
  res_calib$cv$Tcv <- Yhat%*%res_calib$ComDimOPLSModel$Cp[[1]]%*%res_calib$ComDimOPLSModel$Bt[[max.ort+1]]

  res_calib$cv$Q2Yhat <- vector()
  res_calib$cv$Q2YhatVars <- matrix(, nrow = max.ort+1, ncol = length(ssyVarsTot))

  for( ioax in 1:(max.ort+1)){
    #res_calib$cv$Q2YhatVars[[ioax]] <- list()
    #for( ioay in 1:1) {
    #res_calib$cv$Q2Yhat[ioax,ioay] <- 1-pressyTot[ioax,ioay]/ssyTot
    #res_calib$cv$Q2YhatVars[[ioax]][[ioay]] <- 1-pressyVarsTot[[ioax]][[ioay]]/ssyVarsTot
    res_calib$cv$Q2Yhat[ioax] <- 1-pressyTot[ioax,1]/ssyTot
    res_calib$cv$Q2YhatVars[ioax,] <- 1-pressyVarsTot[[ioax]][[1]]/ssyVarsTot
    #}
  }
  colnames(res_calib$cv$Q2YhatVars) <- names(pressyVarsTot[[1]][[1]])

  res_calib$cv$cvTestIndex <- cvTestIndex
  res_calib$cv$cvTrainIndex <- cvTrainIndex

  if(method == 'k-OPLS-DA'){
      # get sens/spec for each y-orth component... eval of model
    # for( i in 1:(max.ort+1)){ #we would have no osc comps for dummy matrix...
    #   if(decisionRule == 'max') {
    #
    #     predClass <- vector()
    #     for(ii in 1:nrow(YhatDaSave[[i]][[1]])){
    #       xx <- which.max(YhatDaSave[[i]][[1]][ii,] == max(YhatDaSave[[i]][[1]][ii,]))
    #       if(length(xx) == 1) {
    #         predClass[ii] <- which(YhatDaSave[[i]][[1]][ii,] == max(YhatDaSave[[i]][[1]][ii,]))
    #       } else { # If there is a tie
    #         predClass[ii] <- NaN
    #       }
    #     }
    #     predClass <- matrix(predClass, nrow = length(predClass))
    #     rm(xx)
    #
    #   } else if(decisionRule == 'fixed') {
    #     predClass <- koplsBasicClassify(YhatDaSave[[i]][[1]],1/nclasses)
    #     predClass <- matrix(rep(NaN, ncol(YhatDaSave[[i]][[1]])*nrow(YhatDaSave[[i]][[1]])),
    #                         nrow = nrow(YhatDaSave[[i]][[1]]),
    #                         ncol = ncol(YhatDaSave[[i]][[1]]))
    #     for(ii in 1:nrow(YhatDaSave[[i]][[1]])){
    #       xx <- which(YhatDaSave[[i]][[1]][ii,] > 1/nclasses)
    #       predClass[ii,1:length(xx)] <- xx
    #     }
    #     rm(xx)
    #   } else {
    #     warning(sprintf('Decision rule given: %s is not valid/implemented', decisionRule))
    #   }
    #
    #   specs <- koplsSensSpec(classVect[cvTestIndex], predClass, y)
    #   da <- list()
    #   da$sensAllOsc[[i]] <- specs$sensvec # This is the sensibility...
    #   da$specAllOsc[[i]] <- specs$specvec
    #   da$classvecAllOsc[[i]] <- specs$classvec
    #   da$tot_sensAllOsc[[i]] <- specs$tot_sens
    #   da$meanSensAllOsc[[i]] <- specs$meanSens
    #   da$meanSpecAllOsc[[i]] <- specs$meanSpec
    # }

    # get sens/spec for max number of oscs.... (hmm redundant).

    if(decisionRule == 'max'){
      predClass <- vector()
      for(i in 1:nrow(Yhat)){
        xx <- which(Yhat[i,]== max(Yhat[i,]))
        if(length(xx) == 1) {
          predClass[i] <- which(Yhat[i,] == max(Yhat[i,]))
        } else {
          predClass[i] <- NaN
        }
      }
      predClass <- as.matrix(predClass)
      rm(xx)
    } else if(decisionRule == 'fixed'){
      predClass <- matrix(rep(NaN, ncol(Yhat)*nrow(Yhat)), nrow = nrow(Yhat), ncol = ncol(Yhat))
      for(ii in 1:nrow(Yhat)){
        xx <- which(Yhat[ii,] > 1/nclasses)
        predClass[ii,1:length(xx)] <- xx
      }
      rm(xx)
    } else {
      warning(sprintf('Decision rule given: %d is not valid/implemnted', decisionRule))
    }

    da <- koplsSensSpec(classVect[cvTestIndex], predClass, y)

    da$trueClass <- classVect[cvTestIndex]
    da$nclasses <- nclasses
    res_calib$da <- da
    res_calib$da$predClass <- colnames(y)[predClass]
    res_calib$da$decisionRule <- decisionRule
    #CHANGE TO ORIGNAL ORDER IF NFOLD CV - for backward
    #compatibility and comparison w/ simca-p etc
    if(cvType == 'nfold'){
      cvOrder <- order(cvTestIndex)
      res_calib$da$predClass <- res_calib$da$predClass[cvOrder]
      res_calib$da$trueClass <- res_calib$da$trueClass[cvOrder]
    }
  }

  #CHANGE TO ORIGNAL ORDER IF NFOLD CV - for backward
  #compatibility and comparison w/ simca-p etc
  if(cvType == 'nfold'){
    cvOrder <- order(cvTestIndex)
    res_calib$cv$Yhat <- res_calib$cv$Yhat[cvOrder,]
    res_calib$cv$Tcv <- res_calib$cv$Tcv[cvOrder,]
  }

  res_calib$args$max.ort <- max.ort
  #res_calib$args$oay <- oay
  res_calib$args$ndim <- ndim
  res_calib$class <- 'koplscv'


  dqq <- matrix(, nrow = max.ort+1, ncol = ncol(y))
  #PRESSD <- vector()
  output <- list()

  if(!(is.data.frame(y))){
    y <- as.data.frame(y)
  }

  if (method == 'k-OPLS-DA'){
    Ylarg <- ncol(y)
    for (i in 0:max.ort) {
      for(j in 1:Ylarg){
        #print(res_calib$cv$AllYhat[,Ylarg*i+1])
        #print(y[,1])
        #print('...')
        output$dq <- DQ2(res_calib$cv$AllYhat[,i*Ylarg+j],y[,j])  # Compute DQ2 index
        dqq[j,i+1] <- output$dq$dqq
        #PRESSD[i+1] <- output$dq$PRESSD[i+1]
        # Search for the optimal model based on DQ2 (minimum size=2)
      }
    }
    index <- ndim
    colnames(dqq) <- colnames(res_calib$cv$Q2YhatVars)

    while(index < max.ort+1 && mean(dqq[index+1,])-mean(dqq[index,]) > 0.01){ # 1 percent DQ2 increase criterion
      index <- index+1
    }

    res_calib$cv$DQ2Yhat <- t(dqq[index,])

    res_calib$cv$OrthoLVsOptimalNum <- index - ndim


  } else if(method == 'k-OPLS-R'){ # Search for the optimal model based on Q2Yhat (minimum size=2)
    Ylarg <- ncol(y)
    index <- ndim

    while(index < max.ort+ndim && (res_calib$cv$Q2Yhat[index+1]-res_calib$cv$Q2Yhat[index]>0.01)){  # 1 percent Q2 increase criterion
      index <- index+ndim
    }

    res_calib$cv$OrthoLVsOptimalNum <- index - ndim
  } else{
    stop("'method' should be either 'k-OPLS-DA' or 'k-OPLS-R'.")
  }

  # Adjust Yhat to the selected model size
  # res_calib.cv.Yhat=res_calib.cv.AllYhat(:,(Ylarg*A)+(res_calib.cv.OrthoLVsOptimalNum*A):(Ylarg*A)+(res_calib.cv.OrthoLVsOptimalNum*A)+Ylarg-1);
  res_calib$cv$Yhat <- res_calib$cv$AllYhat[,(Ylarg*ndim+(res_calib$cv$OrthoLVsOptimalNum*ndim)):((Ylarg*ndim)+(res_calib$cv$OrthoLVsOptimalNum*ndim)+Ylarg-1)]

  # Compute the loadings for the selected model size
  for (ta in 1:ntable){

    for (m in 1:ndim){

      # Uses Tpo, the last vector of Tp
      xx <- (res_calib$ComDimOPLSModel$saliences[ta,m]*
                              t(MB@Data[[ta]]) %*%
                                res_calib$ComDimOPLSModel$Tpo[,m])/
                                  drop(t(res_calib$ComDimOPLSModel$Tpo[,m]) %*%
                                         res_calib$ComDimOPLSModel$Tpo[,m])
      if(m == 1){
        loadings.pred <- xx
      } else{
        loadings.pred <- cbind.data.frame(loadings.pred, xx)
      }
    }
    if(ta == 1){
      loadings.pred2 <- loadings.pred
    } else {
      loadings.pred2 <- rbind.data.frame(loadings.pred2, loadings.pred)
    }
    if (res_calib$cv$OrthoLVsOptimalNum >= 1){
      for (n in 1:res_calib$cv$OrthoLVsOptimalNum){
        # Uses all the To
        xx <- (res_calib$ComDimOPLSModel$saliences[ta,n+m]*
                                    t(MB@Data[[ta]]) %*%
                                    res_calib$ComDimOPLSModel$To[,n])/
                                      drop(t(res_calib$ComDimOPLSModel$To[,n]) %*%
                                        res_calib$ComDimOPLSModel$To[,n])
        if(n == 1){
          loadings.ort <- xx
        } else{
          loadings.ort <- cbind.data.frame(loadings.ort, xx)
        }
      }
      if(ta == 1){
        loadings.ort2 <- loadings.ort
      } else {
        loadings.ort2 <- rbind.data.frame(loadings.ort2, loadings.ort)
      }
    } else {
      loadings.ort2 <- NULL
    }
  }
  loadings.pred <- loadings.pred2
  loadings.ort <- loadings.ort2
  rm(loadings.pred2, loadings.ort2)
  #res_calib$ComDimOPLSModel$loadings <- loadings.pred

  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop('The package pracma is needed.')
  }


  # 07/08/2022 Calculate local scores
  T_Loc_o <- T_Loc <- list()
  for(i in 1:ntable){

    if(res_calib$cv$OrthoLVsOptimalNum >= 1){

      T_Loc_o[[TableLabels[i]]] <- matrix(, nrow = nrowMB, ncol = res_calib$cv$OrthoLVsOptimalNum)
      for(j in 1:res_calib$cv$OrthoLVsOptimalNum){

        tempo <- t(MB@Data[[i]])%*%res_calib$ComDimOPLSModel$To[,j] #Scaled CD 'local' Loadings

        #T_Loc[[TableLabels[i]]][,j] <- temp_tabCalib[[i]]%*%(temp/as.vector(t(temp)%*%temp)) # local Scores
        T_Loc_o[[TableLabels[i]]][,j] <- unlist(MB@Data[[i]]%*%(tempo %*% pracma::pinv(t(tempo)%*%tempo))) # local Scores

        # Deflate each temp_tabCalib (substract orthogonal and predictive information)
        MB@Data[[i]] <- MB@Data[[i]]- res_calib$ComDimOPLSModel$To[,j]%*%t(tempo)
      }
    } else {
      T_Loc_o[[TableLabels[i]]] <- NULL
    }
  }


  for(i in 1:ntable){

    T_Loc[[TableLabels[i]]] <- matrix(, nrow = nrowMB, ncol = ndim)
    for(j in 1:ndim){

      temp <- t(MB@Data[[i]])%*%res_calib$ComDimOPLSModel$Tpo[,j] #Scaled CD 'local' Loadings

      #T_Loc[[TableLabels[i]]][,j] <- temp_tabCalib[[i]]%*%(temp/as.vector(t(temp)%*%temp)) # local Scores
      T_Loc[[TableLabels[i]]][,j] <- unlist(MB@Data[[i]]%*%(temp %*% pracma::pinv(t(temp)%*%temp))) # local Scores

      # Deflate each temp_tabCalib (substract orthogonal and predictive information)
      MB@Data[[i]] <- MB@Data[[i]]- res_calib$ComDimOPLSModel$Tpo[,j]%*%t(temp)
    }
  }
 # End 07/08/2022

  # Define block length
  belong_block <- rep(0, nrow(loadings.pred))
  k <- 1
  for(i in 1:ntable){
    belong_block[k:(k+variable_number[i]-1)] <- TableLabels[i]
    k <- k+variable_number[i]
  }

  end_output <- Sys.time()
  running_time <- (end_output - start_ini_time)
  if(loquace){
    print(sprintf("Analysis finished after : %s seconds", running_time))
  }

  if(res_calib$cv$OrthoLVsOptimalNum >= 1){
    Q.scores.ort <- res_calib$ComDimOPLSModel$To[,1:res_calib$cv$OrthoLVsOptimalNum]
    Saliences.ort <- res_calib$ComDimOPLSModel$saliences[,(ndim+1):(ndim+max.ort)]
  } else {
    Q.scores.ort <- NULL
    Saliences.ort <- NULL
  }

  res_calib$ComDimOPLSModel$saliences <- as.matrix(res_calib$ComDimOPLSModel$saliences[,1:ndim])
  res_calib$cv$Tcv <- as.matrix(res_calib$cv$Tcv) # In case ndim is 1, it would be a vector.

  # Save data in a ComDim structure
  res_calib <- new("ComDim",
                   Method = method, # k-OPLS-DA or k-OPLS-R
                   ndim = ndim,
                   Q.scores = res_calib$cv$Tcv, #Predictive score vector T for all cross-validation rounds.
                   T.scores = T_Loc, #res_calib$ComDimOPLSModel$T,
                   P.loadings = as.matrix(loadings.pred),
                   Saliences = res_calib$ComDimOPLSModel$saliences,
                   Orthogonal = list(
                                nort = res_calib$cv$OrthoLVsOptimalNum,
                                Q.scores.ort = Q.scores.ort,
                                T.scores.ort = T_Loc_o,
                                P.loadings.ort = loadings.ort,
                                Saliences.ort = Saliences.ort,
                                #Y.ort = res_calib$ComDimOPLSModel$Up,
                                R2X = if(res_calib$cv$OrthoLVsOptimalNum == 0){
                                  NULL
                                } else {
                                  res_calib$ComDimOPLSModel$R2XC[2:(res_calib$cv$OrthoLVsOptimalNum+1)] #Explained variation for predictive model components after addition of Y-orthogonal model components.
                                }
                   ),
                   R2X =  res_calib$ComDimOPLSModel$R2X[1:ndim],
                   R2Y =  res_calib$ComDimOPLSModel$R2Y[1:ndim],
                   Q2 = res_calib$cv$Q2YhatVars[ndim+res_calib$cv$OrthoLVsOptimalNum,],
                   DQ2 = if(method == 'k-OPLS-DA'){
                      res_calib$cv$DQ2Yhat
                   } else {
                      vector()
                   },
                   Singular = res_calib$ComDimOPLSModel$Sp,
                   Mean = list(MeanMB = res_calib$MeanMB,
                               MeanY = MeanY),
                   Norm = list(NormMB = res_calib$NormMB),
                   PLS.model = list(), # Parameters not provided because it cannot directly applied without the kernel transformation, as in regular ComDim-PLS
                   cv = list(Type = cvType,
                             cvTrainIndex = res_calib$cv$cvTrainIndex, # Indices for the training set observations during the cross-validation rounds.
                             cvTestIndex = res_calib$cv$cvTestIndex # Indices for the test set observations during the cross-validation rounds.
                   ),
                   Prediction = if(method == 'k-OPLS-DA'){ # To do?
                     list(Y.pred = res_calib$cv$Yhat,
                          decisionRule = decisionRule,
                          trueClass = res_calib$da$trueClass,
                          predClass = res_calib$da$predClass,
                          Sensitivity = res_calib$da$sensvec,
                          Specificity = res_calib$da$specvec,
                          confusionMatrix = res_calib$da$confusionMatrix)
                   }else {
                     list(Y.pred = res_calib$cv$Yhat)
                   },
                   Metadata = if(length(MB@Metadata) != 0){
                      MB@Metadata
                   },
                   variable.block = belong_block,
                   runtime = as.numeric(running_time))

    return(res_calib)
}




DQ2 <- function(Ypred, # Ypred contains predicted Y values
                Y # Y contains class labels (0 and 1)
                ){

# function DQ2 calculates the discriminant Q2 value adjusted for values that are larger then
# 1 for class 1 samples and smaller than 0 for 0 class samples.
# The idea is these values should not be penalized.
# NB: It is assumed y consists of ones and zeros indicating two classes

# PRESSD contains the total squared error for discrimination
# dqq = DQ2
# JAW
# BDA/SILS/UvA, 1.09.2008

# Copyright 2008 Biosystems Data Analysis Group; Universiteit van Amsterdam
# This is a license to use and modify SOFTWARE  & DATA produced by:
# THE BIOSYSTEMS DATA ANALYSIS GROUP OF THE UNIVERSITEIT VAN AMSTERDAM

# If you use, modify or redistribute the SOFTWARE & DATA and/or your source modifications, you agree:
# i)	to use the SOFTWARE & DATA and/or your source modifications solely as part of your research and not in any commercial product;
# ii)	that the SOFTWARE  & DATA and/or your source modifications will not be distributed for profit;
# iii)	all copyright notices and this license note with the SOFTWARE & DATA are retained any  redistribution of the SOFTWARE & DATA, or any portion thereof;
# iv)	to indemnify, hold harmless, and defend the "Biosystems Data Analysis Group of the Universiteit van Amsterdam" from and against any claims or lawsuits that arise or result from the use of the SOFTWARE & DATA or your source modifications.
# v)	to properly reference the SOFTWARE & DATA when used in your reseach in any publication that may result from that research.
# Reserved Rights. The "Biosystems Data Analysis Group of the Universiteit van Amsterdam" retain title and all ownership rights to the SOFTWARE & DATA.

Class0 <- which(Y==0)            # Find values belonging to Class0
Class1 <- which(Y==1)            # Find values belonging to Class1

E0 <- Ypred[Class0]-Y[Class0]   # Calculate Residuals of Class0 samples
E1 <- Ypred[Class1]-Y[Class1]   # Calculate Residuals of Class1 samples

E0count <- which(E0>0)           # Find predictions for Class0 samples larger than 0
E1count <- which(E1<0)           # Find predictions for Class1 samples larger than 1

SSE0 <- t(E0[E0count]) %*% E0[E0count] # Calculate SSE for those samples of Class0
SSE1 <- t(E1[E1count]) %*% E1[E1count] # Calculate SSE for those samples of Class1

PRESS <- t((Y - Ypred) %*% (Y - Ypred)) # PRESS for Q2

PRESSD <- SSE0 + SSE1           # Calculate total SSE for all samples = PRESSD

Ym <- Y - mean(Y)               # Calculate Total sum of squares (over all samples)
TSS <- t(Ym) %*% Ym

dqq <- 1 - (PRESSD/TSS)         # Calculate DQ2

dq <- list()
dq$dqq <- dqq
dq$PRESSD <- PRESSD
dq$Q2 <- 1 - (PRESS/TSS) # The normal Q2

return(dq)
}








# koplsCenterKTrTr <- function(K){
#   ###########################################################################
#   #
#   # Centering function for the training kernel, which is constructed
#   # from the training matrix Xtr as K = <phi(Xtr), phi(Xtr)>
#   # (see 'koplsKernel()' for details on constructing a kernel matrix).
#   #
#   # ** INPUT
#   # K = Training kernel matrix; K = <phi(Xtr), phi(Xtr)>.
#   #
#   # ** OUTPUT
#   # K = The centered kernel matrix.
#   #
#   ###########################################################################
#   #
#   # Authors: Mattias Rantalainen, Imperial College and
#   #   Max Bylesj?, Ume? University
#   # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
#   #
#   ###########################################################################
#   #
#   # This file is part of the K-OPLS package.
#   #
#   # The K-OPLS package is free software; you can redistribute it and/or
#   # modify it under the terms of the GNU General Public License version 2
#   # as published by the Free Software Foundation.
#   #
#   # The K-OPLS package is distributed in the hope that it will be useful,
#   # but WITHOUT ANY WARRANTY; without even the implied warranty of
#   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   # GNU General Public License for more details.
#   #
#   ###########################################################################
#   nrowK <- nrow(K)
#
#   I <- diag(nrowK)
#   I_n <- matrix(rep(1,nrowK), ncol = 1, nrow = nrowK)
#   K <- (I- (1/nrowK)* I_n%*%t(I_n)) %*% K%*%(I-(1/nrowK)*I_n%*%t(I_n))
#   return(K)
# }







# koplsMaxClassify <- function(data){
#   ###########################################################################
#   #
#   # Classification function that assesses class belonging of 'data'
#   # based on the maximum value.
#   #
#   # ** INPUT
#   # data = matrix containing the predicted response matrix Y,
#   #   where columns denote classes and rows observations.
#   #
#   # ** OUTPUT
#   # predClass = the predicted class(es) of 'data'.
#   #
#   ###########################################################################
#   #
#   # Authors: Mattias Rantalainen, Imperial College and
#   #   Max Bylesj?, Ume? University
#   # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
#   #
#   ###########################################################################
#   #
#   # This file is part of the K-OPLS package.
#   #
#   # The K-OPLS package is free software; you can redistribute it and/or
#   # modify it under the terms of the GNU General Public License version 2
#   # as published by the Free Software Foundation.
#   #
#   # The K-OPLS package is distributed in the hope that it will be useful,
#   # but WITHOUT ANY WARRANTY; without even the implied warranty of
#   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   # GNU General Public License for more details.
#   #
#   ###########################################################################
#   #
#
#   predClass <- vector()
#   for(i in 1:length(data[,1])){
#     tmp <- which(data[i,]== max(data[i,]))
#     if(length(tmp) == 1) {
#       predClass[i] <- which(data[i,] == max(data[i,]))
#     } else {
#       predClass[i] <- NaN
#     }
#   }
#
#   predClass <- matrix(predClass, nrow = length(predClass))
#
#   return(predClass)
#
# }


# koplsBasicClassify <- function(data,k){
# ##########################################################################
# #
# # Classification function that assesses class belonging of a predicted
# # response in 'data' based on a fixed threshold 'k'.
# #
# # ** INPUT
# # data = matrix containing the predicted response matrix Y,
# #   where columns denote classes and rows observations.
# # k = threshold value used to assign class categories.
# #
# # ** OUTPUT
# # predClass = the predicted class(es) of 'data' given 'k'.
# #
# ###########################################################################
# #
# # Authors: Mattias Rantalainen, Imperial College and
# #   Max Bylesj?, Ume? University
# # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
# #
# ###########################################################################
# #
# # This file is part of the K-OPLS package.
# #
# # The K-OPLS package is free software; you can redistribute it and/or
# # modify it under the terms of the GNU General Public License version 2
# # as published by the Free Software Foundation.
# #
# # The K-OPLS package is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.
# #
# ###########################################################################
#
# #   k is boundary, data is y_hat
#   predClass <- matrix(rep(NaN, ncol(data)*nrow(data)), nrow = nrow(data), ncol = ncol(data))
#   for(i in 1:length(data[,1])){
#     tmp <- which(data[i,] > k)
#     predClass[i,1:length(tmp)] <- tmp
#   }
#   return(predClass)
# }


koplsSensSpec <- function(classVect, predClass, y){
  ###########################################################################
  #
  # Calculates sensitivity and specificity in a class-wise fashion.
  #
  # ** INPUT
  # classVect = row vector of true class assignments (template).
  # predClass = matrix (or row vector) of class assignments to be compared.
  # y = y-matrix, used to extract the column names
  #
  # ** OUTPUT
  # sensvec = sensitivity for each class.
  # specvec = specificity for each class.
  # classvec = the class identifier corresponding to each column in
  #   sensvec and specvec.
  # meanSens = mean sensitivity over all classes.
  # meanSpec = mean specificity over all classes.
  #
  ###########################################################################
  #
  # Authors: Mattias Rantalainen, Imperial College and
  #   Max Bylesj?, Ume? University
  # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
  #
  ###########################################################################
  #
  # This file is part of the K-OPLS package.
  #
  # The K-OPLS package is free software; you can redistribute it and/or
  # modify it under the terms of the GNU General Public License version 2
  # as published by the Free Software Foundation.
  #
  # The K-OPLS package is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details.
  #
  ###########################################################################

  classy <- colnames(y)
  predClass <- classy[predClass]

  specvec <- sensvec <- rep(NA, length(classy))
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
  }

  #meanSens <- mean(sensvec)
  #meanSpec <-mean(specvec)

  # FP: I had to create the object specs to return several objects as in Matlab.
  specs <- list()
  specs$sensvec <- sensvec
  specs$specvec <- specvec
  specs$classvec <- classVect
  #specs$meanSpec <- meanSpec
  specs$confusionMatrix <- confusionMatrix

  return(specs)
}


# koplsPlotSensSpec <- function(modelFull){
#   ###########################################################################
#   # Plots sensitivity and specificity results from cross-validation
#   # in a bar plot. The produced bars are shown separately for each
#   # class including overall sensitivity and specificity results.
#   #
#   # ** INPUT
#   # modelFull = the K-OPLS cross-validation results from 'koplsCV()'
#   #
#   # ** OUTPUT
#   # res = The resulting sensitivity and specificity measures.
#   #
#   ###########################################################################
#   #
#   # Authors: Mattias Rantalainen, Imperial College and
#   #   Max Bylesj?, Ume? University
#   # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
#   #
#   ###########################################################################
#   #
#   # This file is part of the K-OPLS package.
#   #
#   # The K-OPLS package is free software; you can redistribute it and/or
#   # modify it under the terms of the GNU General Public License version 2
#   # as published by the Free Software Foundation.
#   #
#   # The K-OPLS package is distributed in the hope that it will be useful,
#   # but WITHOUT ANY WARRANTY; without even the implied warranty of
#   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   # GNU General Public License for more details.
#   #
#   ###########################################################################
#
#
#   if (modelFull$class != 'koplscv'){
#     error('Unknown model type (must be of type "koplscv"). Aborting.');
#   }
#
#   labels <- vector()
#
#   num_entries <- length(modelFull$da$meanSensAllOsc)
#   plot_mat <- rbind(rep(0,num_entries), rep(0,num_entries))
#
#   for (i in 1:num_entries){
#     plot_mat[1, i] <- modelFull$da$meanSensAllOsc[[i]]
#     plot_mat[2, i] <- modelFull$da$meanSpecAllOsc[[i]]
#     labels <- c(label, sprintf('To,%d', i-1))
#   }
#
#
#   #[sensvec, specvec, classvec, tot_sens]=koplsSensSpec(modelFull.da.trueClass, modelFull.da.predClass);
#
#   #error('The cross-validation results are not for discriminant analysis. Aborting');
#
#   h <- bar(t(plot_mat), names = labels, ylab = 'Sens. and spec. (#)',
#            title = 'Sensitivity and specificity over cross-validation')
#
#   sensvec <- plot_mat[1,]
#   specvec <- plot_mat[2,]
#   res$sens <- sensvec
#   res$spec <- specvec
#   res$tot_sens <- modelFull$da$tot_sens
#   res$meanSens <- mean(sensvec)
#   res$meanSpec <- mean(specvec)
#
#   plot(h)
#
#   return(res)
# }


# koplsConfusionMatrix <- function(true,pred){
#   ###########################################################################
#   #
#   # Calculates a confusion matrix from classification results.
#   #
#   # ** INPUT
#   # true = true class belonging.
#   # pred = predicted class assignment.
#   #
#   # ** OUTPUT
#   # A = Confusion matrix.
#   #
#   ###########################################################################
#   #
#   # Authors: Mattias Rantalainen, Imperial College and
#   #   Max Bylesj?, Ume? University
#   # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
#   #
#   ###########################################################################
#   #
#   # This file is part of the K-OPLS package.
#   #
#   # The K-OPLS package is free software; you can redistribute it and/or
#   # modify it under the terms of the GNU General Public License version 2
#   # as published by the Free Software Foundation.
#   #
#   # The K-OPLS package is distributed in the hope that it will be useful,
#   # but WITHOUT ANY WARRANTY; without even the implied warranty of
#   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   # GNU General Public License for more details.
#   #
#   ###########################################################################
#
#
#   uniqueClass <- unique(true)
#
#   A <- matrix(rep(0,length(uniqueClass)^2),
#               ncol = length(uniqueClass),
#               nrow = length(uniqueClass))
#
#   for(i in 1:length(uniqueClass)){ #for each class
#
#     indTrue <- which(true == uniqueClass[i])
#
#     for(j in 1:length(indTrue)){
#
#       A[i,which(uniqueClass==pred[indTrue[j]])] <- A[i,which(uniqueClass==pred[indTrue[j]])]+1
#
#     }
#
#     A[i,] <- A[i,]/length(indTrue)
#
#   }
#
#   return(A)
# }


# koplsScale <- function(X,YcenterType, YscaleType){
#
#   if(!(is.matrix(X))){
#     X <- as.matrix(X)
#   }
#
#   scaleS <- list()
#   scaleS$meanV <- colMeans(X)
#   scaleS$stdV <- apply(X, 2, sd)
#
#   if(YcenterType == 'mc'){
#     YScaleObj <- scale(X, center = TRUE, scale = FALSE)
#   } else if (YcenterType == 'no') {
#     YScaleObj <- X
#   } else {
#     stop('YcenterType was not set correctly. Choose "mc" or "no".')
#   }
#
#   if(YscaleType == 'uv'){
#     YScaleObj <- scale(YScaleObj, center = FALSE,
#                        scale = scaleS$stdV)
#   } else if (YscaleType == 'pareto') {
#     YScaleObj <- scale(YScaleObj, center = FALSE,
#                        scale = sqrt(scaleS$stdV))
#   } else if (YscaleType == 'no'){
#     # Do nothing
#   } else {
#     stop('YscaleType was not set correctly. Choose "uv", "pareto", or "no".')
#   }
#   scaleS$scaleType <- YscaleType
#   scaleS$centerType <- YcenterType
#   scaleS$X <- YScaleObj
#
#   return(scaleS)
# }



# koplsScaleApply <- function(X,scaleS){
#
#   scaleSA <- list()
#
#   scaleSA$centerType <- scaleS$centerType
#   scaleSA$scaleType <- scaleS$scaleType
#   scaleSA$meanV <- scaleS$meanV
#   scaleSA$stdV <- scaleS$stdV
#   ncolX <- ncol(X)
#   nrowX <- nrow(X)
#
#   if(scaleS$centerType=='mc'){
#     X <- X - matrix(rep(scaleS$meanV, nrowX), ncol = ncolX, byrow = TRUE)
#   }
#
#   if(scaleS$scaleType=='uv'){
#     X <- X / matrix(rep(scaleS$stdV, nrowX), ncol = ncolX, byrow = TRUE)
#   }
#
#   if(scaleS$scaleType=='pa'){
#     X <- X / matrix(rep(sqrt(scaleS$stdV), nrowX), ncol = ncolX, byrow = TRUE)
#   }
#
#   scaleSA$X <- X
#
#   return(scaleSA)
# }


# koplsCenterKTeTe <- function(KteTe,KteTr,KtrTr){
#   ###########################################################################
#   #
#   # Centering function for the test kernel, which is constructed
#   # from the test matrix Xte as KteTe = <phi(Xte), phi(Xte)>.
#   # Requires additional (un-centered) kernels KteTr and KteTr to
#   # estimate mean values (see 'koplsKernel()' for details on
#   # constructing a kernel matrix).
#   #
#   # ** INPUT
#   # KteTe = Test kernel matrix; KteTe = <phi(Xte), phi(Xte)>.
#   # KteTr = Test/training kernel matrix;
#   #   KteTr = <phi(Xte), phi(Xtr)>.
#   # KtrTr = Training kernel matrix; KtrTr = <phi(Xtr), phi(Xtr)>.
#   #
#   # ** OUTPUT
#   # KteTe = The centered test kernel matrix.
#   #
#   ###########################################################################
#   #
#   # Authors: Mattias Rantalainen, Imperial College and
#   #   Max Bylesj?, Ume? University
#   # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
#   #
#   ###########################################################################
#   #
#   # This file is part of the K-OPLS package.
#   #
#   # The K-OPLS package is free software; you can redistribute it and/or
#   # modify it under the terms of the GNU General Public License version 2
#   # as published by the Free Software Foundation.
#   #
#   # The K-OPLS package is distributed in the hope that it will be useful,
#   # but WITHOUT ANY WARRANTY; without even the implied warranty of
#   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   # GNU General Public License for more details.
#   #
#   ###########################################################################
#
#   nrowKer <- nrow(KteTr)
#   ncolKer <- ncol(KteTr)
#
#   Itrain <- diag(nrowKer)
#   I_nTrain <- as.matrix(rep(1,ncolKer), nrow = ncolKer, ncol = 1)
#   nTrain <- nrowKer
#
#   I <- diag(ncolKer)
#   I_n <- as.matrix(rep(1,nrowKer), nrow = nrowKer, ncol = 1)
#   n <- ncolKer
#   D_te <-  (1/nTrain) * I_n %*% t(I_nTrain)
#
#   KteTe <- KteTe - D_te %*% t(KteTr) - KteTr %*% t(D_te) + D_te %*% KtrTr %*% t(D_te)
#
#   return(KteTe)
# }
#
# koplsCenterKTeTr <- function(KteTr,KtrTr){
#   ###########################################################################
#   #
#   # Centering function for the hybrid test/training kernel, which
#   # is constructed from the test matrix Xte and the training matrix
#   # Xtr as KteTr = <phi(Xte), phi(Xtr)>. Requires additional
#   # (un-centered) training kernel to estimate mean values
#   # (see 'koplsKernel()' for details on constructing a kernel matrix).
#   #
#   # ** INPUT
#   # KteTr = Hybrid test/training kernel matrix;
#   #   KteTr = <phi(Xte), phi(Xtr)>.
#   # KtrTr = Training kernel matrix; Ktrain = <phi(Xtr), phi(Xtr)>.
#   #
#   # ** OUTPUT
#   # KteTr = The centered kernel matrix.
#   #
#   ###########################################################################
#   #
#   # Authors: Mattias Rantalainen, Imperial College and
#   #   Max Bylesj?, Ume? University
#   # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
#   #
#   ###########################################################################
#   #
#   # This file is part of the K-OPLS package.
#   #
#   # The K-OPLS package is free software; you can redistribute it and/or
#   # modify it under the terms of the GNU General Public License version 2
#   # as published by the Free Software Foundation.
#   #
#   # The K-OPLS package is distributed in the hope that it will be useful,
#   # but WITHOUT ANY WARRANTY; without even the implied warranty of
#   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   # GNU General Public License for more details.
#   #
#   ###########################################################################
#
#   nrowKrr <- nrow(KtrTr)
#   ncolKrr <- ncol(KtrTr)
#   nrowKer <- nrow(KteTr)
#   ncolKer <- ncol(KteTr)
#
#
#   Itrain <- diag(nrowKrr)
#   I_nTrain <- as.matrix(rep(1,nrowKrr), nrow = nrowKrr, ncol = 1)
#   nTrain <- nrowKrr
#
#   I <- diag(nrowKer)
#   I_n <- as.matrix(rep(1,nrowKer), nrow = nrowKer, ncol = 1)
#   n <- nrowKer
#
#   KteTr <- (KteTr-(1/nTrain) * I_n %*% t(I_nTrain) %*% KtrTr) %*% (Itrain-(1/nTrain)*I_nTrain%*%t(I_nTrain))
#   return(KteTr)
# }

ComDimOPLSModel <- function(data,Y,ndim,max.ort){
  ###########################################################################
  #
  # ComDim-OPLS modelling - JB 2012
  #
  # Based on KOPLSModel.m from the K-OPLS package
  #
  # Function for training a ComDim OPLS model. The function constructs a
  # predictive regression model for predicting the values of 'Y by using
  # the information in the collection of tables.
  # The explained variation is separated into predictive components
  # (dimensionality is determined by the parameter 'A') and
  # 'Y'-orthogonal components dimensionality determined by the parameter 'nox').
  #
  # ** INPUT
  # K = Kernel matrix (un-centered); K = <phi(Xtr),phi(Xtr)>.
  # (K == collection in ComDim)
  # Y = Response matrix (un-centered/scaled).
  #
  # ndim = Number of predictive components !!!!!!
  # OR
  # # ndim = Number of Y vectors
  #
  # nort = Number of Y-orthogonal components.
  # plotdisp = 1 display - 0 no display
  #
  # ** OUTPUT
  # model = Object with the following entries:
  #   Cp = Y loading matrix.
  #   Sp = Sigma matrix, containing singular values from Y'*K*Y used
  #       for scaling.
  #   Sps = Sp^(-1/2).
  #   Up = Y score matrix.
  #   Tp = Predictive score matrix for all Y-orthogonal components.
  #   T = Predictive score matrix for the final model.
  #   co = Y-orthogonal loading vectors.
  #   so = Eigenvalues from estimation of Y-orthogonal loading vectors.
  #   To = Y-orthogonal score matrix.
  #   toNorm = Norm of the Y-orthogonal score matrix prior to scaling.
  #   Bt = T-U regression coefficients for predictions.
  #
  # A = Number of predictive components (integer) !!!!!!
  # OR
  # # A = Number of Y vectors (intteger)
  #
  #   nox = Number of Y-orthogonal components.
  #   K = The kernel matrix.
  #   EEprime = The deflated kernel matrix for residual statistics.
  #   sstot_K = Total sums of squares in 'K'.
  #   R2X = Cumulative explained variation for all model components.
  #   R2XO = Cumulative explained variation for Y-orthogonal
  #       model components.
  #   R2XC = Explained variation for predictive model components after
  #       addition of Y-orthogonal model components.
  #   sstot_Y = Total sums of squares in Y.
  #   R2Y = Explained variation of Y.
  #   preProc = Pre-processing parameters:
  #  	  K = Pre-processing setting for K = 'preProcK'.
  #	  Y = Pre-processing setting for Y = 'preProcY'.
  #	  paramsY = Scaling parameters for Y.
  #
  ###########################################################################
  #
  # Authors: Mattias Rantalainen, Imperial College and
  #   Max Bylesj?, Ume? University
  # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
  #
  ###########################################################################
  #
  # This file is part of the K-OPLS package.
  #
  # The K-OPLS package is free software; you can redistribute it and/or
  # modify it under the terms of the GNU General Public License version 2
  # as published by the Free Software Foundation.
  #
  # The K-OPLS package is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details.
  #
  ###########################################################################

  if(class(data) == "MultiBlock"){
    block_names <- getBlockNames(data)
    data2 <- list()
    for(i in 1:length(data@Data)){
      datax <- list()
      datax$Data <- as.matrix(data@Data[[i]])
      data2[[ block_names[i] ]] <- datax
    }
    data <- data2
    rm(data2)

    Y2 <- list()
    Y2$X <- Y
    Y <- Y2
    rm(Y2)
  }


  ntable <- length(data)
  nrowMB <- nrow(data[[1]]$Data)
  threshold <- 1E-10

  LAMBDA <- matrix(rep(0,ntable*max.ort), nrow = ntable, ncol = max.ort)
  W_mat <- matrix(rep(0,nrowMB*nrowMB), nrow = nrowMB, ncol = nrowMB)
  # This is the traditional implementation of ComDim...

  # centering and rescaling the tables ***************************
  for (i in 1:ntable){
    average <- colMeans(data[[i]]$Data)
    data[[i]]$Data <- data[[i]]$Data - t(t(rep(1,nrowMB))) %*% t(as.matrix(average))
    data[[i]]$Data <- data[[i]]$Data/norm(data[[i]]$Data,'2')
    temp <- data[[i]]$Data %*% t(data[[i]]$Data)
    W_mat <- W_mat+(1/ntable)*temp
  }

  K <- list()
  K[[1]] <- list()
  K[[1]][[1]] <- W_mat
  W_mat2 <- W_mat

  collection2 <- data

  #--- koplsScale(Y, 'mc', 'no')
  if(!(is.matrix(Y$X))){
    Y$X <- as.matrix(Y$X)
  }
  scaleParams <- list()
  scaleParams$meanV <- colMeans(Y$X)
  scaleParams$stdV <- apply(Y$X, 2, sd)
  scaleParams$X <- scale(Y$X, center = TRUE, scale = FALSE)
  #---
  Y <- scaleParams$X

  ## initiate Yorth related vars
  to <- list()
  co <- list()
  so <- vector()
  toNorm <- vector()

  ## ComDim-OPLS model estimation ------------

  #step 1
  # ndim = Number of predictive components (integer) !!!!!!
  # OR
  ## ndim = Number of Y vectors (integer)
  #
  csv <- svd(t(Y)%*%K[[1]][[1]]%*%Y) # Basis of PLS2
  #Cpt <- csv$u
  #Spt <- csvd
  #Vt <- csv$v
  Cp0 <- csv$u[,1:ndim]
  Sp0 <- csv$d[1:ndim]

  #step2
  Up0 <- Y %*% Cp0
  if(ndim != 1){
    Tp0 <- t(K[[1]][[1]]) %*% Up0 %*% diag(Sp0^(-1/2))
  } else {
    Tp0 <- t(K[[1]][[1]]) %*% Up0 %*% as.matrix(Sp0^(-1/2))
  }

  Bt0 <- pracma::pinv(t(Tp0) %*% Tp0) %*% t(Tp0) %*% Up0
  Lambda_Temp <- vector()

  # step 3
  Cp <- list()
  Sp <- list()
  Up <- list()
  Tp <- list()
  Bt <- list()

  for (i in 1:max.ort){

    lambda <- (1/ntable)*rep(ntable,1)
    previousfit <- 1000
    deltafit <- 1

    while (deltafit>threshold){
      K[[1]][[i]] <- W_mat
      if(length(K) < i){
        K[[i]] <- list()
      }
      K[[i]][[i]] <- W_mat2

      csv <- svd(t(Y) %*% K[[1]][[i]] %*% Y)
      Cp[[i]] <- csv$u[,1:ndim]
      Sp[[i]] <- csv$d[1:ndim]
      Up[[i]] <- Y %*% Cp[[i]]

      #step4
      if(ndim != 1){
        Tp[[i]] <- t(K[[1]][[i]]) %*% Up[[i]] %*% diag(Sp[[i]]^(-1/2))
      } else {
        Tp[[i]] <- t(K[[1]][[i]]) %*% Up[[i]] %*% as.matrix(Sp[[i]]^(-1/2))
      }

      Bt[[i]] <- pracma::pinv(t(Tp[[i]]) %*% Tp[[i]]) %*% t(Tp[[i]]) %*% Up[[i]]

      #step5
      csvo <- svd(t(Tp[[i]]) %*% (K[[i]][[i]]-Tp[[i]] %*% t(Tp[[i]])) %*% Tp[[i]])
      co[[i]] <- csvo$u[,1]
      so[i] <- csvo$d[1]


      #step6
      to[[i]]  <- (K[[i]][[i]] - Tp[[i]] %*% t(Tp[[i]])) %*% Tp[[i]] %*% as.matrix(co[[i]]) %*% t(as.matrix(so[i]))^(-1/2)

      #step7
      toNorm[[i]] <- sqrt(t(to[[i]]) %*% to[[i]])

      #step8
      to[[i]] <- to[[i]]/toNorm[[i]]

      #Saliences
      fit <- 0
      for (ta in 1:ntable) {
        temp <- collection2[[ta]]$Data %*% t(collection2[[ta]]$Data)
        lambda[ta] <- t(to[[i]]) %*% temp %*% to[[i]]
        aux <- temp - lambda[ta] * (to[[i]] %*% t(to[[i]]))
        fit <- fit + sum(sum(aux*aux))
      }

      W_mat2 <- W_mat <- matrix(rep(0,nrowMB*nrowMB),nrow = nrowMB)

      for (ta in 1:ntable){
        temp <- data[[ta]]$Data %*% t(data[[ta]]$Data)
        W_mat <- W_mat+lambda[ta]*temp
        temp2 <- collection2[[ta]]$Data %*% t(collection2[[ta]]$Data)
        W_mat2 <- W_mat2 + lambda[ta] * temp2
      }

      Lambda_Temp <- c(Lambda_Temp,lambda)
      deltafit <- previousfit-fit
      previousfit <- fit

    }
    LAMBDA[,i] <- lambda

    #step9 & 10 - Deflation
    II <- diag(rep(1,nrowMB))
    aux <- II-to[[i]] %*% t(to[[i]])
    for (ta in 1:ntable){
      data[[ta]]$Data <- aux %*% data[[ta]]$Data
      collection2[[ta]]$Data <- t(t(aux %*% collection2[[ta]]$Data) %*% aux)
    }

  }

  #step 11
  W_mat2 <- W_mat <- matrix(rep(0,nrowMB*nrowMB),nrow = nrowMB)
  for (ta in 1:ntable){
    temp <- data[[ta]]$Data %*% t(data[[ta]]$Data)
    W_mat <- W_mat+temp
    temp2 <- collection2[[ta]]$Data %*% t(collection2[[ta]]$Data)
    W_mat2 <- W_mat2+temp2
  }
  K[[1]][[max.ort+1]] <- W_mat
  K[[max.ort+1]] <- list()
  K[[max.ort+1]][[max.ort+1]] <- W_mat2
  csv <- svd(t(Y) %*% K[[1]][[max.ort+1]] %*% Y)
  Cp[[max.ort+1]] <- csv$u[,1:ndim]
  Sp[[max.ort+1]] <- csv$d[1:ndim]
  Up[[max.ort+1]] <- Y %*% Cp[[max.ort+1]]

  #step12
  if(ndim != 1){
    Tp[[max.ort+1]] <- t(K[[1]][[max.ort+1]]) %*% Up[[max.ort+1]] %*% diag(Sp[[max.ort+1]]^(-1/2))
  } else {
    Tp[[max.ort+1]] <- t(K[[1]][[max.ort+1]]) %*% Up[[max.ort+1]] %*% as.matrix(Sp[[max.ort+1]]^(-1/2))
  }


  #Step13
  Bt[[i+1]] <- inv(t(Tp[[max.ort+1]]) %*% Tp[[max.ort+1]]) %*% t(Tp[[max.ort+1]]) %*% Up[[max.ort+1]]

  #---------- extra stuff -----------------
  EEprime <- (K[[max.ort+1]][[max.ort+1]]-Tp[[max.ort+1]] %*% t(Tp[[max.ort+1]]))
  sstot_K <- (sum(diag(K[[1]][[1]])))
  R2XC <- R2XO <- R2X <- vector()

  for (i in 1:(max.ort+1)){
    rss <- sum(diag( K[[i]][[i]]-Tp[[i]] %*% t(Tp[[i]])))
    R2X <- c(R2X, 1 - rss/sstot_K)

    rssc <- sum(diag(K[[1]][[1]] - Tp[[i]] %*% t(Tp[[i]])))
    R2XC <- c(R2XC, 1 - rssc/sstot_K)

    rsso <- sum(diag( K[[i]][[i]]))
    R2XO <- c(R2XO, 1 - rsso/sstot_K)
  }

  # should work but not fully tested (MB 2007-02-19)
  sstot_Y <- sum(sum(Y*Y))
  FF <- Y-Up[[1]] %*% t(Cp[[1]])
  R2Y <- 1 - sum(sum(FF*FF))/sstot_Y
  #----------------------------------------
  #

  lambda_pred <- matrix(rep(0,ndim*ntable), ncol = ntable, nrow = ndim)
  for (ta in 1:ntable){
    for (j in 1:ndim){
      temp <- data[[ta]]$Data %*% t(data[[ta]]$Data)
      lambda_pred[j,ta] <- t(Tp[[max.ort+1]][,j]) %*% temp %*% Tp[[max.ort+1]][,j]
    }
  }

  saliences <- cbind(t(lambda_pred), LAMBDA)

  dimS <- dim(saliences)
  for (j in 1:dimS[2]){
    saliences[,j] <- saliences[,j]/sum(saliences[,j])
  }

  if(!is.null(names(data))){
    rownames(saliences) <- names(data)
  }
  colnames(saliences) <- paste0('CC',1:ncol(saliences))

  model <- list()
  model$Cp <- Cp
  model$Sp <- Sp
  model$Up <- Up
  model$Tp <- Tp
  model$Tpo <- Tp[[max.ort+1]] # In matlab, this is model.T
  model$co <- co
  model$so <- so
  model$to <- to
  model$toNorm <- toNorm
  model$Bt <- Bt
  model$ndim <- ndim
  model$max.ort <- max.ort
  model$K <- K
  model$saliences <- saliences

  #convert struct to matrix
  model$To <- matrix(rep(0, nrow(model$Tpo)*max.ort), ncol = max.ort, nrow = nrow(model$Tpo))
  for (i in 1:length(to)){
    model$To[,i] <- to[[i]]
  }

  #extra stuff
  model$EEprime <- EEprime
  model$sstot_K <- sstot_K
  model$R2X <- R2X
  model$R2XO <- R2XO
  model$R2XC <- R2XC
  model$sstot_Y <- sstot_Y
  model$R2Y <- R2Y
  model$Up0 <- Up0
  model$Cp0 <- Cp0
  model$Sp0 <- Sp0
  model$Tp0 <- Tp0
  model$Bt0 <- Bt0

  #Pre-processing
  model$preProc$K <- 'mc'
  model$preProc$Y <- 'mc'
  model$preProc$paramsY <- scaleParams

  model$class <- 'kopls'

  # if (plotdisp == 1){
  #
  #   iterations <- cbind.data.frame(Iterations = seq(1,length(Lambda_Temp)),
  #                                  Saliences = unlist(Lambda_Temp), # Need to check this is correct.
  #                                  Block = rep(rep(1:ntable),length(Lambda_Temp)))
  #
  #   plot1 <- ggplot(data= iterations, aes(x = Iterations,
  #                                         y = Saliences,
  #                                         group= Block)) +
  #     geom_line()+
  #     geom_point()
  #
  #
  #   plot2 <- ggplot(data = iterations, aes(x = Iterations,
  #                                          y = Block)) +
  #     geom_tile(aes(fill = Saliences))
  #
  #   p <- ggarrange(plot1, plot2, labels = LETTERS[1:2], ncol = 2, nrow = 1)
  #
  #   print(p)
  # }

  return(model)
}



ComDimOPLSPredict <- function(KteTr,Ktest,Ktrain,model,max.ort = NULL,rescaleY =0){

  ###########################################################################
  #
  # Performs prediction of new samples from an existing K-OPLS model
  # (see koplsModel()' to calculate K-OPLS models).
  # The function projects the Y-predictive and Y-orthogonal scores
  # components to predict a value of the response matrix Y.
  # The dimensionality of the parameters is determined from the
  # specified model.
  #
  # ** INPUT
  # KteTr = The hybrid test/training kernel matrix;
  #   KteTr = <phi(Xte),phi(Xtr)>.
  # Ktest = The pure test kernel matrix;
  #   Ktest = <phi(Xte),phi(Xte)>.
  # Ktrain = The training kernel matrix (same as used in
  #   model training); Ktrain = <phi(Xtr),phi(Xtr)>.
  # model = K-OPLS model object.
  # nort = Number of Y-orthogonal components. If not specified, the
  #   number used during model training will be employed.
  # rescaleY = Boolean parameter. If true, predicted values of the
  #   response (Yhat) is rescaled according to the pre-processing
  #   settings of the model. If false, Yhat is not rescaled (default).
  #
  # ** OUTPUT:
  # modelp = Object with the following entries:
  #   Tp = Predicted predictive score matrix for all generations
  #       0:'nort' of Y-orthogonal vectors.
  #   T = Predictive score matrix for the final model with 'nort'
  #       Y-orthogonal vectors.
  #   to = Predicted Y-orthogonal score vectors.
  #   EEprime = Calculated residuals for the test kernel 'Ktest',
  #       useful e.g. for residual statistics.
  #   Yhat = Predicted values of the response matrix.
  #
  ###########################################################################
  #
  # Authors: Mattias Rantalainen, Imperial College and
  #   Max Bylesj?, Ume? University
  # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
  #
  ###########################################################################
  #
  # This file is part of the K-OPLS package.
  #
  # The K-OPLS package is free software; you can redistribute it and/or
  # modify it under the terms of the GNU General Public License version 2
  # as published by the Free Software Foundation.
  #
  # The K-OPLS package is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details.
  #
  ###########################################################################

  if (model$class != 'kopls') {
    stop('Model must be of type "kopls". Aborting.')
  }


  # mean centering of K matrices

  # centering (order of these arg is important...)
  KteTeMc <- Ktest
  if (model$preProc$K == 'mc'){
    #---KteTeMc <- koplsCenterKTeTe(Ktest,KteTr,Ktrain)
    nrowKer <- nrow(KteTr)
    ncolKer <- ncol(KteTr)
    Itrain <- diag(nrowKer)
    I_nTrain <- as.matrix(rep(1,ncolKer), nrow = ncolKer, ncol = 1)
    nTrain <- nrowKer
    I <- diag(ncolKer)
    I_n <- as.matrix(rep(1,nrowKer), nrow = nrowKer, ncol = 1)
    n <- ncolKer
    D_te <-  (1/nTrain) * I_n %*% t(I_nTrain)
    KteTeMc <- Ktest - D_te %*% t(KteTr) - KteTr %*% t(D_te) + D_te %*% Ktrain %*% t(D_te)
    rm(nrowKer, ncolKer, Itrain, I_nTrain, I, I_n, n, D_te)
  }

  KteTe <- list()
  KteTe[[1]] <- list()
  KteTe[[1]][[1]] <- KteTeMc

  KteTrMc <- KteTr
  if (model$preProc$K == 'mc'){
    #KteTrMc <- koplsCenterKTeTr(KteTr,Ktrain)
    nrowKrr <- nrow(Ktrain)
    ncolKrr <- ncol(Ktrain)
    nrowKer <- nrow(KteTr)
    ncolKer <- ncol(KteTr)
    Itrain <- diag(nrowKrr)
    I_nTrain <- as.matrix(rep(1,nrowKrr), nrow = nrowKrr, ncol = 1)
    nTrain <- nrowKrr
    I <- diag(nrowKer)
    I_n <- as.matrix(rep(1,nrowKer), nrow = nrowKer, ncol = 1)
    n <- nrowKer
    KteTrMc <- (KteTr-(1/nTrain) * I_n %*% t(I_nTrain) %*% Ktrain) %*% (Itrain-(1/nTrain)*I_nTrain%*%t(I_nTrain))
    rm(nrowKrr, ncolKrr, nrowKer, ncolKer, I_nTrain, nTrain, I, I_n, n)
  }

  KteTr <- list()
  KteTr[[1]] <- list()
  KteTr[[1]][[1]] <- KteTrMc

  ## init of Y-orth comps
  to <- list()

  ## check if last arg is number of components to use in prediction:

  if(is.null(max.ort)){
    max.ort <- model$max.ort
  } else if(max.ort > model$max.ort){
    warning('Number of Y-orthogonal components to use is higher than in model - setting number of Yorth to max in model');
    max.ort <- model$max.ort
  }

  ## KOPLS prediction
  Tp <- list()

  if(max.ort != 0){
    for(i in 1:max.ort){ #step1

      #step2
      xx <- length(unlist(model$Sp[[i]]^(-1/2)))
      if(xx != 1){
        Tp[[i]] <- KteTr[[i]][[1]]%*%model$Up[[i]] %*% diag(model$Sp[[i]]^(-1/2))
      } else {
        Tp[[i]] <- KteTr[[i]][[1]]%*%model$Up[[i]] %*% as.matrix(model$Sp[[i]]^(-1/2))
      }
      #Yhat[[i]] <- Tp{i}*model.Bt{i}*model.Cp';

      #step3
      to[[i]] <- (KteTr[[i]][[i]]-Tp[[i]]%*%t(model$Tp[[i]]))%*%model$Tp[[i]]%*%model$co[[i]] %*% model$so[[i]]^(-1/2)

      #step4
      to[[i]] <- to[[i]]/model$toNorm[[i]]

      #step 4.5 deflate KteTe. (this is an EXTRA feature - not in alg. in paper )
      if(length(KteTe) < i+1){
         KteTe[[i+1]] <- list()
      }
      KteTe[[i+1]][[i+1]] <- KteTe[[i]][[i]] -
        KteTr[[i]][[i]] %*% model$to[[i]] %*% t(to[[i]]) -
        to[[i]] %*% t(model$to[[i]]) %*% t(KteTr[[i]][[i]]) +
        to[[i]] %*% t(model$to[[i]]) %*% model$K[[i]][[i]] %*% model$to[[i]] %*% t(to[[i]])

      #step5

      if(length(KteTr) < i+1){
        KteTr[[i+1]] <- list()
      }
      KteTr[[i+1]][[1]] <- KteTr[[i]][[1]] - to[[i]] %*% t(model$to[[i]]) %*% t(model$K[[1]][[i]])

      #step6
      KteTr[[i+1]][[i+1]] <- KteTr[[i]][[i]] - KteTr[[i]][[i]]%*%model$to[[i]]%*%t(model$to[[i]])-to[[i]]%*%t(model$to[[i]])%*%model$K[[i]][[i]] + to[[i]] %*% t(model$to[[i]]) %*% model$K[[i]][[i]]%*%model$to[[i]]%*%t(model$to[[i]])
    }
  }

  #step7
  if(max.ort==0){
    i <- 0
  }

  xx <- length(unlist(model$Sp0^(-1/2)))
  if(xx != 1){
    Tp[[i+1]] <- KteTr[[i+1]][[1]] %*% model$Up0 %*% diag(model$Sp0^(-1/2))
  } else {
    Tp[[i+1]] <- KteTr[[i+1]][[1]] %*% model$Up0 %*% as.matrix(model$Sp0^(-1/2))
  }

  Yhat <- Tp[[i+1]] %*% model$Bt0 %*% t(model$Cp0)


  #---- Extra stuff ----------------------------------
  #this appears to be correct - but does not match previous code...
  EEprime <- KteTe[[i+1]][[i+1]]-Tp[[i+1]]%*%t(Tp[[i+1]])
  #--------------------------------------------------

  modelp <- list()
  modelp$Tp <- Tp
  modelp$to <- to
  modelp$EEprime <- EEprime
  modelp$Yhat <- Yhat

  return(modelp)
}

# koplsRescale <- function(scaleS = NULL,varargin = NULL){
#   ###########################################################################
#   #
#   # Scales a matrix based on pre-defined parameters from a scaling
#   # object.
#   #
#   # ** INPUT
#   # scaleS = An object containing scaling parameters
#   #   (see 'koplsScale()').
#   # varargin = If defined, this matrix will be scaled and returned.
#   #	Otherwise the original data set in the scaleS object will be
#   #   scaled and returned.
#   #
#   # ** OUTPUT
#   # scaleS = An object containing the following entries:
#   #   centerType = 'mc' (mean-centering) or 'no' (no centering).
#   #   scaleType = 'uv' (unit variance), 'pa' (pareto) or 'no'
#   #       (no scaling).
#   #   meanV = vector with mean values for all columns in X.
#   #   stdV = vector with standard deviations for all columns in X.
#   #   X = Scaled version of 'varargin', if defined, otherwise,
#   #       scaled version of scaleS.X from input. Scaling is done
#   #       according to 'centerType' and 'scaleType'.
#   #
#   ###########################################################################
#   #
#   # Authors: Mattias Rantalainen, Imperial College and
#   #   Max Bylesj?, Ume? University
#   # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
#   #
#   ###########################################################################
#   #
#   # This file is part of the K-OPLS package.
#   #
#   # The K-OPLS package is free software; you can redistribute it and/or
#   # modify it under the terms of the GNU General Public License version 2
#   # as published by the Free Software Foundation.
#   #
#   # The K-OPLS package is distributed in the hope that it will be useful,
#   # but WITHOUT ANY WARRANTY; without even the implied warranty of
#   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   # GNU General Public License for more details.
#   #
#   ###########################################################################
#
#   if(!(is.null(varargin))){
#     X <- varargin
#   }  else {
#     X <- scaleS$X
#   }
#
#   mX <- nrow(X)
#   nX <- ncol(X)
#
#
#   if(scaleS$scaleType=='uv'){
#     X <- X * matrix(rep(scaleS$stdV, mX), nrow = mX, byrow = TRUE)
#   }
#
#   if(scaleS$scaleType=='pareto') {
#     X <- X * matrix(rep(sqrt(scaleS$stdV), mX), nrow = mX, byrow = TRUE)
#   }
#
#   if(scaleS$centerType=='mc') {
#     X <- X + matrix(rep(1,mX), nrow = mX, ncol = 1) %*% scaleS$meanV
#   }
#
#   scaleS$centerType <- 'no'
#   scaleS$scaleType <- 'no'
#   scaleS$X <- X
#
#   return(scaleS)
# }


# koplsDummy <- function(class, numClasses = length(unique(class))){
#   ###########################################################################
#   #
#   # Converts integer vector to binary matrix (dummy matrix).
#   #
#   # ** INPUT
#   # class = vector with class belongings (integer).
#   # numClasses = pre-defines the number of classes in the output
#   #	(if undefined, the number of unique entries in 'class' will be
#   #   used).
#   #
#   # ** OUTPUT
#   # dummy = A matrix with rows corresponding to observations and
#   #   columns to classes. Each element in matrix is either one
#   #   (observation belongs to class) or zero (observation does not
#   #   belong to class).
#   # labels_sorted = The class labels that are found in class in
#   #	sorted order.
#   #
#   ###########################################################################
#   #
#   # Authors: Mattias Rantalainen, Imperial College and
#   #   Max Bylesj?, Ume? University
#   # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
#   #
#   ###########################################################################
#   #
#   # This file is part of the K-OPLS package.
#   #
#   # The K-OPLS package is free software; you can redistribute it and/or
#   # modify it under the terms of the GNU General Public License version 2
#   # as published by the Free Software Foundation.
#   #
#   # The K-OPLS package is distributed in the hope that it will be useful,
#   # but WITHOUT ANY WARRANTY; without even the implied warranty of
#   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   # GNU General Public License for more details.
#   #
#   ###########################################################################
#
#   labels <- unique(class) #the set of classlabels in class
#   labels_sorted <- sort(labels) #sort labels in ascending order
#
#   len_class <- length(class)#number of samples
#
#   dummy <- matrix(rep(0,len_class*numClasses), nrow = len_class, ncol = numClasses)
#   #dummy matrix initialized as a zero matrix
#
#   for (i in 1:numClasses){ #for each class label
#     ind <- which(class==labels_sorted[i])
#     #find the rows (samples) that belongs to the current class, labels_sorted(i
#     dummy[ind,i] <- 1 #!!!write ones in the positions where we have the current class....
#   }
#
#   return(dummy)
# }


# koplsReDummy <-function(Y){
#   ###########################################################################
#   #
#   # Reconstructs a (integer) class vector from a binary (dummy) matrix.
#   #
#   # ** INPUT
#   # Y = Dummy matrix. See 'koplsDummy()' for details.
#   #
#   # ** OUTPUT
#   # classVect = The reconstructed integer class vector.
#   #
#   ###########################################################################
#   #
#   # Authors: Mattias Rantalainen, Imperial College and
#   #   Max Bylesj?, Ume? University
#   # Copyright (c) 2007-2008 Mattias Rantalainen and Max Bylesj?
#   #
#   ###########################################################################
#   #
#   # This file is part of the K-OPLS package.
#   #
#   # The K-OPLS package is free software; you can redistribute it and/or
#   # modify it under the terms of the GNU General Public License version 2
#   # as published by the Free Software Foundation.
#   #
#   # The K-OPLS package is distributed in the hope that it will be useful,
#   # but WITHOUT ANY WARRANTY; without even the implied warranty of
#   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   # GNU General Public License for more details.
#   #
#   ###########################################################################
#
#   if(!is.matrix(Y)){
#     Y <- as.matrix(Y)
#   }
#
#   classVect <- rep(NaN, nrow(Y))
#   for(i in 1:ncol(Y)){
#     ind <- which(Y[,i]==1)
#     classVect[ind] <- i
#   }
#
#   return(classVect)
# }
