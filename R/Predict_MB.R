Predict_MB <- function(MB = MB, y, model = model,
                       normalise = FALSE, loquace = TRUE) {


  #' ComDim - Finding common dimensions in multi-block datasets
  #'
  #' Predict_MB applies an existent ComDim model to a MultiBlock object.
  #' @param MB A MultiBlock object for which prediction is desired.
  #' @param y An y block to calculate DQ2 in ComDim-PLS and ComDim-OPLS models (optional).
  #' @param model The ComDim model used for prediction.
  #' @param normalise To apply normalisation. FALSE == no (default), TRUE == yes.
  #' @param loquace If TRUE, it returns a message with the ComDim elements to be predicted (i.e. Q.scores, T.scores...)
  #' @return A ComDim object.

  #' @export

  if(class(MB) != 'MultiBlock'){
    stop("'MB' is  not a MultiBlock.")
  }

  if(class(model) != 'ComDim'){
    stop("'model' is  not class ComDim.")
  }

  if(!(all(names(model@T.scores) %in% getBlockNames(MB)))){
    stop('Not all the model blocks are included in MB.')
  }

  # Delete blocks not included in the model.
  if(any(!(getBlockNames(MB) %in% names(model@T.scores)))){
    dif_names <- setdiff(getBlockNames(MB),names(model@T.scores))
    for(i in rev(dif_names)){
      MB@Data[[i]] <- NULL
      MB@Variables[[i]] <- NULL
    }
  }


  for(i in names(model@T.scores)){

    if(all(rownames(model@P.loadings)[model@variable.block == i] %in% MB@Variables[[i]])){
      MB@Data[[i]] <- MB@Data[[i]][,match(rownames(model@P.loadings)[model@variable.block == i], MB@Variables[[i]]) ]
      MB@Variables[[i]] <- MB@Variables[[i]][match(rownames(model@P.loadings)[model@variable.block == i], MB@Variables[[i]])]
    } else {
      stop('Some variables are missing in MB.')
    }
  }

  if(hasArg(y)){

    if(model@Method == 'PLS-DA' || model@Method == 'OPLS-DA' ||
       model@Method == 'PLS-R' || model@Method == 'OPLS-R'){

      if(is.vector(y)){
        if(length(y) != length(MB@Samples)){
          stop("'y' does not have the appropiate length.")
        }
      } else {
        if(nrow(y) != length(MB@Samples)){
          stop("'y' must be of the same length than MB.")
        }
      }
    }

    if(model@Method == 'PLS-DA' || model@Method == 'OPLS-DA'){

      tmp <- unique(as.vector(y))
      if(all(tmp %in% c(0,1))){ # Is it dummy?

        if(!is.matrix(y)){
          y <- as.matrix(y)
        }

        classVect <- rep(NA, nrow(y))

        if(ncol(y) == 1){
          classVect <- as.vector(y)
        } else {

          if(is.null(colnames(y))){
            stop("Column names for 'y' must be provided.")
          }

          if(any(rowSums(y) != 1)){
            if(any(rowSums(y) == 0)){
              stop('At least one sample was not assigned to a class.')
            } else if (any(colSums(y) > 1)){
              stop('At least one sample was assigned to more than one class.')
            }
          }
        }

        for(i in 1:ncol(y)){
          classVect[which(y[,i] == 1)] <- colnames(y)[i]
        }

      } else { # If not dummy
        if((is.matrix(y) || is.data.frame(y)) && ncol(y) > 1){
          stop('ComDim-PLS can only be applied if Y is a dummy matrix or a class vector')
        }

        # If y is not dummy (need to convert to dummy)
        classy <- sort(unique(as.vector(y)))
        if(any(!(classy %in% colnames(model@PLS.model$Y)))){
          stop("Not all the classes in 'y' can be predicted with the given model.")
        }
        #classVect <- colnames(model@PLS.model$Y)
        classVect <- as.vector(y)

        # Now generate the dummy matrix from y
        y <- matrix(rep(0, length(classy)*length(classVect)),
                    ncol = length(classy), nrow =length(classVect))
        for(i in 1:length(classy)){
          pos <- which(classVect == classy[i])
          if(length(pos) != 0){
            y[pos,i] <- 1
          }
        }
        colnames(y) <- classy
        rm(classy)
      }

    } else if(model@Method == 'PLS-R' || model@Method == 'OPLS-R'){
      if((is.matrix(y) || is.data.frame(y)) && ncol(y) > 1){
        stop('ComDim-PLS can only be applied if Y is a dummy matrix or a class vector')
      } else if(is.vector(y)){
        y <- as.matrix(y)
      }
    }
  }

  # Create concatenated MB and normalised as in the original data
  names_table = names(model@T.scores) # block names
  nrowMB = length(getSampleNames(MB))

  Xnorm <- list()
  Xnorm_mat <- matrix(, nrow = nrowMB, ncol = nrow(model@P.loadings))
  T_Loc <- list()
  W <- matrix(, nrow = nrowMB, ncol = 0)

  for(i in names_table){ # I use block names instead of integers because I don't know if block order is kept in the new MB.
    if (normalise){
    # Normalise original blocks
      X_mean <- MB@Data[[i]] - matrix(data = rep(1,nrowMB), ncol = 1, nrow = nrowMB) %*% model@Mean$MeanMB[[i]]
      XX <- X_mean*X_mean
      Norm_X <- sqrt(sum(XX))
      X_Normed <- X_mean/Norm_X
      Xnorm[[i]] <- X_Normed
      Xnorm_mat[,model@variable.block == i] <- Xnorm[[i]]

    } else {
      Xnorm[[i]] <- MB@Data[[i]]
      Xnorm_mat[,model@variable.block == i] <- Xnorm[[i]]
    }
  }
  print(Xnorm_mat[1:3,1:3])



  Q <- matrix(, ncol = model@ndim, nrow  = nrowMB)

  # Application of the models to the new MBs.
  if(model@Method == 'PCA' || model@Method == 'PLS-DA' || model@Method == 'PLS-R') {

    for(j in 1:model@ndim){

      for(i in names_table){
        # Q.d are orthonormal in ICA & PCA

        temp <- t(t(model@P.loadings[model@variable.block == i,j])) # Local loadings

        if(i == names_table[1]){
          Q[,j] <- pracma::pinv(t(Xnorm[[i]])) %*% temp
        }

        if(j == 1) {
          T_Loc[[i]] <- matrix(,ncol = model@ndim, nrow = nrowMB)
        }

        T_Loc[[i]][,j] <- Xnorm[[i]]%*%(temp %*% pracma::pinv(t(temp)%*%temp)) # local Scores

        # Deflate each temp_tabCalib
        Xnorm[[i]] <- Xnorm[[i]]-Q[,j]%*%t(temp)
      }
    }

    model@Q.scores <- Q
    model@T.scores <- T_Loc
    cat("'Q.scores' and 'T.scores' were predicted from the ComDim model.\n")

  } else if (model@Method == 'OPLS-DA' || model@Method == 'OPLS-R'){

    # Substract orthogonal part
    Q.ort <- NULL # Initialization
    T_Loc_o <- list() # Initialization
    if(model@Orthogonal$nort != 0){
      Q.ort <- matrix(, ncol = model@Orthogonal$nort, nrow = nrowMB)
      for(j in 1:model@Orthogonal$nort){

        for(i in names_table){

          tempo <- t(t(model@Orthogonal$P.loadings.ort[model@variable.block == i,j])) # Local loadings

          if(i == names_table[1]){
            Q.ort[,j] <- pracma::pinv(t(Xnorm[[i]])) %*% tempo
          }

          if(j == 1) {
            T_Loc_o[[i]] <- matrix(,ncol = model@Orthogonal$nort, nrow = nrowMB)
          }

          T_Loc_o[[i]][,j] <- Xnorm[[i]]%*%(tempo %*% pracma::pinv(t(tempo)%*%tempo)) # local Scores

          # Deflate each temp_tabCalib
          Xnorm[[i]] <- Xnorm[[i]]-Q.ort[,j]%*%t(tempo)
        }
      }
      for(i in 1:length(names_table)){
        if(i == 1){
          Xnorm_mat <- Xnorm[[i]]
        } else {
          Xnorm_mat <- cbind(Xnorm_mat, Xnorm[[i]])
        }
      }
    }
    model@Orthogonal$Q.scores.ort <- Q.ort
    model@Orthogonal$T.scores.ort <- T_Loc_o
    cat("'Q.scores.ort' and 'T.scores.ort' were predicted from the ComDim model.\n")

    # Substract predictive part
    for(j in 1:model@ndim){

      for(i in names_table){
        # Q.d are orthonormal in ICA & PCA

        temp <- t(t(model@P.loadings[model@variable.block == i,j])) # Local loadings

        if(i == names_table[1]){
          Q[,j] <- pracma::pinv(t(Xnorm[[i]])) %*% temp
        }

        if(j == 1) {
          T_Loc[[i]] <- matrix(,ncol = model@ndim, nrow = nrowMB)
        }

        T_Loc[[i]][,j] <- Xnorm[[i]]%*%(temp %*% pracma::pinv(t(temp)%*%temp)) # local Scores

        # Deflate each temp_tabCalib
        Xnorm[[i]] <- Xnorm[[i]]-Q[,j]%*%t(temp)
      }
    }
    model@Q.scores <- Q
    model@T.scores <- T_Loc
    cat("'Q.scores' and 'T.scores' were predicted from the ComDim model.\n")


  } else {
    stop('The method used in ComDim is not a valid ComDim method.\n')
  }

  if(model@Method == 'PLS-DA' || model@Method == 'PLS-R' ||
     model@Method == 'OPLS-DA' || model@Method == 'OPLS-R'){

    # Calculate y predicted
    Ypred <- Xnorm_mat %*% model@PLS.model$B

    for(i in 1:ncol(Ypred)){
      Ypred[,i] <- unlist(Ypred[,i]) + model@PLS.model$B0[i]
    }

    if(hasArg(y)){

      # Calculate Q2
      if(model@Method == 'PLS-R' || model@Method == 'OPLS-R'){
        PRESS <- t(Ypred-y) %*% (Ypred-y)
        Ym <- y - mean(y)
        TSS <- t(Ym) %*% Ym
        model@Prediction$Q2 <- 1 - PRESS/TSS
        cat("'Q2' was predicted from the ComDim model.\n")
      } else if(model@Method == 'PLS-DA' || model@Method == 'OPLS-DA'){

        predClass <- rep(NA, nrow(y))

        if(model@Prediction$decisionRule == 'max'){

          for(i in 1:nrow(Ypred)){
            xx <- which(Ypred[i,] == max(Ypred[i,], na.rm = TRUE))
            if(length(xx) == 1) {
              predClass[i] <- colnames(y)[xx] # Y stores the Class names
            } else {
              predClass[i] <- NaN # There is a tie.
              warning(sprintf('Class for sample %s could not be predicted.\n', i))
            }
          }

        } else if(model@Prediction$decisionRule == 'fixed'){
          for(i in 1:nrow(Ypred)){
            xx <- which(Ypred[,i] > 1/ncol(y))
            if(length(xx) == 1) {
              predClass[i] <- colnames(y)[xx] # Y stores the Class names
            } else {
              predClass[i] <- NaN # There is a tie.
              warning(sprintf('Class for sample %s could not be predicted.\n', i))
            }
          }
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
          colnames(cm) <- c("predClass1", "predClass0")
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
        model@Q2 <- Q2
        model@DQ2 <- DQ2
        model@Prediction$trueClass = classVect
        model@Prediction$predClass = predClass
        model@Prediction$Sensitivity = sensvec
        model@Prediction$Specificity = specvec
        model@Prediction$confusionMatrix = confusionMatrix
        cat("Q2, DQ2, classVect, predClass, Sensitivity, Specificity, and confusionMatrix were predicted from the ComDim model.\n")
      }
    }

    # for(i in 1:ncol(Ypred)){
    #   Ypred[,i] <- Ypred[,i] * model@Mean$MeanY[i] # Re-scale Ypred.
    # }
    model@Prediction$Y.pred <- Ypred
    cat("'Y.pred' was predicted from the ComDim model.\n")

  }

  return(model)
}
