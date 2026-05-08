#' ComDim_OPLS
#'
#' Finding common dimensions in multi-block datasets using OPLS. Also known as ConsensusOPLS
#' (ComDim-OPLS) for multiblock structures: orthogonal components uncorrelated
#' with Y are extracted from all blocks simultaneously before the predictive
#' components are computed.
#'
#' This function is a wrapper around \code{\link[ConsensusOPLS]{ConsensusOPLS}}.
#' The core kernel-OPLS extraction is delegated to that package; all ComDim
#' output slots (local scores, loadings, VIP, sensitivity, confusion matrix,
#' etc.) are computed from the returned model objects.
#' @param MB A MultiBlock object.
#' @param y The Y-block. A class vector or dummy matrix for OPLS-DA, or a
#'   numeric matrix/vector for OPLS-R.
#' @param ndim Number of predictive Common Dimensions. Default is 1.
#' @param nort Maximum number of orthogonal Common Dimensions. Default is 1.
#'   The actual number used is determined by ConsensusOPLS cross-validation and
#'   may be less than this value.
#' @param method 'OPLS-DA' for discriminant analysis or 'OPLS-R' for regression.
#' @param decisionRule Only used if method is 'OPLS-DA'. If 'fixed', samples are
#'   assigned to the class with Y-hat above 1/nclasses. If 'max', samples are
#'   assigned to the class with the highest Y-hat.
#' @param normalise To apply block normalisation. FALSE == no (default), TRUE == yes.
#' @param loquace To display the calculation times. TRUE == yes, FALSE == no (default).
#' @param cv.k Number of folds for k-fold cross-validation (default 7). Set
#'   to 0 to skip CV output. ConsensusOPLS always performs internal CV to select
#'   the optimal number of orthogonal components; when \code{cv.k >= 2} the
#'   resulting \code{Q2} and \code{DQ2} reflect that cross-validation.
#' @return A \code{ComDim} object.  All slots are populated.  Key slots:
#' \describe{
#'   \item{\code{Method}}{\code{"OPLS-DA"} or \code{"OPLS-R"}.}
#'   \item{\code{ndim}}{Number of predictive Common Dimensions.}
#'   \item{\code{Q.scores}}{Predictive global scores matrix (\eqn{n \times ndim}).}
#'   \item{\code{T.scores}}{Named list of block-specific predictive local scores.}
#'   \item{\code{P.loadings}}{Global predictive loadings.}
#'   \item{\code{Saliences}}{Predictive block salience matrix (\eqn{ntable \times ndim}).}
#'   \item{\code{Orthogonal}}{List with orthogonal component outputs: \code{nort},
#'     \code{Q.scores}, \code{T.scores}, \code{P.loadings.ort}, \code{Saliences.ort}.}
#'   \item{\code{R2X}}{Named vector (length \eqn{ndim + nort}) of X-variance fractions.}
#'   \item{\code{R2Y}}{Named vector (length \eqn{ndim + nort}) of Y-variance fractions.}
#'   \item{\code{Q2}}{Cross-validated Q2 per class/response (when \code{cv.k >= 2};
#'     otherwise training-set fit).}
#'   \item{\code{DQ2}}{(OPLS-DA only) Cross-validated discriminant Q2 per class.}
#'   \item{\code{VIP}}{Global total VIP (named vector, length \eqn{p_{tot}}).}
#'   \item{\code{VIP.block}}{Named list (one data.frame per block) with columns
#'     \code{p}, \code{o}, \code{tot}.}
#'   \item{\code{PLS.model}}{KOPLS regression objects: \code{W}, \code{B}, \code{B0}, \code{Y}.}
#'   \item{\code{cv}}{Cross-validation results when \code{cv.k >= 2}: \code{k},
#'     \code{Ypred}, \code{Q2}, \code{DQ2}.}
#'   \item{\code{Prediction}}{Training-set predictions: \code{Y.pred}; for OPLS-DA
#'     also \code{decisionRule}, \code{trueClass}, \code{predClass},
#'     \code{Sensitivity}, \code{Specificity}, \code{confusionMatrix}.}
#'   \item{\code{Mean}}{List with \code{MeanMB} and \code{MeanY}.}
#'   \item{\code{Norm}}{List with \code{NormMB}, \code{FrobNorms}, \code{RVweights}.}
#'   \item{\code{variable.block}}{Block membership of each variable.}
#'   \item{\code{runtime}}{Total computation time in seconds.}
#' }
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50)
#' b2 <- matrix(rnorm(800), 10, 80)
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#' y <- rep(c("A", "B"), 5)
#' results <- ComDim_OPLS(mb, y, ndim = 1, nort = 1, method = "OPLS-DA")
#' @references Boccard J, Rutledge DN (2013).
#'   A consensus OPLS-DA strategy for multiblock Omics data fusion.
#'   \emph{Analytica Chimica Acta}, 769, 30--39.
#'   \doi{10.1016/j.aca.2013.01.022}
#' @export
ComDim_OPLS <- function(MB = MB, y = y, ndim = 1, nort = 1,
                        method = c("OPLS-DA", "OPLS-R"),
                        decisionRule = c("fixed", "max")[2],
                        normalise = FALSE, #threshold = 1e-10,
                        loquace = FALSE,
                        cv.k = 7) {

  if (!requireNamespace("ConsensusOPLS", quietly = TRUE)) {
    stop("Package 'ConsensusOPLS' is required. Install with: install.packages('ConsensusOPLS')")
  }
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("The package pracma is needed.")
  }

  # --- Input validation ---
  if (!(method == "OPLS-DA" || method == "OPLS-R")) {
    stop("'method' must be 'OPLS-DA' or 'OPLS-R'.")
  }
  if (nort < 1) stop("'nort' must be at least 1.")
  if (ndim < 1) stop("'ndim' must be at least 1.")

  progress_bar <- utils::txtProgressBar(min = 0, max = 80, style = 3, char = "=")
  start_ini_time <- Sys.time()

  if (!inherits(MB, "MultiBlock")) stop("'MB' is not a MultiBlock.")

  ntable         <- length(blockNames(MB))
  nrowMB         <- length(sampleNames(MB))
  variable_number <- ncol(MB)
  give_error     <- 0

  # --- Zero-variance removal ---
  numvar   <- 0
  numblock0 <- vector()
  for (i in 1:ntable) {
    var0 <- which(apply(MB@Data[[i]], 2, var) == 0)
    if (length(var0) != 0) {
      numvar <- numvar + length(var0)
      if (length(var0) == length(MB@Variables[[i]])) numblock0 <- append(numblock0, i)
      MB@Data[[i]]      <- MB@Data[[i]][, setdiff(1:variable_number[i], var0)]
      MB@Variables[[i]] <- MB@Variables[[i]][setdiff(1:variable_number[i], var0)]
      variable_number[i] <- length(MB@Variables[[i]])
    }
  }
  if (numvar != 0)
    warning(sprintf("Number of variables excluded from the analysis because of their zero variance: %d", numvar))
  if (length(numblock0) != 0) {
    numblock0 <- rev(numblock0)
    for (i in seq_along(numblock0)) {
      MB@Data[[numblock0[i]]] <- NULL
      MB@Variables[[numblock0[i]]] <- NULL
      if (names(MB@Data)[numblock0[i]] %in% names(MB@Batch))    MB@Batch[[numblock0[i]]]    <- NULL
      if (names(MB@Data)[numblock0[i]] %in% names(MB@Metadata)) MB@Metadata[[numblock0[i]]] <- NULL
    }
    warning(sprintf("Number of blocks excluded from the analysis because of their zero variance: %d",
                    length(numblock0)))
  }

  ntable          <- length(blockNames(MB))
  variable_number <- ncol(MB)

  if (length(MB@Data) == 0) { warning("All blocks had 0 variance"); give_error <- 1 }

  for (i in 1:ntable) {
    if (any(is.na(MB@Data[[i]])))       stop("The MB contains NA values. They must be removed first for ComDim analysis.")
    if (any(is.infinite(MB@Data[[i]]))) stop("The MB contains infinite values. They must be removed first for ComDim analysis.")
  }
  if (give_error) stop("The data is not ready for ComDim.") else message("The data can be used for ComDim.")

  # --- Y preparation ---
  if (method == "OPLS-DA") {
    tmp <- unique(as.vector(y))
    if (all(tmp %in% c(0, 1))) {
      if (!is.matrix(y)) y <- as.matrix(y)
      classVect <- rep(NA, nrow(y))
      if (ncol(y) == 1) {
        classVect <- as.vector(y)
      } else {
        if (any(rowSums(y) != 1)) {
          if (any(rowSums(y) == 0))   stop("At least one sample was not assigned to a class.")
          else if (any(colSums(y) > 1)) stop("At least one sample was assigned to more than one class.")
        }
      }
      if (is.null(colnames(y))) colnames(y) <- paste0("Class", seq_len(ncol(y)))
      for (i in seq_len(ncol(y))) classVect[which(y[, i] == 1)] <- colnames(y)[i]
    } else {
      if ((is.matrix(y) || is.data.frame(y)) && ncol(y) > 1)
        stop("ComDim-OPLS can only be applied if Y is a dummy matrix or a class vector")
      classVect        <- as.vector(y)
      classVect_sorted <- sort(unique(classVect))
      y <- matrix(0, nrow = length(classVect), ncol = length(classVect_sorted))
      for (i in seq_along(classVect_sorted)) y[which(classVect == classVect_sorted[i]), i] <- 1
      colnames(y) <- classVect_sorted
      rm(classVect_sorted)
    }
  }
  if (!is.matrix(y)) y <- as.matrix(y)

  pieceBar <- 80 / (4 + 2 * ntable + ndim + nort)
  total_progress <- pieceBar

  DimLabels   <- paste0("CC",  seq_len(ndim))
  OrtLabels   <- paste0("ort", seq_len(nort))
  TableLabels <- blockNames(MB)

  if (loquace)
    message(sprintf("Initialisation finished after : %s millisecs",
                    (Sys.time() - start_ini_time) * 1000))
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  # --- Block preparation (normalisation) ---
  temp_tabCalib <- vector("list", ntable)
  MeanMB  <- list()
  NormMB  <- list()
  meanY   <- colMeans(y)

  for (i in 1:ntable) {
    MeanMB[[TableLabels[i]]] <- colMeans(MB@Data[[i]])
    names(MeanMB[[TableLabels[i]]]) <- MB@Variables[[i]]

    if (normalise) {
      X_mean  <- MB@Data[[i]] - matrix(rep(1, nrowMB), ncol = 1) %*% MeanMB[[TableLabels[i]]]
      Norm_X  <- sqrt(sum(X_mean^2))
      NormMB[[TableLabels[i]]] <- setNames(rep(Norm_X, variable_number[i]), MB@Variables[[i]])
      temp_tabCalib[[i]] <- X_mean / Norm_X
    } else {
      NormMB[[TableLabels[i]]] <- setNames(rep(1, variable_number[i]), MB@Variables[[i]])
      temp_tabCalib[[i]] <- MB@Data[[i]]
    }

    total_progress <- total_progress + pieceBar
    utils::setTxtProgressBar(progress_bar, value = total_progress)
  }

  # Frobenius norms of block kernels (needed for B_kopls)
  frob_norms <- setNames(numeric(ntable), TableLabels)
  for (i in seq_len(ntable)) {
    Xi <- temp_tabCalib[[i]]
    if (normalise) { ss <- sum(Xi^2); if (ss > 1e-14) Xi <- Xi / sqrt(ss) }
    Ki      <- tcrossprod(Xi)
    frob_Ki <- sqrt(sum(Ki^2))
    frob_norms[i] <- if (frob_Ki > 1e-14) frob_Ki else 1
  }

  # --- Build data list for ConsensusOPLS ---
  data_for_cOpls <- setNames(temp_tabCalib, TableLabels)
  for (i in seq_len(ntable)) {
    rownames(data_for_cOpls[[i]]) <- MB@Samples
    colnames(data_for_cOpls[[i]]) <- MB@Variables[[i]]
  }

  cOpls <- ConsensusOPLS::ConsensusOPLS(
    data      = data_for_cOpls,
    Y         = if (method == "OPLS-DA") classVect else y,
    maxPcomp  = ndim,
    maxOcomp  = nort,
    modelType = if (method == "OPLS-DA") "da" else "reg",
    nperm     = 0L,
    cvType    = "nfold",
    nfold     = if (cv.k >= 2L) as.integer(cv.k) else 5L,
    verbose   = loquace
  )

  # Use component counts actually selected by ConsensusOPLS
  ndim <- cOpls@nPcomp
  nort <- cOpls@nOcomp
  DimLabels <- paste0("CC",  seq_len(ndim))
  OrtLabels <- paste0("ort", seq_len(nort))

  # Extract koplsModel internals
  kModel      <- cOpls@model$koplsModel
  Tp          <- kModel$scoresP
  Q_ort       <- kModel$scoresO
  Cp          <- kModel$Cp
  Up          <- kModel$Up
  Sp          <- if (is.matrix(kModel$Sp)) kModel$Sp else diag(as.vector(kModel$Sp), nrow = ndim)
  Bt          <- kModel$Bt[[length(kModel$Bt)]]   # regression after all ort deflations
  rv_weights  <- setNames(cOpls@model$RV, TableLabels)
  normKernels <- cOpls@model$normKernels

  colnames(Tp)    <- DimLabels; rownames(Tp)    <- MB@Samples
  colnames(Q_ort) <- OrtLabels; rownames(Q_ort) <- MB@Samples
  colnames(Cp)    <- DimLabels; rownames(Cp)    <- colnames(y)

  total_progress <- total_progress + pieceBar * ntable
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  # --- Saliences from normalised block kernels ---
  saliences     <- matrix(NA_real_, nrow = ntable, ncol = ndim,
                          dimnames = list(TableLabels, DimLabels))
  saliences_ort <- matrix(NA_real_, nrow = ntable, ncol = nort,
                          dimnames = list(TableLabels, OrtLabels))
  for (i in seq_len(ntable)) {
    Ki <- normKernels[[i]]
    for (j  in seq_len(ndim))  saliences[i, j]      <- as.vector(crossprod(Tp[, j],     Ki) %*% Tp[, j])
    for (jj in seq_len(nort))  saliences_ort[i, jj]  <- as.vector(crossprod(Q_ort[, jj], Ki) %*% Q_ort[, jj])
  }

  # --- Orthogonal local scores, orthogonal loadings, and ort-deflated X ---
  # Save pre-deflation blocks for VIP
  X_orig_blocks <- temp_tabCalib

  T_Loc_ort <- setNames(lapply(seq_len(ntable), function(i)
    matrix(NA_real_, nrow = nrowMB, ncol = nort,
           dimnames = list(MB@Samples, OrtLabels))), TableLabels)
  orthoP <- NULL

  for (jj in seq_len(nort)) {
    for (i in seq_len(ntable)) {
      tempo <- t(temp_tabCalib[[i]]) %*% Q_ort[, jj]
      T_Loc_ort[[TableLabels[i]]][, jj] <-
        temp_tabCalib[[i]] %*% (tempo %*% pracma::pinv(t(tempo) %*% tempo))
      temp_tabCalib[[i]] <- temp_tabCalib[[i]] - Q_ort[, jj] %*% t(tempo)
      tempo2 <- if (i == 1) tempo else append(tempo2, tempo)
    }
    orthoP <- if (jj == 1) as.matrix(tempo2) else cbind(orthoP, tempo2)
  }
  colnames(orthoP) <- OrtLabels

  # Rebuild orthogonally-deflated Xnorm_mat
  Xnorm_mat <- do.call(cbind, temp_tabCalib)

  # --- Predictive local scores ---
  T_Loc <- setNames(lapply(seq_len(ntable), function(i)
    matrix(NA_real_, nrow = nrowMB, ncol = ndim,
           dimnames = list(MB@Samples, DimLabels))), TableLabels)

  for (j in seq_len(ndim)) {
    q_j_ss   <- sum(Tp[, j]^2)
    q_j_unit <- if (q_j_ss > 1e-14) Tp[, j] / sqrt(q_j_ss) else Tp[, j]
    for (i in seq_len(ntable)) {
      temp <- t(temp_tabCalib[[i]]) %*% q_j_unit
      T_Loc[[TableLabels[i]]][, j] <-
        temp_tabCalib[[i]] %*% (temp %*% pracma::pinv(t(temp) %*% temp))
      temp_tabCalib[[i]] <- temp_tabCalib[[i]] - q_j_unit %*% t(temp)
    }
  }

  varnames <- unlist(lapply(seq_len(ntable), function(i) MB@Variables[[i]]))

  L_CD_Vec <- t(Xnorm_mat) %*% Tp
  colnames(L_CD_Vec) <- DimLabels; rownames(L_CD_Vec) <- varnames
  rownames(orthoP)   <- varnames

  total_progress <- total_progress + pieceBar
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  # --- B_kopls (variable-space regression coefficients) ---
  UpSp   <- Up %*% Sp
  W_pred <- matrix(NA_real_, nrow = sum(variable_number), ncol = ndim)
  k_idx  <- 1L
  for (i in seq_len(ntable)) {
    vi   <- variable_number[i]
    Xb_ort <- Xnorm_mat[, k_idx:(k_idx + vi - 1L), drop = FALSE]
    W_pred[k_idx:(k_idx + vi - 1L), ] <- t(Xb_ort) %*% UpSp * rv_weights[i] / frob_norms[i]
    k_idx <- k_idx + vi
  }
  colnames(W_pred) <- DimLabels; rownames(W_pred) <- varnames
  B_kopls  <- W_pred %*% Bt %*% t(Cp)
  B0_kopls <- colMeans(y) - colMeans(Xnorm_mat) %*% B_kopls

  # --- R2X and R2Y from koplsModel cumulative vectors ---
  # kModel$R2X:  cumulative X-variance explained by predictive components
  # kModel$R2XO: cumulative X-variance explained by orthogonal components
  # kModel$R2Yhat: cumulative Y-variance explained by predictive components
  r2x_pc <- kModel$R2X[seq_len(ndim)]
  r2x_oc <- kModel$R2XO[seq_len(nort)]
  r2x <- setNames(c(c(r2x_pc[1], if (ndim > 1) diff(r2x_pc)),
                    c(r2x_oc[1], if (nort > 1) diff(r2x_oc))),
                  c(DimLabels, OrtLabels))

  r2y_pc <- kModel$R2Yhat[seq_len(ndim)]
  r2y <- setNames(c(c(r2y_pc[1], if (ndim > 1) diff(r2y_pc)),
                    rep(0, nort)),
                  c(DimLabels, OrtLabels))

  SingVal            <- setNames(colSums(Tp^2), DimLabels)
  Sum_saliences_Dim  <- setNames(colSums(saliences), DimLabels)
  Sum_saliences_Tab  <- setNames(rowSums(saliences), TableLabels)

  total_progress <- total_progress + pieceBar * (1 + ndim)
  utils::setTxtProgressBar(progress_bar, value = min(total_progress, 79))

  # --- VIP (same formula as the direct implementation) ---
  .vip_for_scores <- function(scores_mat, Y_raw, X_b, p_b) {
    ncomp <- ncol(scores_mat)
    ss_T  <- colSums(scores_mat^2)
    if (any(ss_T < .Machine$double.eps)) return(rep(NA_real_, p_b))
    Qs <- crossprod(Y_raw, scores_mat) %*% diag(1 / ss_T, nrow = ncomp)
    ss_Q <- colSums(Qs^2)
    if (any(ss_Q < .Machine$double.eps)) return(rep(NA_real_, p_b))
    Us <- Y_raw %*% Qs %*% diag(1 / ss_Q, nrow = ncomp)
    ss_U <- colSums(Us^2)
    if (any(ss_U < .Machine$double.eps)) return(rep(NA_real_, p_b))
    Ws <- crossprod(X_b, Us) %*% diag(1 / ss_U, nrow = ncomp)
    Ws <- apply(Ws, 2, function(w) { nrm <- sqrt(sum(w^2)); if (nrm > 1e-14) w / nrm else w })
    s_a   <- ss_T * ss_Q
    sum_s <- sum(s_a)
    if (sum_s < .Machine$double.eps) return(rep(NA_real_, p_b))
    as.vector(sqrt(p_b * (Ws^2 %*% s_a) / sum_s))
  }

  .vip_ort_for_loadings <- function(ort_scores_mat, ort_loadings_b) {
    p_b  <- nrow(ort_loadings_b)
    ss_T <- colSums(ort_scores_mat^2)
    if (any(ss_T < .Machine$double.eps)) return(rep(NA_real_, p_b))
    P_norm <- apply(ort_loadings_b, 2, function(pp) {
      nrm <- sqrt(sum(pp^2)); if (nrm > 1e-14) pp / nrm else pp
    })
    if (is.null(dim(P_norm))) P_norm <- matrix(P_norm, ncol = 1)
    sum_s <- sum(ss_T)
    if (sum_s < .Machine$double.eps) return(rep(NA_real_, p_b))
    as.vector(sqrt(p_b * (P_norm^2 %*% ss_T) / sum_s))
  }

  VIP_blocks <- setNames(vector("list", ntable), TableLabels)
  k_ort <- 1L
  for (i in seq_len(ntable)) {
    p_b     <- variable_number[i]
    vip_p   <- .vip_for_scores(Tp, y, X_orig_blocks[[i]], p_b)
    ort_P_b <- orthoP[k_ort:(k_ort + p_b - 1L), , drop = FALSE]
    vip_o   <- .vip_ort_for_loadings(Q_ort, ort_P_b)
    vip_t   <- sqrt((vip_p^2 + vip_o^2) / 2)
    VIP_blocks[[TableLabels[i]]] <- data.frame(p = vip_p, o = vip_o, tot = vip_t,
                                               row.names = MB@Variables[[i]])
    k_ort <- k_ort + p_b
  }
  VIP_global <- setNames(unlist(lapply(VIP_blocks, `[[`, "tot")), varnames)

  # --- variable.block ---
  belong_block <- character(sum(variable_number))
  k <- 1L
  for (i in seq_len(ntable)) {
    belong_block[k:(k + variable_number[i] - 1L)] <- TableLabels[i]
    k <- k + variable_number[i]
  }

  # --- Training-set Y prediction ---
  Ypred <- Xnorm_mat %*% B_kopls +
    matrix(as.vector(B0_kopls), nrow = nrowMB, ncol = ncol(y), byrow = TRUE)
  rownames(Ypred) <- MB@Samples
  if (!is.null(colnames(y))) colnames(Ypred) <- colnames(y)

  # --- Class prediction, sensitivity, specificity, confusion matrix ---
  if (method == "OPLS-DA") {
    predClass <- rep(NA_character_, nrow(y))
    if (decisionRule == "max") {
      for (i in seq_len(nrow(Ypred))) {
        xx <- which(Ypred[i, ] == max(Ypred[i, ], na.rm = TRUE))
        if (length(xx) == 1) predClass[i] <- colnames(y)[xx]
        else { predClass[i] <- NaN; warning(sprintf("Class for sample %s could not be predicted.", i)) }
      }
    } else if (decisionRule == "fixed") {
      for (i in seq_len(nrow(Ypred))) {
        xx <- which(Ypred[, i] > 1 / ncol(y))
        if (length(xx) == 1) predClass[i] <- colnames(y)[xx]
        else { predClass[i] <- NaN; warning(sprintf("Class for sample %s could not be predicted.", i)) }
      }
    }

    classy <- colnames(y)
    Q2 <- DQ2 <- specvec <- sensvec <- setNames(rep(NA_real_, length(classy)), classy)
    confusionMatrix <- list()

    for (i in seq_along(classy)) {
      pos <- which(classVect == classy[i])
      neg <- which(classVect != classy[i])
      tp  <- length(which(classVect == classy[i] & predClass == classy[i]))
      tn  <- length(which(classVect != classy[i] & predClass != classy[i]))
      if (length(tp) == 0) tp <- 0
      if (length(tn) == 0) tn <- 0

      sensvec[i] <- tp / length(pos)
      specvec[i] <- tn / length(neg)

      cm <- matrix(c(tp,
                     length(which(classVect != classy[i] & predClass == classy[i])),
                     length(which(classVect == classy[i] & predClass != classy[i])),
                     tn), 2, 2)
      rownames(cm) <- c("predClass1", "predClass0")
      colnames(cm) <- c("trueClass1", "trueClass0")
      confusionMatrix[[classy[i]]] <- cm

      E0   <- Ypred[neg, i] - y[neg, i]
      E1   <- Ypred[pos, i] - y[pos, i]
      Ym   <- y[, i] - mean(y[, i])
      TSS  <- sum(Ym^2)
      DQ2[i] <- 1 - (sum(E0[E0 > 0]^2) + sum(E1[E1 < 0]^2)) / TSS
      Q2[i]  <- 1 - (sum(E0^2) + sum(E1^2)) / TSS
    }
  } else {
    Ym <- y - mean(y)
    Q2 <- 1 - sum((Ypred - y)^2) / sum(Ym^2)
  }

  # Override training-set Q2/DQ2 with cross-validated values when available
  if (cv.k >= 2) {
    q2_cv  <- cOpls@Q2
    dq2_cv <- if (method == "OPLS-DA") cOpls@DQ2 else NULL
    if (!is.null(q2_cv)  && length(q2_cv)  > 0) Q2  <- q2_cv
    if (!is.null(dq2_cv) && length(dq2_cv) > 0) DQ2 <- dq2_cv
  }

  # --- CV result slot ---
  if (cv.k >= 2) {
    Ypred_cv <- NULL
    cv_obj   <- cOpls@cv
    if (!is.null(cv_obj) && !is.null(cv_obj$AllYhat)) {
      yhat_list <- cv_obj$AllYhat
      last_yhat <- if (is.list(yhat_list)) yhat_list[[length(yhat_list)]] else yhat_list
      if (!is.null(cv_obj$cvTestIndex) && is.matrix(last_yhat)) {
        test_order <- unlist(cv_obj$cvTestIndex)
        if (length(test_order) == nrowMB) {
          Ypred_cv <- matrix(NA_real_, nrow = nrowMB, ncol = ncol(y),
                             dimnames = list(MB@Samples, colnames(y)))
          Ypred_cv[test_order, ] <- last_yhat
        }
      }
      if (is.null(Ypred_cv)) Ypred_cv <- last_yhat
    }
    cv_result <- list(
      k     = cv.k,
      Ypred = Ypred_cv,
      Q2    = Q2,
      DQ2   = if (method == "OPLS-DA") DQ2 else NULL
    )
  } else {
    if (cv.k != 0)
      warning("'cv.k' must be >= 2 or 0 (skip CV). No cross-validation performed.")
    cv_result <- list()
  }

  # --- Assemble ComDim object ---
  end_output   <- Sys.time()
  running_time <- as.numeric(end_output - start_ini_time)
  if (loquace) message(sprintf("Analysis finished after : %s millisecs", running_time * 1000))

  progress_bar <- utils::txtProgressBar(min = 0, max = 80, style = 3, char = "=")
  utils::setTxtProgressBar(progress_bar, value = 80)
  close(progress_bar)

  metadata <- list()
  if (length(MB@Metadata) != 0) metadata[[1]] <- MB@Metadata

  new("ComDim",
    Method     = method,
    ndim       = ndim,
    Q.scores   = Tp,
    T.scores   = T_Loc,
    P.loadings = L_CD_Vec,
    Saliences  = saliences,
    Orthogonal = list(
      nort            = nort,
      Q.scores        = Q_ort,
      T.scores        = T_Loc_ort,
      P.loadings.ort  = orthoP,
      Saliences.ort   = saliences_ort
    ),
    VIP        = VIP_global,
    VIP.block  = VIP_blocks,
    R2X        = r2x,
    R2Y        = r2y,
    Q2         = Q2,
    DQ2        = if (method == "OPLS-DA") DQ2 else vector(),
    Singular   = SingVal,
    Mean       = list(MeanMB = MeanMB, MeanY = meanY),
    Norm       = list(NormMB = NormMB, FrobNorms = frob_norms, RVweights = rv_weights),
    PLS.model  = list(W = W_pred, B = B_kopls, B0 = B0_kopls, Y = y),
    cv         = cv_result,
    Prediction = if (method == "OPLS-DA") {
      list(
        Y.pred          = Ypred,
        decisionRule    = decisionRule,
        trueClass       = classVect,
        predClass       = data.frame(predClass = predClass, row.names = MB@Samples),
        Sensitivity     = sensvec,
        Specificity     = specvec,
        confusionMatrix = confusionMatrix
      )
    } else {
      list(Y.pred = Ypred)
    },
    Metadata       = metadata,
    variable.block = belong_block,
    runtime        = running_time
  )
}


# ==============================================================================
# ARCHIVED DIRECT IMPLEMENTATION
# The code below is the original self-contained KOPLS implementation of
# ComDim_OPLS, kept here in case ConsensusOPLS becomes unavailable.
# To restore it: uncomment the function definition and comment out the wrapper
# above.
# ==============================================================================

# ComDim_OPLS <- function(MB = MB, y = y, ndim = 1, nort = 1,
#                         method = c("OPLS-DA", "OPLS-R"),
#                         decisionRule = c("fixed", "max")[2],
#                         normalise = FALSE, threshold = 1e-10,
#                         loquace = FALSE,
#                         cv.k = 7) {
#   # INITIALISATION
#
#   if (!(method == "OPLS-DA" || method == "OPLS-R")) {
#     stop("'method' must be 'OPLS-DA' or 'OPLS-R'.")
#   }
#
#   if (nort < 1) {
#     stop("'nort' must be at least 1.")
#   }
#
#   if (ndim < 1) {
#     stop("'ndim' must be at least 1.")
#   }
#
#   progress_bar <- utils::txtProgressBar(min = 0, max = 80, style = 3, char = "=")
#
#   start_ini_time <- Sys.time() # To start counting calculation time for the initialization
#
#   if (!inherits(MB, "MultiBlock")) {
#     stop("'MB' is  not a MultiBlock.")
#   }
#
#   ntable <- length(blockNames(MB)) # number of blocks
#   nrowMB <- length(sampleNames(MB)) # number of samples
#   variable_number <- ncol(MB) # Number of variables per block
#
#   give_error <- 0
#
#   # Remove all variables with 0 variance. If there are no remaining variables, give_error
#   numvar <- 0
#   numblock0 <- vector()
#   for (i in 1:ntable) {
#     var0 <- apply(MB@Data[[i]], 2, function(x) {
#       var(x)
#     })
#     var0 <- which(var0 == 0)
#
#     if (length(var0) != 0) {
#       numvar <- numvar + length(var0)
#       if (length(var0) == length(MB@Variables[[i]])) {
#         numblock0 <- append(numblock0, i)
#       }
#       MB@Data[[i]] <- MB@Data[[i]][, setdiff(1:variable_number[i], var0)]
#       MB@Variables[[i]] <- MB@Variables[[i]][setdiff(1:variable_number[i], var0)]
#       variable_number[i] <- length(MB@Variables[[i]])
#     }
#   }
#   if (numvar != 0) {
#     warning(sprintf("Number of variables excluded from the analysis because of their zero variance: %d", numvar))
#   }
#   if (length(numblock0) != 0) {
#     numblock0 <- rev(numblock0) # To sort in decreasing order.
#     for (i in seq_along(numblock0)) {
#       MB@Data[[numblock0[i]]] <- NULL
#       MB@Variables[[numblock0[i]]] <- NULL
#       if (names(MB@Data)[numblock0[i]] %in% names(MB@Batch)) {
#         MB@Batch[[numblock0[i]]] <- NULL
#       }
#       if (names(MB@Data)[numblock0[i]] %in% names(MB@Metadata)) {
#         MB@Metadata[[numblock0[i]]] <- NULL
#       }
#     }
#     warning(sprintf("Number of blocks excluded from the analysis because of their zero variance: %d", length(numblock0)))
#   }
#
#   ntable <- length(blockNames(MB)) # number of blocks
#   variable_number <- ncol(MB) # In case the number of blocks has changed.
#
#   if (length(MB@Data) == 0) { # In case all blocks had 0 variance...
#     warning("All blocks had 0 variance")
#     give_error <- 1
#   }
#
#   # Check for infinite or missing values (they cannot be handled with svd)
#   for (i in 1:ntable) {
#     if (any(is.na(MB@Data[[i]]))) {
#       stop("The MB contains NA values. They must be removed first for ComDim analysis.")
#     }
#     if (any(is.infinite(MB@Data[[i]]))) {
#       stop("The MB contains infinite values. They must be removed first for ComDim analysis.")
#     }
#   }
#
#   if (give_error) {
#     stop("The data is not ready for ComDim.")
#   } else {
#     message("The data can be used for ComDim.")
#   }
#
#   ## Check y matrix (convert to dummy matrix if it is not already.)
#   if (method == "OPLS-DA") {
#     tmp <- unique(as.vector(y))
#     if (all(tmp %in% c(0, 1))) { # Is it dummy?
#
#       if (!is.matrix(y)) { # If it's a vector and the method is OPLS-R
#         y <- as.matrix(y)
#       }
#
#       classVect <- rep(NA, nrow(y))
#
#       if (ncol(y) == 1) {
#         classVect <- as.vector(y)
#       } else {
#         if (any(rowSums(y) != 1)) {
#           if (any(rowSums(y) == 0)) {
#             stop("At least one sample was not assigned to a class.")
#           } else if (any(colSums(y) > 1)) {
#             stop("At least one sample was assigned to more than one class.")
#           }
#         }
#       }
#       if (is.null(colnames(y))) {
#         colnames(y) <- paste0("Class", seq_len(ncol(y)))
#       }
#       for (i in seq_len(ncol(y))) {
#         classVect[which(y[, i] == 1)] <- colnames(y)[i]
#       }
#     } else { # If not dummy
#       if ((is.matrix(y) || is.data.frame(y)) && ncol(y) > 1) {
#         stop("ComDim-OPLS can only be applied if Y is a dummy matrix or a class vector")
#       }
#       classVect <- as.vector(y)
#       classVect_sorted <- sort(unique(classVect))
#
#       # Now generate the dummy matrix from y
#       y <- matrix(rep(0, length(classVect_sorted) * length(classVect)),
#         ncol = length(classVect_sorted), nrow = length(classVect)
#       )
#       for (i in seq_along(classVect_sorted)) {
#         y[which(classVect == classVect_sorted[i]), i] <- 1
#       }
#       colnames(y) <- classVect_sorted
#       rm(classVect_sorted)
#     }
#   }
#
#   if (!is.matrix(y)) { # If it's a vector and the method is OPLS-R
#     y <- as.matrix(y)
#   }
#
#   pieceBar <- 4 + 2 * ntable + ndim + nort # Number of updates in the progress bar.
#   pieceBar <- 80 / pieceBar
#   total_progress <- pieceBar
#
#   DimLabels <- paste0("CC", 1:ndim) # One label per predictive component.
#   OrtLabels <- paste0("ort", 1:nort) # One label per orthogonal component.
#   TableLabels <- blockNames(MB) # One label per block.
#
#   end_ini_time <- Sys.time() # To end the count of the calculation time.
#
#   if (loquace) {
#     message(sprintf("Initialisation finished after : %s millisecs", (end_ini_time - start_ini_time) * 1000))
#   }
#
#   utils::setTxtProgressBar(progress_bar, value = total_progress)
#
#   # NORMALISATION
#
#   X_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))
#   Xnorm_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))
#
#   res_calib <- list()
#   temp_tabCalib <- list()
#   s_n <- list()
#   res_calib$SingVal <- as.vector(NULL)
#   res_calib$NormMB <- as.vector(NULL)
#   res_calib$MeanMB <- list()
#
#   for (i in 1:ntable) {
#     res_calib$MeanMB[[TableLabels[i]]] <- colMeans(MB@Data[[i]])
#     names(res_calib$MeanMB[[TableLabels[i]]]) <- MB@Variables[[i]]
#
#     if (normalise) {
#       # Normalise original blocks
#       X_mean <- MB@Data[[i]] - matrix(data = rep(1, nrowMB), ncol = 1, nrow = nrowMB) %*% res_calib$MeanMB[[i]]
#       XX <- X_mean * X_mean
#       Norm_X <- sqrt(sum(XX))
#       X_Normed <- X_mean / Norm_X
#       res_calib$NormMB[i] <- Norm_X
#       temp_tabCalib[[i]] <- X_Normed
#       s_n[[i]] <- X_Normed
#     } else {
#       res_calib$NormMB[[TableLabels[i]]] <- rep(1, length(MB@Variables[[i]]))
#       names(res_calib$NormMB[[TableLabels[i]]]) <- MB@Variables[[i]]
#       temp_tabCalib[[i]] <- MB@Data[[i]]
#       s_n[[i]] <- MB@Data[[i]]
#     }
#
#     if (i == 1) {
#       X_mat[, 1:variable_number[1]] <- MB@Data[[i]]
#       Xnorm_mat[, 1:variable_number[1]] <- temp_tabCalib[[i]]
#     } else {
#       beg <- sum(variable_number[1:(i - 1)]) + 1
#       ending <- sum(variable_number[1:i])
#       X_mat[, beg:ending] <- MB@Data[[i]]
#       Xnorm_mat[, beg:ending] <- temp_tabCalib[[i]]
#     }
#
#     # PLS addition
#     meanX <- colMeans(Xnorm_mat)
#     meanY <- colMeans(y)
#     # End PLS addition
#
#     norm_comdim <- Sys.time()
#     if (loquace) {
#       message(sprintf("Normalization of block %s finished after : %s millisecs", i, (norm_comdim - start_ini_time) * 1000))
#     }
#     total_progress <- total_progress + pieceBar
#     utils::setTxtProgressBar(progress_bar, value = total_progress)
#   }
#
#   names(res_calib$NormMB) <- TableLabels
#
#   nR <- nrow(Xnorm_mat)
#   nC <- ncol(Xnorm_mat)
#
#   # =========================================================================
#   # KOPLS KERNEL-BASED EXTRACTION
#   # =========================================================================
#
#   end_comp_time <- Sys.time()
#   if (loquace) {
#     message(sprintf(
#       "Kernel preparation finished after : %s millisecs",
#       (end_comp_time - start_ini_time) * 1000
#     ))
#   }
#   total_progress <- total_progress + pieceBar * ntable
#   utils::setTxtProgressBar(progress_bar, value = total_progress)
#
#   # Helper: modified RV coefficient (zeros diagonal before computing)
#   .RVmod <- function(X, Y) {
#     AA <- tcrossprod(X)
#     diag(AA) <- 0
#     BB <- tcrossprod(Y)
#     diag(BB) <- 0
#     denom <- sqrt(sum(AA^2)) * sqrt(sum(BB^2))
#     if (denom < .Machine$double.eps) {
#       return(0)
#     }
#     sum(diag(crossprod(AA, BB))) / denom
#   }
#
#   K_blocks <- vector("list", ntable)
#   frob_norms <- numeric(ntable) # Frobenius norm of each raw block kernel
#   rv_weights <- numeric(ntable) # RV weight for each block
#
#   for (i in seq_len(ntable)) {
#     Xi <- temp_tabCalib[[i]]
#     if (normalise) {
#       ss <- sum(Xi^2)
#       if (ss > 1e-14) Xi <- Xi / sqrt(ss)
#     }
#     Ki <- tcrossprod(Xi)
#     frob_Ki <- sqrt(sum(Ki^2))
#     frob_norms[i] <- if (frob_Ki > 1e-14) frob_Ki else 1
#     K_norm_i <- Ki / frob_norms[i]
#     rv_weights[i] <- (.RVmod(K_norm_i, y) + 1) / 2
#     K_blocks[[i]] <- K_norm_i
#   }
#   names(rv_weights) <- TableLabels
#
#   # RV-weighted combined kernel
#   W_mat <- Reduce("+", lapply(
#     seq_len(ntable),
#     function(i) rv_weights[i] * K_blocks[[i]]
#   ))
#
#   # Double-centre: K_combined = H * W_mat * H  (H = I - 1/n * J)
#   K_combined <- scale(t(scale(W_mat, scale = FALSE)), scale = FALSE)
#
#   rm(s_n)
#
#   # --- 2. Centre Y; SVD of Y' K Y -------------------------------------------
#   Y_c <- sweep(y, 2, colMeans(y)) # n x ncol(y)
#   YKY <- crossprod(Y_c, K_combined %*% Y_c) # ncol(y) x ncol(y)
#
#   A <- ndim
#   CSV <- svd(YKY, nu = A, nv = A)
#   Cp <- CSV$u[, 1:A, drop = FALSE] # ncol(y) x ndim  Y loadings
#   Sp <- diag(CSV$d[1:A]^(-0.5), nrow = A) # ndim x ndim     scale
#   Up <- Y_c %*% Cp # n x ndim        Y scores
#
#   colnames(Cp) <- DimLabels
#   rownames(Cp) <- colnames(y)
#   colnames(Up) <- DimLabels
#   rownames(Up) <- MB@Samples
#
#   # --- 3. Orthogonal extraction (KOPLS loop) --------------------------------
#   Q_ort <- matrix(NA_real_, nrow = nrowMB, ncol = nort)
#   saliences_ort <- matrix(NA_real_, nrow = ntable, ncol = nort)
#
#   to_list <- vector("list", nort) # unit-norm ort scores
#   co_list <- vector("list", nort) # ort direction in Y-score space
#   so_list <- vector("list", nort) # singular value
#   toNorm_list <- vector("list", nort) # norm before unit-normalisation
#
#   # Per-component storage for R2X and R2Y
#   # Output will have ndim + nort values: CC1..CC{ndim}, ort1..ort{nort}
#   sstot_K      <- sum(diag(K_combined))    # trace(K_original) = SS(X)
#   Kdiag_levels <- numeric(nort + 1)        # trace(K_both) before/after each ort deflation
#   Kdiag_levels[1] <- sstot_K              # before any ort deflation
#   Tp_levels    <- vector("list", nort + 1) # predictive scores at each ort level
#
#   K_right <- K_combined # n x n  right-only deflated kernel
#   K_both <- K_combined # n x n  both-sides deflated kernel
#
#   for (jj in seq_len(nort)) {
#     Tp_curr <- crossprod(K_right, Up %*% Sp) # = t(K_right) %*% Up %*% Sp
#     Tp_levels[[jj]] <- Tp_curr              # Tp using K deflated by (jj-1) ort components
#     K_resid <- K_both - tcrossprod(Tp_curr) # K - Tp Tp'
#
#     sv_tmp <- svd(crossprod(Tp_curr, K_resid %*% Tp_curr), nu = 1, nv = 1)
#     co_j <- sv_tmp$u[, 1, drop = FALSE] # ndim x 1
#     so_j <- sv_tmp$d[1]
#
#     to_raw <- K_resid %*% Tp_curr %*% co_j /
#       sqrt(max(so_j, .Machine$double.eps)) # n x 1
#     toNorm_j <- sqrt(sum(to_raw^2))
#     if (toNorm_j < 1e-14) toNorm_j <- 1
#     to_j <- to_raw / toNorm_j # unit-norm
#
#     to_list[[jj]] <- to_j
#     co_list[[jj]] <- co_j
#     so_list[[jj]] <- so_j
#     toNorm_list[[jj]] <- toNorm_j
#     Q_ort[, jj] <- as.vector(to_j)
#
#     for (i in seq_len(ntable)) {
#       saliences_ort[i, jj] <-
#         as.vector(crossprod(to_j, K_blocks[[i]]) %*% to_j)
#     }
#
#     # Kernel deflation
#     I_to <- diag(nrowMB) - tcrossprod(to_j)
#     K_right <- K_right %*% I_to # right only
#     K_both <- I_to %*% K_both %*% I_to # both sides
#     Kdiag_levels[jj + 1] <- sum(diag(K_both)) # trace after jj ort deflations
#
#     ort_comdim <- Sys.time()
#     if (loquace) {
#       message(sprintf(
#         "Orthogonal component %s determined after : %s millisecs",
#         jj, (ort_comdim - start_ini_time) * 1000
#       ))
#     }
#     total_progress <- total_progress + pieceBar
#     utils::setTxtProgressBar(progress_bar, value = total_progress)
#   }
#
#   colnames(Q_ort) <- OrtLabels
#   rownames(Q_ort) <- MB@Samples
#   colnames(saliences_ort) <- OrtLabels
#   rownames(saliences_ort) <- TableLabels
#
#   # --- 4. Final predictive scores -------------------------------------------
#   Tp <- crossprod(K_right, Up %*% Sp) # n x ndim
#   Tp_levels[[nort + 1]] <- Tp         # final predictive scores (all ort deflations applied)
#
#   # Inner regression coefficient
#   Bt <- pracma::pinv(crossprod(Tp)) %*% crossprod(Tp, Up) # ndim x ndim
#
#   # Predictive block saliences
#   saliences <- matrix(NA_real_, nrow = ntable, ncol = ndim)
#   for (i in seq_len(ntable)) {
#     for (j in seq_len(ndim)) {
#       tp_j <- Tp[, j]
#       saliences[i, j] <-
#         as.vector(crossprod(tp_j, K_blocks[[i]]) %*% tp_j)
#     }
#   }
#   colnames(saliences) <- DimLabels
#   rownames(saliences) <- TableLabels
#
#   # SingVal: SS of predictive scores
#   res_calib$SingVal <- setNames(colSums(Tp^2), DimLabels)
#
#   # R2X per component (non-cumulative)
#   r2x_pred        <- setNames(colSums(Tp^2) / sstot_K, DimLabels)
#   r2x_ort         <- setNames(-diff(Kdiag_levels) / sstot_K, OrtLabels)
#   res_calib$explained <- c(r2x_pred, r2x_ort)
#
#   # Store predictive scores as Q
#   Q <- Tp
#   colnames(Q) <- DimLabels
#   rownames(Q) <- MB@Samples
#   res_calib$Q <- Q
#
#   ## Adding metadata
#   res_calib$metadata <- list()
#   if (length(MB@Metadata) != 0) {
#     res_calib$metadata[[1]] <- MB@Metadata
#   }
#
#   iter_comdim <- Sys.time()
#   if (loquace) {
#     message(sprintf(
#       "KOPLS scores finished after : %s millisecs",
#       (iter_comdim - start_ini_time) * 1000
#     ))
#   }
#   total_progress <- total_progress + pieceBar * (1 + ndim)
#   utils::setTxtProgressBar(progress_bar, value = min(total_progress, 79))
#
#   rm(K_combined, K_right, K_both)
#
#   # --- 5. R2Y per component --------------------------------------------------
#   sstot_Y_r2 <- sum(Y_c^2)
#
#   Tp0 <- Tp_levels[[1]]  # pred scores before any ort deflation
#   r2y_pred_cum <- numeric(ndim)
#   for (.j in seq_len(ndim)) {
#     Tp_j   <- Tp0[, seq_len(.j), drop = FALSE]
#     Bt_j   <- pracma::pinv(crossprod(Tp_j)) %*% crossprod(Tp_j, Up[, seq_len(.j), drop = FALSE])
#     Yhat_j <- Tp_j %*% Bt_j %*% t(Cp[, seq_len(.j), drop = FALSE])
#     r2y_pred_cum[.j] <- if (sstot_Y_r2 < .Machine$double.eps) 0 else
#       1 - sum((Yhat_j - Y_c)^2) / sstot_Y_r2
#   }
#   r2y_pred <- setNames(c(r2y_pred_cum[1], diff(r2y_pred_cum)), DimLabels)
#
#   r2y_lev <- numeric(nort + 1)
#   for (.i in seq_len(nort + 1)) {
#     Tp_i      <- Tp_levels[[.i]]
#     Bt_i      <- pracma::pinv(crossprod(Tp_i)) %*% crossprod(Tp_i, Up)
#     Yhat_i    <- Tp_i %*% Bt_i %*% t(Cp)
#     r2y_lev[.i] <- if (sstot_Y_r2 < .Machine$double.eps) 0 else
#       1 - sum((Yhat_i - Y_c)^2) / sstot_Y_r2
#   }
#   r2y_ort <- setNames(diff(r2y_lev), OrtLabels)
#   r2y <- c(r2y_pred, r2y_ort)
#
#   ## Calculate Sums of saliences for each Dim
#   Sum_saliences_Dim <- colSums(saliences)
#   names(Sum_saliences_Dim) <- DimLabels
#   res_calib$Sum_saliences_Dim <- Sum_saliences_Dim
#   rm(Sum_saliences_Dim)
#
#   # Calculate Sums of saliences for each Table
#   Sum_saliences_Tab <- rowSums(saliences)
#   names(Sum_saliences_Tab) <- TableLabels
#   res_calib$Sum_saliences_Tab <- Sum_saliences_Tab
#   rm(Sum_saliences_Tab)
#
#   res_calib$saliences <- saliences
#   colnames(res_calib$saliences) <- DimLabels
#   rownames(res_calib$saliences) <- TableLabels
#
#   ## LOADING COMPUTATION
#
#   L_X <- list()
#   T_Loc_ort <- list()
#   T_Loc <- list()
#   for (i in 1:ntable) { # Prepare lists
#     T_Loc_ort[[TableLabels[i]]] <- matrix(, ncol = nort, nrow = nrowMB)
#     T_Loc[[TableLabels[i]]] <- matrix(, ncol = ndim, nrow = nrowMB)
#   }
#   b <- matrix(, ncol = ndim, nrow = ntable)
#
#   if (!requireNamespace("pracma", quietly = TRUE)) {
#     stop("The package pracma is needed.")
#   }
#
#   # Save pre-deflation blocks for VIP computation
#   X_orig_blocks <- temp_tabCalib
#
#   # Orthogonal loadings
#   orthoP <- NULL
#   for (jj in 1:nort) {
#     for (i in 1:ntable) {
#       tempo <- t(temp_tabCalib[[i]]) %*% Q_ort[, jj] # Scaled CD 'local' Loadings
#       T_Loc_ort[[TableLabels[i]]][, jj] <- temp_tabCalib[[i]] %*% (tempo %*% pracma::pinv(t(tempo) %*% tempo)) # local Scores
#       # Deflate each temp_tabCalib by orthogonal component
#       temp_tabCalib[[i]] <- temp_tabCalib[[i]] - Q_ort[, jj] %*% t(tempo)
#       if (i == 1) {
#         tempo2 <- tempo
#       } else {
#         tempo2 <- append(tempo2, tempo)
#       }
#     }
#     if (jj == 1) {
#       orthoP <- tempo2
#     } else {
#       orthoP <- cbind(orthoP, tempo2)
#     }
#   }
#   orthoP <- as.matrix(orthoP)
#   colnames(orthoP) <- OrtLabels
#
#   for (i in 1:ntable) {
#     rownames(T_Loc_ort[[TableLabels[i]]]) <- MB@Samples
#     colnames(T_Loc_ort[[TableLabels[i]]]) <- OrtLabels
#   }
#
#   # Rebuild Xnorm_mat from orthogonally-deflated temp_tabCalib.
#   for (i in 1:ntable) {
#     if (i == 1) {
#       Xnorm_mat <- temp_tabCalib[[i]]
#     } else {
#       Xnorm_mat <- cbind(Xnorm_mat, temp_tabCalib[[i]])
#     }
#   }
#
#   # Predictive loadings
#   for (j in 1:ndim) {
#     T_mat <- matrix(, nrow = nrowMB, ncol = 0)
#
#     q_j <- res_calib$Q[, j]
#     q_j_ss <- sum(q_j^2)
#     q_j_unit <- if (q_j_ss > 1e-14) q_j / sqrt(q_j_ss) else q_j
#
#     for (i in 1:ntable) {
#       temp <- t(temp_tabCalib[[i]]) %*% q_j_unit # Scaled CD 'local' Loadings
#
#       T_mat <- cbind(T_mat, temp_tabCalib[[i]] %*% (temp %*% pracma::pinv(t(temp) %*% temp))) # local Scores
#
#       T_Loc[[TableLabels[i]]][, j] <- temp_tabCalib[[i]] %*% (temp %*% pracma::pinv(t(temp) %*% temp)) # local Scores
#
#       # Deflate each temp_tabCalib
#       temp_tabCalib[[i]] <- temp_tabCalib[[i]] - q_j_unit %*% t(temp)
#     }
#
#     # For each CC
#     # MLR b-coefficients between Local and Global Scores
#     b[, j] <- pracma::pinv(t(T_mat) %*% T_mat) %*% t(T_mat) %*% res_calib$Q[, j]
#   }
#
#   for (i in 1:ntable) {
#     rownames(T_Loc[[TableLabels[i]]]) <- rownames(res_calib$Q)
#     colnames(T_Loc[[TableLabels[i]]]) <- colnames(res_calib$Q)
#   }
#
#   # Calculate Global Loadings from orthogonally-deflated X
#   L_CD_Vec <- t(Xnorm_mat) %*% res_calib$Q # Scaled CD 'global' Loadings
#
#   load_comdim <- Sys.time()
#   if (loquace) {
#     message(sprintf("Loadings finished after : %s millisecs",
#                     (load_comdim - start_ini_time) * 1000))
#   }
#
#   total_progress <- total_progress + pieceBar
#   utils::setTxtProgressBar(progress_bar, value = total_progress)
#
#   rm(X_mat)
#   rm(T_mat)
#   rm(temp_tabCalib)
#   rm(temp)
#
#   ## Output
#   end_function <- Sys.time()
#
#   ## KOPLS B matrix (variable-space regression coefficient)
#   UpSp <- Up %*% Sp # n x ndim
#   W_pred <- matrix(NA_real_, nrow = sum(variable_number), ncol = ndim)
#   k_idx <- 1
#   for (i in seq_len(ntable)) {
#     vi <- variable_number[i]
#     Xb_ort <- Xnorm_mat[, k_idx:(k_idx + vi - 1), drop = FALSE]
#     W_pred[k_idx:(k_idx + vi - 1), ] <- t(Xb_ort) %*% UpSp * rv_weights[i] / frob_norms[i]
#     k_idx <- k_idx + vi
#   }
#   B_kopls <- W_pred %*% Bt %*% t(Cp) # p x ncol(y)
#   B0_kopls <- colMeans(y) - colMeans(Xnorm_mat) %*% B_kopls # 1 x ncol(y)
#   colnames(W_pred) <- DimLabels
#   rm(UpSp)
#
#   res_calib$b <- b
#   colnames(res_calib$b) <- DimLabels
#   rownames(res_calib$b) <- TableLabels
#
#   res_calib$T_Loc <- T_Loc
#   rm(T_Loc)
#
#   varnames <- vector()
#   for (i in 1:ntable) {
#     varnames <- append(varnames, MB@Variables[[i]])
#   }
#
#   res_calib$P <- L_CD_Vec
#   colnames(res_calib$P) <- DimLabels
#   rownames(res_calib$P) <- varnames
#   rownames(orthoP) <- varnames
#
#   rm(L_CD_Vec)
#
#   ## VIP
#   .vip_for_scores <- function(scores_mat, Y_raw, X_b, p_b) {
#     ncomp <- ncol(scores_mat)
#     ss_T <- colSums(scores_mat^2)
#     if (any(ss_T < .Machine$double.eps)) {
#       return(rep(NA_real_, p_b))
#     }
#     Qs <- crossprod(Y_raw, scores_mat) %*% diag(1 / ss_T, nrow = ncomp)
#     ss_Q <- colSums(Qs^2)
#     if (any(ss_Q < .Machine$double.eps)) {
#       return(rep(NA_real_, p_b))
#     }
#     Us <- Y_raw %*% Qs %*% diag(1 / ss_Q, nrow = ncomp)
#     ss_U <- colSums(Us^2)
#     if (any(ss_U < .Machine$double.eps)) {
#       return(rep(NA_real_, p_b))
#     }
#     Ws <- crossprod(X_b, Us) %*% diag(1 / ss_U, nrow = ncomp)
#     Ws <- apply(Ws, 2, function(w) {
#       nrm <- sqrt(sum(w^2))
#       if (nrm > 1e-14) w / nrm else w
#     })
#     s_a <- ss_T * ss_Q
#     sum_s <- sum(s_a)
#     if (sum_s < .Machine$double.eps) {
#       return(rep(NA_real_, p_b))
#     }
#     as.vector(sqrt(p_b * (Ws^2 %*% s_a) / sum_s))
#   }
#
#   .vip_ort_for_loadings <- function(ort_scores_mat, ort_loadings_b) {
#     ncomp <- ncol(ort_scores_mat)
#     p_b   <- nrow(ort_loadings_b)
#     ss_T  <- colSums(ort_scores_mat^2)
#     if (any(ss_T < .Machine$double.eps)) return(rep(NA_real_, p_b))
#     P_norm <- apply(ort_loadings_b, 2, function(pp) {
#       nrm <- sqrt(sum(pp^2)); if (nrm > 1e-14) pp / nrm else pp
#     })
#     if (is.null(dim(P_norm))) P_norm <- matrix(P_norm, ncol = 1)
#     sum_s <- sum(ss_T)
#     if (sum_s < .Machine$double.eps) return(rep(NA_real_, p_b))
#     as.vector(sqrt(p_b * (P_norm^2 %*% ss_T) / sum_s))
#   }
#
#   VIP_blocks <- setNames(vector("list", ntable), TableLabels)
#   k_ort <- 1L  # running row index into orthoP for the current block
#   for (i in seq_len(ntable)) {
#     X_b   <- X_orig_blocks[[i]]
#     p_b   <- variable_number[i]
#     vip_p <- .vip_for_scores(Tp, y, X_b, p_b)
#     ort_P_b <- orthoP[k_ort:(k_ort + p_b - 1L), , drop = FALSE]
#     vip_o <- .vip_ort_for_loadings(Q_ort, ort_P_b)
#     vip_t <- sqrt((vip_p^2 + vip_o^2) / 2)
#     vdf <- data.frame(
#       p = vip_p, o = vip_o, tot = vip_t,
#       row.names = MB@Variables[[i]]
#     )
#     VIP_blocks[[TableLabels[i]]] <- vdf
#     k_ort <- k_ort + p_b
#   }
#
#   VIP_global <- as.vector(unlist(lapply(VIP_blocks, function(b) b$tot)))
#   names(VIP_global) <- varnames
#
#   # SingVal and R2X are computed from scores (set above).
#
#   rm(saliences)
#
#   # Define block length
#   belong_block <- rep(0, nrow(res_calib$P))
#   k <- 1
#   for (i in 1:ntable) {
#     belong_block[k:(k + variable_number[i] - 1)] <- TableLabels[i]
#     k <- k + variable_number[i]
#   }
#
#   rm(Q)
#
#   ## Y prediction using the orthogonally-deflated X and KOPLS B matrix
#   Ypred <- Xnorm_mat %*% B_kopls +
#     matrix(as.vector(B0_kopls), nrow = nrowMB, ncol = ncol(y), byrow = TRUE)
#   rownames(Ypred) <- MB@Samples
#   if (!is.null(colnames(y))) colnames(Ypred) <- colnames(y)
#
#   ## Predict classes
#   if (method == "OPLS-DA") {
#     predClass <- rep(NA, nrow(y))
#
#     if (decisionRule == "max") {
#       for (i in seq_len(nrow(Ypred))) {
#         xx <- which(Ypred[i, ] == max(Ypred[i, ], na.rm = TRUE))
#         if (length(xx) == 1) {
#           predClass[i] <- colnames(y)[xx] # Y stores the Class names
#         } else {
#           predClass[i] <- NaN # There is a tie.
#           warning(sprintf("Class for sample %s could not be predicted.", i))
#         }
#       }
#     } else if (decisionRule == "fixed") {
#       for (i in seq_len(nrow(Ypred))) {
#         xx <- which(Ypred[, i] > 1 / ncol(y))
#         if (length(xx) == 1) {
#           predClass[i] <- colnames(y)[xx] # Y stores the Class names
#         } else {
#           predClass[i] <- NaN # There is a tie.
#           warning(sprintf("Class for sample %s could not be predicted.", i))
#         }
#       }
#     }
#
#     if (is.matrix(y)) {
#       meanY <- colMeans(y)
#     } else {
#       meanY <- mean(y)
#     }
#
#     classy <- colnames(y)
#     Q2 <- DQ2 <- specvec <- sensvec <- rep(NA, length(classy))
#     confusionMatrix <- list()
#
#     for (i in seq_along(classy)) {
#       pos <- which(classVect == classy[i])
#       neg <- which(classVect != classy[i])
#       tp <- length(which(classVect == classy[i] & predClass == classy[i]))
#       tn <- length(which(classVect != classy[i] & predClass != classy[i]))
#       fp <- length(which(classVect != classy[i] & predClass == classy[i]))
#       fn <- length(which(classVect == classy[i] & predClass != classy[i]))
#
#       if (length(tp) == 0) {
#         tp <- 0
#       }
#
#       if (length(tn) == 0) {
#         tn <- 0
#       }
#
#       if (length(fp) == 0) {
#         fp <- 0
#       }
#
#       if (length(fn) == 0) {
#         fn <- 0
#       }
#
#       sensvec[i] <- tp / length(pos)
#       specvec[i] <- tn / length(neg)
#
#       cm <- matrix(c(tp, fp, fn, tn), ncol = 2, nrow = 2)
#       rownames(cm) <- c("predClass1", "predClass0")
#       colnames(cm) <- c("trueClass1", "trueClass0")
#       confusionMatrix[[classy[i]]] <- cm
#
#       ## DQ2
#       E0 <- Ypred[neg, i] - y[neg, i] # Calculate Residuals of Class0 samples
#       E1 <- Ypred[pos, i] - y[pos, i] # Calculate Residuals of Class1 samples
#
#       E0count <- which(E0 > 0) # Find predictions for Class0 samples larger than 0
#       E1count <- which(E1 < 0) # Find predictions for Class1 samples larger than 1
#
#       SSE0_all <- t(E0) %*% E0
#       SSE1_all <- t(E1) %*% E1
#
#       SSE0 <- t(E0[E0count]) %*% E0[E0count]
#       SSE1 <- t(E1[E1count]) %*% E1[E1count]
#
#       PRESSD <- SSE0 + SSE1
#       PRESS <- SSE0_all + SSE1_all
#
#       Ym <- y[, i] - mean(y[, i])
#       TSS <- t(Ym) %*% Ym
#
#       DQ2[i] <- 1 - (PRESSD / TSS)
#       Q2[i] <- 1 - PRESS / TSS
#     }
#     names(Q2)      <- classy
#     names(DQ2)     <- classy
#     names(sensvec) <- classy
#     names(specvec) <- classy
#   } else if (method == "OPLS-R") {
#     PRESS <- t(Ypred - y) %*% (Ypred - y)
#     Ym <- y - mean(y)
#     TSS <- t(Ym) %*% Ym
#     Q2 <- 1 - PRESS / TSS
#   }
#
#
#   # --- Internal kernel centering helpers ---
#   .cv_center_KtrTr <- function(K) scale(t(scale(K, scale = FALSE)), scale = FALSE)
#
#   .cv_center_KteTr <- function(KteTr, KtrTr) {
#     ntr <- nrow(KtrTr)
#     nte <- nrow(KteTr)
#     sc <- matrix(1 / ntr, nte, ntr)
#     (KteTr - sc %*% KtrTr) %*% (diag(ntr) - (1 / ntr) * matrix(1, ntr, ntr))
#   }
#
#   # --- Internal KOPLS model fit on a centred kernel ---
#   .cv_kopls_fit <- function(K_cc, Yc, a, nox) {
#     n <- nrow(K_cc)
#     CSV <- svd(crossprod(Yc, K_cc %*% Yc), nu = a, nv = a)
#     Cp <- CSV$u[, seq_len(a), drop = FALSE]
#     Sp <- diag(CSV$d[seq_len(a)]^(-0.5), nrow = a)
#     Up <- Yc %*% Cp
#     UpSp <- Up %*% Sp
#
#     K_right <- K_cc
#     K_both <- K_cc
#     to_l <- co_l <- so_l <- toNorm_l <- Tp_tr_l <- list()
#     K_right_l <- list()
#     K_right_l[[1]] <- K_cc
#
#     for (jj in seq_len(nox)) {
#       Tp_j <- crossprod(K_right, UpSp)
#       K_resid <- K_both - tcrossprod(Tp_j)
#       sv <- svd(crossprod(Tp_j, K_resid %*% Tp_j), nu = 1, nv = 1)
#       co_j <- sv$u[, 1, drop = FALSE]
#       so_j <- sv$d[1]
#       to_raw <- K_resid %*% Tp_j %*% co_j / sqrt(max(so_j, .Machine$double.eps))
#       tn <- sqrt(sum(to_raw^2))
#       if (tn < 1e-14) tn <- 1
#       to_j <- to_raw / tn
#
#       Tp_tr_l[[jj]] <- Tp_j
#       to_l[[jj]] <- to_j
#       co_l[[jj]] <- co_j
#       so_l[[jj]] <- so_j
#       toNorm_l[[jj]] <- tn
#
#       I_to <- diag(n) - tcrossprod(to_j)
#       K_right <- K_right %*% I_to
#       K_both <- I_to %*% K_both %*% I_to
#       K_right_l[[jj + 1]] <- K_right
#     }
#
#     Tp_f <- crossprod(K_right, UpSp)
#     Bt_f <- pracma::pinv(crossprod(Tp_f)) %*% crossprod(Tp_f, Up)
#
#     list(
#       Cp = Cp, Sp = Sp, Up = Up, UpSp = UpSp, Bt = Bt_f,
#       to = to_l, co = co_l, so = so_l, toNorm = toNorm_l,
#       Tp_tr = Tp_tr_l, K_right = K_right_l
#     )
#   }
#
#   # --- Internal KOPLS predict: returns Yhat in centred Y space ---
#   .cv_kopls_predict <- function(KteTr_c, m, nox) {
#     KteTr_curr <- KteTr_c
#
#     for (i in seq_len(nox)) {
#       Tp_te_i <- KteTr_curr %*% m$UpSp
#       K_res_tetr <- KteTr_curr - tcrossprod(Tp_te_i, m$Tp_tr[[i]])
#       to_te_i <- (K_res_tetr %*% m$Tp_tr[[i]] %*% m$co[[i]]) /
#         sqrt(max(m$so[[i]], .Machine$double.eps))
#       to_te_i <- to_te_i / m$toNorm[[i]]
#       KteTr_curr <- KteTr_curr -
#         tcrossprod(to_te_i, m$K_right[[i]] %*% m$to[[i]])
#     }
#
#     Tp_te <- KteTr_curr %*% m$UpSp
#     Tp_te %*% m$Bt %*% t(m$Cp)
#   }
#
#   if (cv.k >= 2) {
#     n_cv   <- nrowMB
#     n_lev  <- nort + 1L   # levels: nox = 0, 1, ..., nort
#
#     Ypred_cv_all <- vector("list", n_lev)
#     for (.lev in seq_len(n_lev))
#       Ypred_cv_all[[.lev]] <- matrix(NA_real_, nrow = n_cv, ncol = ncol(y),
#                                      dimnames = list(NULL, colnames(y)))
#     pressy_cls_lev <- matrix(0, nrow = n_lev, ncol = ncol(y))
#     TSS_cv_cls     <- numeric(ncol(y))
#     cv_test_order  <- integer(0)
#
#     for (fold in seq_len(cv.k)) {
#       pred_idx  <- seq.int(fold, n_cv, by = cv.k)
#       train_idx <- setdiff(seq_len(n_cv), pred_idx)
#
#       KtrTr_raw <- W_mat[train_idx, train_idx, drop = FALSE]
#       KteTr_raw <- W_mat[pred_idx,  train_idx, drop = FALSE]
#
#       KtrTr_c <- .cv_center_KtrTr(KtrTr_raw)
#       KteTr_c <- .cv_center_KteTr(KteTr_raw, KtrTr_raw)
#
#       Y_tr_mean <- colMeans(y[train_idx, , drop = FALSE])
#       Y_train_c <- sweep(y[train_idx, , drop = FALSE], 2, Y_tr_mean)
#       Y_test_c  <- sweep(y[pred_idx,  , drop = FALSE], 2, Y_tr_mean)
#
#       for (nox in 0L:nort) {
#         .lev  <- nox + 1L
#         m_lev <- .cv_kopls_fit(KtrTr_c, Y_train_c, ndim, nox)
#         Yhat  <- .cv_kopls_predict(KteTr_c, m_lev, nox)
#         Ypred_cv_all[[.lev]][pred_idx, ] <- sweep(Yhat, 2, Y_tr_mean, "+")
#         pressy_cls_lev[.lev, ] <- pressy_cls_lev[.lev, ] + colSums((Y_test_c - Yhat)^2)
#       }
#
#       TSS_cv_cls    <- TSS_cv_cls + colSums(Y_test_c^2)
#       cv_test_order <- c(cv_test_order, pred_idx)
#     }
#
#     Q2_cv <- setNames(1 - pressy_cls_lev[n_lev, ] / TSS_cv_cls, colnames(y))
#     q2_names  <- c("p", paste0("po", seq_len(nort)))
#     Q2_levels <- setNames(1 - rowSums(pressy_cls_lev) / sum(TSS_cv_cls), q2_names)
#
#     if (method == "OPLS-DA") {
#       DQ2_cv      <- setNames(rep(NA_real_, ncol(y)), colnames(y))
#       Ypred_final <- Ypred_cv_all[[n_lev]]
#       for (i in seq_len(ncol(y))) {
#         y_obs  <- y[cv_test_order, i]
#         y_pred <- Ypred_final[cv_test_order, i]
#         E0     <- y_pred[y_obs == 0]
#         E1     <- y_pred[y_obs == 1] - 1
#         PRESSD <- sum(E0[E0 > 0]^2) + sum(E1[E1 < 0]^2)
#         TSS_i  <- sum((y_obs - mean(y_obs))^2)
#         DQ2_cv[i] <- 1 - PRESSD / TSS_i
#       }
#     } else {
#       DQ2_cv <- vector()
#     }
#
#     Q2 <- Q2_cv
#     if (method == "OPLS-DA") DQ2 <- DQ2_cv
#
#     Ypred_cv   <- Ypred_cv_all[[n_lev]]
#     Ypred_cv_p <- Ypred_cv_all[[1L]]
#
#     cv_result <- list(
#       k          = cv.k,
#       Ypred      = Ypred_cv,
#       Ypred.p    = Ypred_cv_p,
#       Q2         = Q2_cv,
#       Q2.levels  = Q2_levels,
#       DQ2        = if (method == "OPLS-DA") DQ2_cv else NULL
#     )
#   } else {
#     if (cv.k != 0) {
#       warning("'cv.k' must be >= 2 or 0 (skip CV). No cross-validation performed.")
#     }
#     cv_result <- list()
#   }
#
#   end_output <- Sys.time()
#   running_time <- (end_output - start_ini_time)
#   if (loquace) {
#     message(sprintf("Analysis finished after : %s millisecs", running_time * 1000))
#   }
#
#   progress_bar <- utils::txtProgressBar(min = 0, max = 80, style = 3, char = "=")
#   utils::setTxtProgressBar(progress_bar, value = 80)
#
#   close(progress_bar)
#
#   res_calib <- new("ComDim",
#     Method = method,
#     ndim = ndim,
#     Q.scores = res_calib$Q,
#     T.scores = res_calib$T_Loc,
#     P.loadings = res_calib$P,
#     Saliences = res_calib$saliences,
#     Orthogonal = list(
#       nort = nort,
#       Q.scores = Q_ort,
#       T.scores = T_Loc_ort,
#       P.loadings.ort = orthoP,
#       Saliences.ort = saliences_ort
#     ),
#     VIP = VIP_global,
#     VIP.block = VIP_blocks,
#     R2X = res_calib$explained,
#     R2Y = r2y,
#     Q2 = Q2,
#     DQ2 = if (method == "OPLS-DA") {
#       DQ2
#     } else {
#       vector()
#     },
#     Singular = res_calib$SingVal,
#     Mean = list(
#       MeanMB = res_calib$MeanMB,
#       MeanY = meanY
#     ),
#     Norm = list(
#       NormMB = res_calib$NormMB,
#       FrobNorms = frob_norms,
#       RVweights = rv_weights
#     ),
#     PLS.model = list(
#       W = W_pred,
#       B = B_kopls,
#       B0 = B0_kopls,
#       Y = y
#     ),
#     cv = cv_result,
#     Prediction = if (method == "OPLS-DA") {
#       list(
#         Y.pred = Ypred,
#         decisionRule = decisionRule,
#         trueClass = classVect,
#         predClass = data.frame(predClass = predClass, row.names = MB@Samples),
#         Sensitivity = sensvec,
#         Specificity = specvec,
#         confusionMatrix = confusionMatrix
#       )
#     } else {
#       list(Y.pred = Ypred)
#     },
#     Metadata = res_calib$metadata,
#     variable.block = belong_block,
#     runtime = as.numeric(running_time)
#   )
#
#   return(res_calib)
# }


#' OPLS_NIPALS_DNR
#'
#' One NIPALS OPLS step on a (lambda-weighted) concatenated block matrix.
#' Computes one predictive component and one orthogonal component in a single
#' pass. This function is the recommended building block for constructing an
#' OPLS-based \code{FUN} argument for \code{ComDim_y()} when \code{nort > 0}.
#' For \code{nort = 0} (plain PLS) a simpler PLS wrapper is sufficient.
#' @param W Numeric matrix (n x p): the concatenated, lambda-weighted blocks
#'   as passed by \code{ComDim_y()}.
#' @param y Numeric matrix (n x q): the response block (dummy matrix for
#'   discriminant analysis, numeric matrix for regression). Only the first
#'   column drives the NIPALS u-score iteration; all columns are used to
#'   compute the Y-loading \code{q}.
#' @param threshold Convergence threshold for the u-score update (default
#'   \code{1e-10}).
#' @return A named list:
#' \describe{
#'   \item{t_pred}{Predictive X-score (length n).}
#'   \item{w_pred}{Predictive X-weight (length p), L2-normalised.}
#'   \item{p}{X-loading (length p).}
#'   \item{q}{Y-loading (length q = ncol(y)).}
#'   \item{u}{Y-score (length n).}
#'   \item{t_ort}{Orthogonal X-score (length n).}
#'   \item{w_ort}{Orthogonal X-weight (length p), L2-normalised.}
#'   \item{p_ort}{Orthogonal X-loading (length p).}
#' }
#' @seealso \code{\link{ComDim_y}} for the multi-block OPLS wrapper that uses
#'   this function.
#' @export
OPLS_NIPALS_DNR <- function(W = W, y = y, threshold = 1e-10) {
  u <- y[, 1]
  u_old <- u + threshold * 10 # Ensures the while loop runs at least once
  iters <- 0

  while (norm(u - u_old, "2") > threshold && iters < 100) {
    iters <- iters + 1
    u_old <- u

    # Predictive X-weight and score
    w_pred <- t(W) %*% u
    w_pred_len <- norm(w_pred, "2")
    if (w_pred_len > 1e-12) w_pred <- w_pred / w_pred_len
    t_pred <- W %*% w_pred

    # Y-loading and score
    q_y <- t(y) %*% t_pred / as.vector(t(t_pred) %*% t_pred)
    u <- y %*% q_y / as.vector(t(q_y) %*% q_y)
  }

  # X-loading for predictive direction
  p <- t(W) %*% t_pred / as.vector(t(t_pred) %*% t_pred)

  # Orthogonal weight: component of p orthogonal to w_pred
  w_ort <- p - w_pred * as.vector(t(w_pred) %*% p) / as.vector(t(w_pred) %*% w_pred)
  w_ort_len <- norm(w_ort, "2")
  if (w_ort_len > 1e-12) w_ort <- w_ort / w_ort_len

  # Orthogonal score and loading
  t_ort <- W %*% w_ort
  p_ort <- t(W) %*% t_ort / as.vector(t(t_ort) %*% t_ort)

  return(list(
    t_pred = t_pred,
    w_pred = w_pred,
    p = p,
    q = q_y,
    u = u,
    t_ort = t_ort,
    w_ort = w_ort,
    p_ort = p_ort
  ))
}
