#' ComDim_OPLS
#'
#' Finding common dimensions in multi-block datasets using OPLS. Also known as ConsensusOPLS
#' (ComDim-OPLS) for multiblock structures: orthogonal components uncorrelated
#' with Y are extracted from all blocks simultaneously before the predictive
#' components are computed.
#' @param MB A MultiBlock object.
#' @param y The Y-block. A class vector or dummy matrix for OPLS-DA, or a
#'   numeric matrix/vector for OPLS-R.
#' @param ndim Number of predictive Common Dimensions. Default is 1.
#' @param nort Number of orthogonal Common Dimensions. Default is 1.
#' @param method 'OPLS-DA' for discriminant analysis or 'OPLS-R' for regression.
#' @param decisionRule Only used if method is 'OPLS-DA'. If 'fixed', samples are
#'   assigned to the class with Y-hat above 1/nclasses. If 'max', samples are
#'   assigned to the class with the highest Y-hat.
#' @param normalise To apply block normalisation. FALSE == no (default), TRUE == yes.
#' @param threshold The threshold limit to stop the iterations (default 1e-10).
#' @param loquace To display the calculation times. TRUE == yes, FALSE == no (default).
#' @param cv.k Number of folds for k-fold cross-validation (default 7). Set
#'   to 0 to skip CV. When \code{cv.k >= 2}, \code{Q2} and \code{DQ2} in the
#'   output reflect cross-validated predictive ability; otherwise they reflect
#'   training-set fit. For each fold, test samples are first ort-deflated using
#'   \code{P.loadings.ort} from the training fold model before prediction.
#' @return A \code{ComDim} object.  All slots are populated.  Key slots:
#' \describe{
#'   \item{\code{Method}}{\code{"OPLS-DA"} or \code{"OPLS-R"}.}
#'   \item{\code{ndim}}{Number of predictive Common Dimensions.}
#'   \item{\code{Q.scores}}{Predictive global scores matrix (\eqn{n \times
#'     ndim}).  These are the KOPLS predictive scores \eqn{\mathbf{T}_p}
#'     computed from the RV-weighted double-centred kernel after orthogonal
#'     deflation.}
#'   \item{\code{T.scores}}{Named list of block-specific predictive local
#'     scores (\eqn{n \times ndim} each), derived from the orthogonally-
#'     deflated blocks.}
#'   \item{\code{P.loadings}}{Global predictive loadings (\eqn{p_{tot}
#'     \times ndim}):
#'     \eqn{\mathbf{P} = \tilde{\mathbf{X}}'\mathbf{Q}},
#'     where \eqn{\tilde{\mathbf{X}}} is the orthogonally-deflated
#'     concatenated X.}
#'   \item{\code{Saliences}}{Predictive block salience matrix (\eqn{ntable
#'     \times ndim}): RV-weighted kernel-based saliences per block per
#'     predictive component.}
#'   \item{\code{Orthogonal}}{List with orthogonal component outputs:
#'     \describe{
#'       \item{\code{nort}}{Number of orthogonal components.}
#'       \item{\code{Q.scores}}{Global orthogonal scores (\eqn{n \times nort}),
#'         unit-norm columns.}
#'       \item{\code{T.scores}}{Named list of block-specific orthogonal local
#'         scores (\eqn{n \times nort} each).}
#'       \item{\code{P.loadings.ort}}{Orthogonal loadings (\eqn{p_{tot} \times
#'         nort}): concatenated block-slice local ort loadings.}
#'       \item{\code{Saliences.ort}}{Orthogonal block saliences (\eqn{ntable
#'         \times nort}).}
#'     }
#'   }
#'   \item{\code{R2X}}{Named vector (length \eqn{ndim + nort}) combining
#'     predictive and orthogonal X-variance fractions, using the KOPLS
#'     kernel \eqn{\mathbf{K} = \mathbf{H}
#'     \bigl(\sum_b \mathrm{RV}_b\,\mathbf{X}_b\mathbf{X}_b'/
#'     \|\mathbf{X}_b\mathbf{X}_b'\|_F\bigr) \mathbf{H}} (RV-weighted,
#'     double-centred, \eqn{\mathbf{H} = \mathbf{I} - n^{-1}\mathbf{J}}):
#'     \describe{
#'       \item{Predictive \code{CC_a}:}{
#'         \eqn{R2X_{p,a} = \|\mathbf{T}_p[,a]\|^2 / \mathrm{tr}(\mathbf{K})},
#'         the fraction of X kernel sum-of-squares captured by predictive
#'         score \eqn{a}.}
#'       \item{Orthogonal \code{ort_a}:}{
#'         \eqn{R2X_{o,a} =
#'         (\mathrm{tr}(\mathbf{K}_{a-1}) - \mathrm{tr}(\mathbf{K}_a)) /
#'         \mathrm{tr}(\mathbf{K}_0)},
#'         the fraction of kernel SS removed by the \eqn{a}-th orthogonal
#'         deflation.}
#'     }
#'   }
#'   \item{\code{R2Y}}{Named vector (length \eqn{ndim + nort}): marginal
#'     Y-variance fractions.  For predictive components the first entry is the
#'     cumulative R2Y with one predictive component (before any ort deflation);
#'     subsequent entries are marginal increments.  Orthogonal entries
#'     are the change in R2Y when each ort deflation is applied (by
#'     construction \eqn{\approx 0} since ort components are orthogonal
#'     to Y).}
#'   \item{\code{Q2}}{Cross-validated Q2 per class/response (when
#'     \code{cv.k >= 2}; otherwise training-set fit), named by class:
#'     \deqn{Q2 = 1 - PRESS / TSS_Y,}
#'     where \eqn{PRESS = \sum_i (\hat{y}_i - y_i)^2} and
#'     \eqn{TSS_Y = \sum_i (y_i - \bar{y})^2}, computed from the full-model
#'     CV level (all ort deflations applied).}
#'   \item{\code{DQ2}}{(OPLS-DA only) Cross-validated discriminant Q2 per
#'     class, using only penalising residuals:
#'     \deqn{DQ2 = 1 - PRESSD / TSS_Y,}
#'     where \eqn{PRESSD} sums \eqn{\hat{y}_i^2} for class-0 samples with
#'     \eqn{\hat{y}_i > 0}, and \eqn{(\hat{y}_i - 1)^2} for class-1 samples
#'     with \eqn{\hat{y}_i < 1}.}
#'   \item{\code{VIP}}{Global total VIP (named vector, length \eqn{p_{tot}}):
#'     concatenation of \code{VIP.block[[b]]$tot} across blocks:
#'     \eqn{VIPtot_j = \sqrt{(VIPp_j^2 + VIPo_j^2)/2}}.}
#'   \item{\code{VIP.block}}{Named list (one \code{data.frame} per block)
#'     with columns:
#'     \describe{
#'       \item{\code{p}}{Predictive VIP (Y-based).
#'         Let
#'         \eqn{\mathbf{Q}_s = \mathbf{Y}'\mathbf{T}_p\,
#'         \mathrm{diag}(1/\|\mathbf{T}_p[,a]\|^2)},
#'         \eqn{\mathbf{U}_s = \mathbf{Y}\mathbf{Q}_s\,
#'         \mathrm{diag}(1/\|\mathbf{Q}_s[,a]\|^2)},
#'         and \eqn{\tilde{\mathbf{W}}_s} = column-L2-normalised
#'         \eqn{\mathbf{X}_b'\mathbf{U}_s\,
#'         \mathrm{diag}(1/\|\mathbf{U}_s[,a]\|^2)};
#'         with \eqn{s_a = \|\mathbf{T}_p[,a]\|^2\|\mathbf{Q}_s[,a]\|^2}:
#'         \deqn{VIPp_j = \sqrt{p_b \cdot
#'           \frac{\sum_a s_a \tilde{W}_{s,j,a}^2}{\sum_a s_a}}.}}
#'       \item{\code{o}}{Orthogonal VIP (loadings-based).
#'         Let \eqn{\tilde{\mathbf{P}}_o} be the column-L2-normalised
#'         block-slice of the orthogonal loadings and
#'         \eqn{s_a = \|\mathbf{q}_{ort}[,a]\|^2} (\eqn{= 1} for unit-norm
#'         ort scores):
#'         \deqn{VIPo_j = \sqrt{p_b \cdot
#'           \frac{\sum_a s_a \tilde{P}_{o,j,a}^2}{\sum_a s_a}}.}}
#'       \item{\code{tot}}{Total VIP:
#'         \eqn{VIPtot_j = \sqrt{(VIPp_j^2 + VIPo_j^2)/2}}.}
#'     }
#'     Row names are variable names.
#'   }
#'   \item{\code{PLS.model}}{KOPLS regression objects: \code{W}
#'     (variable-space prediction weight matrix, \eqn{p_{tot} \times ndim},
#'     \eqn{\mathbf{W}_{pred}[b] =
#'     \tilde{\mathbf{X}}_b'\mathbf{U}_p\mathbf{S}_p\,
#'     \mathrm{RV}_b / \|\mathbf{X}_b\mathbf{X}_b'\|_F});
#'     \code{B} (regression coefficients,
#'     \eqn{\mathbf{B} =
#'     \mathbf{W}_{pred}\mathbf{B}_t\mathbf{C}_p'},
#'     where \eqn{\mathbf{B}_t = (\mathbf{T}_p'\mathbf{T}_p)^{-1}
#'     \mathbf{T}_p'\mathbf{U}_p} and \eqn{\mathbf{C}_p} are Y loadings
#'     from SVD(\eqn{\mathbf{Y}'\mathbf{K}\mathbf{Y}}));
#'     \code{B0} (intercept,
#'     \eqn{\mathbf{B}_0 = \bar{\mathbf{y}} -
#'     \overline{\tilde{\mathbf{x}}}\mathbf{B}});
#'     \code{Y} (response matrix as supplied).
#'     Training-set predictions:
#'     \eqn{\hat{\mathbf{Y}} =
#'     \tilde{\mathbf{X}}\mathbf{B} + \mathbf{B}_0},
#'     where \eqn{\tilde{\mathbf{X}}} is the orthogonally-deflated X.}
#'   \item{\code{cv}}{Cross-validation results when \code{cv.k >= 2} (empty
#'     list otherwise): \code{k} (number of folds), \code{Ypred}
#'     (\eqn{n \times ncol(Y)} full-model CV predictions, all ort deflations
#'     applied), \code{Ypred.p} (pred-only CV predictions, no ort deflation),
#'     \code{Q2} (CV Q2 per class/response), \code{Q2.levels} (overall Q2
#'     per ort level \code{"p"}, \code{"po1"}, \ldots, for diagnostics),
#'     \code{DQ2} (CV DQ2 per class, OPLS-DA only).}
#'   \item{\code{Prediction}}{Training-set predictions: \code{Y.pred}; for
#'     OPLS-DA also \code{decisionRule}, \code{trueClass}, \code{predClass},
#'     \code{Sensitivity}, \code{Specificity}, \code{confusionMatrix}.}
#'   \item{\code{Mean}}{List with \code{MeanMB} (column means per block) and
#'     \code{MeanY} (column means of Y).}
#'   \item{\code{Norm}}{List with \code{NormMB} (block Frobenius norms),
#'     \code{FrobNorms} (Frobenius norms of block kernels, used in KOPLS
#'     weight computation), and \code{RVweights} (per-block modified RV
#'     coefficients used to combine block kernels).}
#'   \item{\code{variable.block}}{Block membership of each row in
#'     \code{P.loadings} and each element of \code{VIP}.}
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
                        normalise = FALSE, threshold = 1e-10,
                        loquace = FALSE,
                        cv.k = 7) {
  # INITIALISATION

  if (!(method == "OPLS-DA" || method == "OPLS-R")) {
    stop("'method' must be 'OPLS-DA' or 'OPLS-R'.")
  }

  if (nort < 1) {
    stop("'nort' must be at least 1.")
  }

  if (ndim < 1) {
    stop("'ndim' must be at least 1.")
  }

  progress_bar <- utils::txtProgressBar(min = 0, max = 80, style = 3, char = "=")

  start_ini_time <- Sys.time() # To start counting calculation time for the initialization

  if (!inherits(MB, "MultiBlock")) {
    stop("'MB' is  not a MultiBlock.")
  }

  ntable <- length(blockNames(MB)) # number of blocks
  nrowMB <- length(sampleNames(MB)) # number of samples
  variable_number <- ncol(MB) # Number of variables per block

  give_error <- 0

  # Remove all variables with 0 variance. If there are no remaining variables, give_error
  numvar <- 0
  numblock0 <- vector()
  for (i in 1:ntable) {
    var0 <- apply(MB@Data[[i]], 2, function(x) {
      var(x)
    })
    var0 <- which(var0 == 0)

    if (length(var0) != 0) {
      numvar <- numvar + length(var0)
      if (length(var0) == length(MB@Variables[[i]])) {
        numblock0 <- append(numblock0, i)
      }
      MB@Data[[i]] <- MB@Data[[i]][, setdiff(1:variable_number[i], var0)]
      MB@Variables[[i]] <- MB@Variables[[i]][setdiff(1:variable_number[i], var0)]
      variable_number[i] <- length(MB@Variables[[i]])
    }
  }
  if (numvar != 0) {
    warning(sprintf("Number of variables excluded from the analysis because of their zero variance: %d", numvar))
  }
  if (length(numblock0) != 0) {
    numblock0 <- rev(numblock0) # To sort in decreasing order.
    for (i in seq_along(numblock0)) {
      MB@Data[[numblock0[i]]] <- NULL
      MB@Variables[[numblock0[i]]] <- NULL
      if (names(MB@Data)[numblock0[i]] %in% names(MB@Batch)) {
        MB@Batch[[numblock0[i]]] <- NULL
      }
      if (names(MB@Data)[numblock0[i]] %in% names(MB@Metadata)) {
        MB@Metadata[[numblock0[i]]] <- NULL
      }
    }
    warning(sprintf("Number of blocks excluded from the analysis because of their zero variance: %d", length(numblock0)))
  }

  ntable <- length(blockNames(MB)) # number of blocks
  variable_number <- ncol(MB) # In case the number of blocks has changed.

  if (length(MB@Data) == 0) { # In case all blocks had 0 variance...
    warning("All blocks had 0 variance")
    give_error <- 1
  }

  # Check for infinite or missing values (they cannot be handled with svd)
  for (i in 1:ntable) {
    if (any(is.na(MB@Data[[i]]))) {
      stop("The MB contains NA values. They must be removed first for ComDim analysis.")
    }
    if (any(is.infinite(MB@Data[[i]]))) {
      stop("The MB contains infinite values. They must be removed first for ComDim analysis.")
    }
  }

  if (give_error) {
    stop("The data is not ready for ComDim.")
  } else {
    print("The data can be used for ComDim.")
  }

  ## Check y matrix (convert to dummy matrix if it is not already.)
  if (method == "OPLS-DA") {
    tmp <- unique(as.vector(y))
    if (all(tmp %in% c(0, 1))) { # Is it dummy?

      if (!is.matrix(y)) { # If it's a vector and the method is OPLS-R
        y <- as.matrix(y)
      }

      classVect <- rep(NA, nrow(y))

      if (ncol(y) == 1) {
        classVect <- as.vector(y)
      } else {
        if (any(rowSums(y) != 1)) {
          if (any(rowSums(y) == 0)) {
            stop("At least one sample was not assigned to a class.")
          } else if (any(colSums(y) > 1)) {
            stop("At least one sample was assigned to more than one class.")
          }
        }
      }
      if (is.null(colnames(y))) {
        colnames(y) <- paste0("Class", seq_len(ncol(y)))
      }
      for (i in seq_len(ncol(y))) {
        classVect[which(y[, i] == 1)] <- colnames(y)[i]
      }
    } else { # If not dummy
      if ((is.matrix(y) || is.data.frame(y)) && ncol(y) > 1) {
        stop("ComDim-OPLS can only be applied if Y is a dummy matrix or a class vector")
      }
      classVect <- as.vector(y)
      classVect_sorted <- sort(unique(classVect))

      # Now generate the dummy matrix from y
      y <- matrix(rep(0, length(classVect_sorted) * length(classVect)),
        ncol = length(classVect_sorted), nrow = length(classVect)
      )
      for (i in seq_along(classVect_sorted)) {
        y[which(classVect == classVect_sorted[i]), i] <- 1
      }
      colnames(y) <- classVect_sorted
      rm(classVect_sorted)
    }
  }

  if (!is.matrix(y)) { # If it's a vector and the method is OPLS-R
    y <- as.matrix(y)
  }

  pieceBar <- 4 + 2 * ntable + ndim + nort # Number of updates in the progress bar.
  pieceBar <- 80 / pieceBar
  total_progress <- pieceBar

  DimLabels <- paste0("CC", 1:ndim) # One label per predictive component.
  OrtLabels <- paste0("ort", 1:nort) # One label per orthogonal component.
  TableLabels <- blockNames(MB) # One label per block.

  end_ini_time <- Sys.time() # To end the count of the calculation time.

  if (loquace) {
    print(sprintf("Initialisation finished after : %s millisecs", (end_ini_time - start_ini_time) * 1000))
  }

  utils::setTxtProgressBar(progress_bar, value = total_progress)

  # NORMALISATION

  X_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))
  Xnorm_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))

  res_calib <- list()
  temp_tabCalib <- list()
  s_n <- list()
  res_calib$SingVal <- as.vector(NULL)
  res_calib$NormMB <- as.vector(NULL)
  res_calib$MeanMB <- list()

  for (i in 1:ntable) {
    res_calib$MeanMB[[TableLabels[i]]] <- colMeans(MB@Data[[i]])
    names(res_calib$MeanMB[[TableLabels[i]]]) <- MB@Variables[[i]]

    if (normalise) {
      # Normalise original blocks
      X_mean <- MB@Data[[i]] - matrix(data = rep(1, nrowMB), ncol = 1, nrow = nrowMB) %*% res_calib$MeanMB[[i]]
      XX <- X_mean * X_mean
      Norm_X <- sqrt(sum(XX))
      X_Normed <- X_mean / Norm_X
      res_calib$NormMB[i] <- Norm_X
      temp_tabCalib[[i]] <- X_Normed
      s_n[[i]] <- X_Normed
    } else {
      res_calib$NormMB[[TableLabels[i]]] <- rep(1, length(MB@Variables[[i]]))
      names(res_calib$NormMB[[TableLabels[i]]]) <- MB@Variables[[i]]
      temp_tabCalib[[i]] <- MB@Data[[i]]
      s_n[[i]] <- MB@Data[[i]]
    }

    if (i == 1) {
      X_mat[, 1:variable_number[1]] <- MB@Data[[i]]
      Xnorm_mat[, 1:variable_number[1]] <- temp_tabCalib[[i]]
    } else {
      beg <- sum(variable_number[1:(i - 1)]) + 1
      ending <- sum(variable_number[1:i])
      X_mat[, beg:ending] <- MB@Data[[i]]
      Xnorm_mat[, beg:ending] <- temp_tabCalib[[i]]
    }

    # PLS addition
    meanX <- colMeans(Xnorm_mat)
    meanY <- colMeans(y)
    # End PLS addition

    norm_comdim <- Sys.time()
    if (loquace) {
      print(sprintf("Normalization of block %s finished after : %s millisecs", i, (norm_comdim - start_ini_time) * 1000))
    }
    total_progress <- total_progress + pieceBar
    utils::setTxtProgressBar(progress_bar, value = total_progress)
  }

  names(res_calib$NormMB) <- TableLabels

  nR <- nrow(Xnorm_mat)
  nC <- ncol(Xnorm_mat)

  # =========================================================================
  # KOPLS KERNEL-BASED EXTRACTION 
  # =========================================================================

  end_comp_time <- Sys.time()
  if (loquace) {
    print(sprintf(
      "Kernel preparation finished after : %s millisecs",
      (end_comp_time - start_ini_time) * 1000
    ))
  }
  total_progress <- total_progress + pieceBar * ntable
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  # Helper: modified RV coefficient (zeros diagonal before computing)
  .RVmod <- function(X, Y) {
    AA <- tcrossprod(X)
    diag(AA) <- 0
    BB <- tcrossprod(Y)
    diag(BB) <- 0
    denom <- sqrt(sum(AA^2)) * sqrt(sum(BB^2))
    if (denom < .Machine$double.eps) {
      return(0)
    }
    sum(diag(crossprod(AA, BB))) / denom
  }

  K_blocks <- vector("list", ntable)
  frob_norms <- numeric(ntable) # Frobenius norm of each raw block kernel
  rv_weights <- numeric(ntable) # RV weight for each block

  for (i in seq_len(ntable)) {
    Xi <- temp_tabCalib[[i]]
    if (normalise) {
      ss <- sum(Xi^2)
      if (ss > 1e-14) Xi <- Xi / sqrt(ss)
    }
    Ki <- tcrossprod(Xi)
    frob_Ki <- sqrt(sum(Ki^2))
    frob_norms[i] <- if (frob_Ki > 1e-14) frob_Ki else 1
    K_norm_i <- Ki / frob_norms[i]
    rv_weights[i] <- (.RVmod(K_norm_i, y) + 1) / 2
    K_blocks[[i]] <- K_norm_i
  }
  names(rv_weights) <- TableLabels

  # RV-weighted combined kernel
  W_mat <- Reduce("+", lapply(
    seq_len(ntable),
    function(i) rv_weights[i] * K_blocks[[i]]
  ))

  # Double-centre: K_combined = H * W_mat * H  (H = I - 1/n * J)
  K_combined <- scale(t(scale(W_mat, scale = FALSE)), scale = FALSE)

  rm(s_n)

  # --- 2. Centre Y; SVD of Y' K Y -------------------------------------------
  Y_c <- sweep(y, 2, colMeans(y)) # n x ncol(y)
  YKY <- crossprod(Y_c, K_combined %*% Y_c) # ncol(y) x ncol(y)

  A <- ndim
  CSV <- svd(YKY, nu = A, nv = A)
  Cp <- CSV$u[, 1:A, drop = FALSE] # ncol(y) x ndim  Y loadings
  Sp <- diag(CSV$d[1:A]^(-0.5), nrow = A) # ndim x ndim     scale
  Up <- Y_c %*% Cp # n x ndim        Y scores

  colnames(Cp) <- DimLabels
  rownames(Cp) <- colnames(y)
  colnames(Up) <- DimLabels
  rownames(Up) <- MB@Samples

  # --- 3. Orthogonal extraction (KOPLS loop) --------------------------------
  Q_ort <- matrix(NA_real_, nrow = nrowMB, ncol = nort)
  saliences_ort <- matrix(NA_real_, nrow = ntable, ncol = nort)

  to_list <- vector("list", nort) # unit-norm ort scores
  co_list <- vector("list", nort) # ort direction in Y-score space
  so_list <- vector("list", nort) # singular value
  toNorm_list <- vector("list", nort) # norm before unit-normalisation

  # Per-component storage for R2X and R2Y
  # Output will have ndim + nort values: CC1..CC{ndim}, ort1..ort{nort}
  sstot_K      <- sum(diag(K_combined))    # trace(K_original) = SS(X)
  Kdiag_levels <- numeric(nort + 1)        # trace(K_both) before/after each ort deflation
  Kdiag_levels[1] <- sstot_K              # before any ort deflation
  Tp_levels    <- vector("list", nort + 1) # predictive scores at each ort level

  K_right <- K_combined # n x n  right-only deflated kernel
  K_both <- K_combined # n x n  both-sides deflated kernel

  for (jj in seq_len(nort)) {
    Tp_curr <- crossprod(K_right, Up %*% Sp) # = t(K_right) %*% Up %*% Sp
    Tp_levels[[jj]] <- Tp_curr              # Tp using K deflated by (jj-1) ort components
    K_resid <- K_both - tcrossprod(Tp_curr) # K - Tp Tp'

    sv_tmp <- svd(crossprod(Tp_curr, K_resid %*% Tp_curr), nu = 1, nv = 1)
    co_j <- sv_tmp$u[, 1, drop = FALSE] # ndim x 1
    so_j <- sv_tmp$d[1]

    to_raw <- K_resid %*% Tp_curr %*% co_j /
      sqrt(max(so_j, .Machine$double.eps)) # n x 1
    toNorm_j <- sqrt(sum(to_raw^2))
    if (toNorm_j < 1e-14) toNorm_j <- 1
    to_j <- to_raw / toNorm_j # unit-norm

    to_list[[jj]] <- to_j
    co_list[[jj]] <- co_j
    so_list[[jj]] <- so_j
    toNorm_list[[jj]] <- toNorm_j
    Q_ort[, jj] <- as.vector(to_j)

    for (i in seq_len(ntable)) {
      saliences_ort[i, jj] <-
        as.vector(crossprod(to_j, K_blocks[[i]]) %*% to_j)
    }

    # Kernel deflation
    I_to <- diag(nrowMB) - tcrossprod(to_j)
    K_right <- K_right %*% I_to # right only
    K_both <- I_to %*% K_both %*% I_to # both sides
    Kdiag_levels[jj + 1] <- sum(diag(K_both)) # trace after jj ort deflations

    ort_comdim <- Sys.time()
    if (loquace) {
      print(sprintf(
        "Orthogonal component %s determined after : %s millisecs",
        jj, (ort_comdim - start_ini_time) * 1000
      ))
    }
    total_progress <- total_progress + pieceBar
    utils::setTxtProgressBar(progress_bar, value = total_progress)
  }

  colnames(Q_ort) <- OrtLabels
  rownames(Q_ort) <- MB@Samples
  colnames(saliences_ort) <- OrtLabels
  rownames(saliences_ort) <- TableLabels

  # --- 4. Final predictive scores -------------------------------------------
  # Use t(K_right) (= crossprod) because K_right is the right-only deflated
  # kernel and is no longer symmetric after deflation. 
  Tp <- crossprod(K_right, Up %*% Sp) # n x ndim
  Tp_levels[[nort + 1]] <- Tp         # final predictive scores (all ort deflations applied)

  # Inner regression coefficient (pinv handles rank-deficient Tp when ndim > rank(Y))
  Bt <- pracma::pinv(crossprod(Tp)) %*% crossprod(Tp, Up) # ndim x ndim

  # Predictive block saliences
  saliences <- matrix(NA_real_, nrow = ntable, ncol = ndim)
  for (i in seq_len(ntable)) {
    for (j in seq_len(ndim)) {
      tp_j <- Tp[, j]
      saliences[i, j] <-
        as.vector(crossprod(tp_j, K_blocks[[i]]) %*% tp_j)
    }
  }
  colnames(saliences) <- DimLabels
  rownames(saliences) <- TableLabels

  # SingVal: SS of predictive scores (kept for compatibility)
  res_calib$SingVal <- setNames(colSums(Tp^2), DimLabels)

  # R2X per component (non-cumulative):
  #   CC_j  = fraction of X kernel SS captured by j-th predictive score
  #           = sum(Tp[:,j]^2) / trace(K_original)
  #   ort_j = fraction of X kernel SS removed by j-th ort deflation
  #           = (trace(K_before_j) - trace(K_after_j)) / trace(K_original)
  r2x_pred        <- setNames(colSums(Tp^2) / sstot_K, DimLabels)
  r2x_ort         <- setNames(-diff(Kdiag_levels) / sstot_K, OrtLabels)
  res_calib$explained <- c(r2x_pred, r2x_ort)

  # Store predictive scores as Q
  Q <- Tp
  colnames(Q) <- DimLabels
  rownames(Q) <- MB@Samples
  res_calib$Q <- Q

  ## Adding metadata
  res_calib$metadata <- list()
  if (length(MB@Metadata) != 0) {
    res_calib$metadata[[1]] <- MB@Metadata
  }

  iter_comdim <- Sys.time()
  if (loquace) {
    print(sprintf(
      "KOPLS scores finished after : %s millisecs",
      (iter_comdim - start_ini_time) * 1000
    ))
  }
  total_progress <- total_progress + pieceBar * (1 + ndim)
  utils::setTxtProgressBar(progress_bar, value = min(total_progress, 79))

  rm(K_combined, K_right, K_both)

  # --- 5. R2Y per component --------------------------------------------------
  # Predictive CC_j: incremental Y variance explained by j successive predictive
  #   components using Tp BEFORE any ort deflation (Tp_levels[[1]]).
  # Orthogonal ort_j: marginal change in R2Y when applying j-th ort deflation
  #   (by design ~0, since ort components are orthogonal to Y).
  sstot_Y_r2 <- sum(Y_c^2)

  # Cumulative R2Y using the first j columns of Tp before ort deflation; take diffs
  Tp0 <- Tp_levels[[1]]  # pred scores before any ort deflation
  r2y_pred_cum <- numeric(ndim)
  for (.j in seq_len(ndim)) {
    Tp_j   <- Tp0[, seq_len(.j), drop = FALSE]
    Bt_j   <- pracma::pinv(crossprod(Tp_j)) %*% crossprod(Tp_j, Up[, seq_len(.j), drop = FALSE])
    Yhat_j <- Tp_j %*% Bt_j %*% t(Cp[, seq_len(.j), drop = FALSE])
    r2y_pred_cum[.j] <- if (sstot_Y_r2 < .Machine$double.eps) 0 else
      1 - sum((Yhat_j - Y_c)^2) / sstot_Y_r2
  }
  r2y_pred <- setNames(c(r2y_pred_cum[1], diff(r2y_pred_cum)), DimLabels)

  # R2Y at each ort level (all ndim pred components); marginal diff per ort step
  r2y_lev <- numeric(nort + 1)
  for (.i in seq_len(nort + 1)) {
    Tp_i      <- Tp_levels[[.i]]
    Bt_i      <- pracma::pinv(crossprod(Tp_i)) %*% crossprod(Tp_i, Up)
    Yhat_i    <- Tp_i %*% Bt_i %*% t(Cp)
    r2y_lev[.i] <- if (sstot_Y_r2 < .Machine$double.eps) 0 else
      1 - sum((Yhat_i - Y_c)^2) / sstot_Y_r2
  }
  r2y_ort <- setNames(diff(r2y_lev), OrtLabels)

  r2y <- c(r2y_pred, r2y_ort)

  ## Calculate Sums of saliences for each Dim
  Sum_saliences_Dim <- colSums(saliences)
  names(Sum_saliences_Dim) <- DimLabels
  res_calib$Sum_saliences_Dim <- Sum_saliences_Dim
  rm(Sum_saliences_Dim)

  # Calculate Sums of saliences for each Table
  Sum_saliences_Tab <- rowSums(saliences)
  names(Sum_saliences_Tab) <- TableLabels
  res_calib$Sum_saliences_Tab <- Sum_saliences_Tab
  rm(Sum_saliences_Tab)

  res_calib$saliences <- saliences
  colnames(res_calib$saliences) <- DimLabels
  rownames(res_calib$saliences) <- TableLabels

  ## LOADING COMPUTATION

  L_X <- list()
  T_Loc_ort <- list()
  T_Loc <- list()
  for (i in 1:ntable) { # Prepare lists
    T_Loc_ort[[TableLabels[i]]] <- matrix(, ncol = nort, nrow = nrowMB)
    T_Loc[[TableLabels[i]]] <- matrix(, ncol = ndim, nrow = nrowMB)
  }
  b <- matrix(, ncol = ndim, nrow = ntable)

  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("The package pracma is needed.")
  }

  # Save pre-deflation blocks for VIP computation
  X_orig_blocks <- temp_tabCalib

  # Orthogonal loadings: deflate temp_tabCalib block by block, compute local
  # orthogonal scores and the concatenated orthogonal loading matrix (orthoP).
  orthoP <- NULL
  for (jj in 1:nort) {
    for (i in 1:ntable) {
      tempo <- t(temp_tabCalib[[i]]) %*% Q_ort[, jj] # Scaled CD 'local' Loadings
      T_Loc_ort[[TableLabels[i]]][, jj] <- temp_tabCalib[[i]] %*% (tempo %*% pracma::pinv(t(tempo) %*% tempo)) # local Scores
      # Deflate each temp_tabCalib by orthogonal component
      temp_tabCalib[[i]] <- temp_tabCalib[[i]] - Q_ort[, jj] %*% t(tempo)
      if (i == 1) {
        tempo2 <- tempo
      } else {
        tempo2 <- append(tempo2, tempo)
      }
    }
    if (jj == 1) {
      orthoP <- tempo2
    } else {
      orthoP <- cbind(orthoP, tempo2)
    }
  }
  orthoP <- as.matrix(orthoP)
  colnames(orthoP) <- OrtLabels

  for (i in 1:ntable) {
    rownames(T_Loc_ort[[TableLabels[i]]]) <- MB@Samples
    colnames(T_Loc_ort[[TableLabels[i]]]) <- OrtLabels
  }

  # Rebuild Xnorm_mat from orthogonally-deflated temp_tabCalib.
  # This deflated X is used for global loadings and Y prediction.
  for (i in 1:ntable) {
    if (i == 1) {
      Xnorm_mat <- temp_tabCalib[[i]]
    } else {
      Xnorm_mat <- cbind(Xnorm_mat, temp_tabCalib[[i]])
    }
  }

  # Predictive loadings: compute local predictive scores and b-coefficients
  # from the orthogonally-deflated temp_tabCalib.
  for (j in 1:ndim) {
    T_mat <- matrix(, nrow = nrowMB, ncol = 0)

    # Unit-norm predictive score used for variable-space deflation.
    # KOPLS Tp has eigenvalue-based scale; unit-normalise so that the
    # deflation X_new = X - q_unit*(X'q_unit)' is an orthogonal projection.
    q_j <- res_calib$Q[, j]
    q_j_ss <- sum(q_j^2)
    q_j_unit <- if (q_j_ss > 1e-14) q_j / sqrt(q_j_ss) else q_j

    for (i in 1:ntable) {
      temp <- t(temp_tabCalib[[i]]) %*% q_j_unit # Scaled CD 'local' Loadings

      T_mat <- cbind(T_mat, temp_tabCalib[[i]] %*% (temp %*% pracma::pinv(t(temp) %*% temp))) # local Scores

      T_Loc[[TableLabels[i]]][, j] <- temp_tabCalib[[i]] %*% (temp %*% pracma::pinv(t(temp) %*% temp)) # local Scores

      # Deflate each temp_tabCalib (unit-norm direction for correct projection)
      temp_tabCalib[[i]] <- temp_tabCalib[[i]] - q_j_unit %*% t(temp)
    }

    # For each CC
    # MLR b-coefficients between Local and Global Scores
    # b=pinv(X'*X)*X'*Y;
    b[, j] <- pracma::pinv(t(T_mat) %*% T_mat) %*% t(T_mat) %*% res_calib$Q[, j]
  }

  for (i in 1:ntable) {
    rownames(T_Loc[[TableLabels[i]]]) <- rownames(res_calib$Q)
    colnames(T_Loc[[TableLabels[i]]]) <- colnames(res_calib$Q)
  }

  # Calculate Global Loadings from orthogonally-deflated X
  L_CD_Vec <- t(Xnorm_mat) %*% res_calib$Q # Scaled CD 'global' Loadings

  load_comdim <- Sys.time()
  if (loquace) {
    print(sprintf("Loadings finished after : %s millisecs", (load_comdim - start_ini_time) * 1000))
  }

  total_progress <- total_progress + pieceBar
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  rm(X_mat)
  rm(T_mat)
  rm(temp_tabCalib)
  rm(temp)

  ## Output
  end_function <- Sys.time()

  ## KOPLS B matrix (variable-space regression coefficient)
  # W_pred is computed block-by-block with the Frobenius normalization that was
  # applied to each block kernel: W_b = t(X_ort_b) * Up * Sp / frob_norms[b].
  # This ensures that X_test %*% W_pred = K_test_train_norm %*% Up %*% Sp, i.e.,
  # B  = W_pred * Bt * Cp'    [p x ncol(y)]
  # B0 = colMeans(y) - colMeans(X_ort) * B
  UpSp <- Up %*% Sp # n x ndim
  W_pred <- matrix(NA_real_, nrow = sum(variable_number), ncol = ndim)
  k_idx <- 1
  for (i in seq_len(ntable)) {
    vi <- variable_number[i]
    Xb_ort <- Xnorm_mat[, k_idx:(k_idx + vi - 1), drop = FALSE]
    W_pred[k_idx:(k_idx + vi - 1), ] <- t(Xb_ort) %*% UpSp * rv_weights[i] / frob_norms[i]
    k_idx <- k_idx + vi
  }
  B_kopls <- W_pred %*% Bt %*% t(Cp) # p x ncol(y)
  B0_kopls <- colMeans(y) - colMeans(Xnorm_mat) %*% B_kopls # 1 x ncol(y)
  colnames(W_pred) <- DimLabels
  rm(UpSp)

  res_calib$b <- b
  colnames(res_calib$b) <- DimLabels
  rownames(res_calib$b) <- TableLabels

  res_calib$T_Loc <- T_Loc
  rm(T_Loc)

  varnames <- vector()
  for (i in 1:ntable) {
    varnames <- append(varnames, MB@Variables[[i]])
  }

  res_calib$P <- L_CD_Vec
  colnames(res_calib$P) <- DimLabels
  rownames(res_calib$P) <- varnames
  rownames(orthoP) <- varnames

  rm(L_CD_Vec)

  ## -------------------------------------------------------------------------
  ## VIP (Variable Importance in Projection)
  ##
  ## Predictive VIP ("p") — Y-based formula:
  ##   Qs  = t(Y) %*% Tp %*% diag(1/colSums(Tp^2))
  ##   Us  = Y   %*% Qs  %*% diag(1/colSums(Qs^2))
  ##   Ws  = t(X_b) %*% Us %*% diag(1/colSums(Us^2)), column-L2-normalised
  ##   s_a = colSums(Tp^2) * colSums(Qs^2)
  ##   VIPp_j = sqrt( p_b * sum(s_a * Ws_j^2) / sum(s_a) )
  ##
  ## Orthogonal VIP ("o") — loadings-based formula (orthogonal scores ⊥ Y so
  ##   the Y-based formula collapses; instead weight each variable by its
  ##   contribution to the orthogonal loadings):
  ##   s_a = colSums(Q_ort^2)  [= 1 for unit-norm ort scores]
  ##   P_norm = column-L2-normalised block-slice of orthoP
  ##   VIPo_j = sqrt( p_b * sum(s_a * P_norm_j^2) / sum(s_a) )
  ##
  ## Total VIP ("tot") = sqrt((p^2 + o^2) / 2)
  ## Global VIP (@VIP) = "tot" concatenated across all blocks.
  ## -------------------------------------------------------------------------

  .vip_for_scores <- function(scores_mat, Y_raw, X_b, p_b) {
    ncomp <- ncol(scores_mat)
    ss_T <- colSums(scores_mat^2)
    if (any(ss_T < .Machine$double.eps)) {
      return(rep(NA_real_, p_b))
    }
    Qs <- crossprod(Y_raw, scores_mat) %*% diag(1 / ss_T, nrow = ncomp)
    ss_Q <- colSums(Qs^2)
    if (any(ss_Q < .Machine$double.eps)) {
      return(rep(NA_real_, p_b))
    }
    Us <- Y_raw %*% Qs %*% diag(1 / ss_Q, nrow = ncomp)
    ss_U <- colSums(Us^2)
    if (any(ss_U < .Machine$double.eps)) {
      return(rep(NA_real_, p_b))
    }
    Ws <- crossprod(X_b, Us) %*% diag(1 / ss_U, nrow = ncomp)
    Ws <- apply(Ws, 2, function(w) {
      nrm <- sqrt(sum(w^2))
      if (nrm > 1e-14) w / nrm else w
    })
    s_a <- ss_T * ss_Q
    sum_s <- sum(s_a)
    if (sum_s < .Machine$double.eps) {
      return(rep(NA_real_, p_b))
    }
    as.vector(sqrt(p_b * (Ws^2 %*% s_a) / sum_s))
  }

  # Orthogonal VIP uses ort-loading magnitudes (not Y) since ort scores ⊥ Y.
  .vip_ort_for_loadings <- function(ort_scores_mat, ort_loadings_b) {
    ncomp <- ncol(ort_scores_mat)
    p_b   <- nrow(ort_loadings_b)
    ss_T  <- colSums(ort_scores_mat^2)
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
  k_ort <- 1L  # running row index into orthoP for the current block
  for (i in seq_len(ntable)) {
    X_b   <- X_orig_blocks[[i]]
    p_b   <- variable_number[i]
    vip_p <- .vip_for_scores(Tp, y, X_b, p_b)
    ort_P_b <- orthoP[k_ort:(k_ort + p_b - 1L), , drop = FALSE]
    vip_o <- .vip_ort_for_loadings(Q_ort, ort_P_b)
    vip_t <- sqrt((vip_p^2 + vip_o^2) / 2)
    vdf <- data.frame(
      p = vip_p, o = vip_o, tot = vip_t,
      row.names = MB@Variables[[i]]
    )
    VIP_blocks[[TableLabels[i]]] <- vdf
    k_ort <- k_ort + p_b
  }

  VIP_global <- as.vector(unlist(lapply(VIP_blocks, function(b) b$tot)))
  names(VIP_global) <- varnames
  ## -------------------------------------------------------------------------

  # SingVal and R2X are computed from scores (set above).

  rm(saliences)

  # Define block length
  belong_block <- rep(0, nrow(res_calib$P))
  k <- 1
  for (i in 1:ntable) {
    belong_block[k:(k + variable_number[i] - 1)] <- TableLabels[i]
    k <- k + variable_number[i]
  }

  rm(Q)

  ## Y prediction using the orthogonally-deflated X and KOPLS B matrix
  Ypred <- Xnorm_mat %*% B_kopls +
    matrix(as.vector(B0_kopls), nrow = nrowMB, ncol = ncol(y), byrow = TRUE)
  rownames(Ypred) <- MB@Samples
  if (!is.null(colnames(y))) colnames(Ypred) <- colnames(y)

  ## Predict classes
  if (method == "OPLS-DA") {
    predClass <- rep(NA, nrow(y))

    if (decisionRule == "max") {
      for (i in seq_len(nrow(Ypred))) {
        xx <- which(Ypred[i, ] == max(Ypred[i, ], na.rm = TRUE))
        if (length(xx) == 1) {
          predClass[i] <- colnames(y)[xx] # Y stores the Class names
        } else {
          predClass[i] <- NaN # There is a tie.
          warning(sprintf("Class for sample %s could not be predicted.", i))
        }
      }
    } else if (decisionRule == "fixed") {
      for (i in seq_len(nrow(Ypred))) {
        xx <- which(Ypred[, i] > 1 / ncol(y))
        if (length(xx) == 1) {
          predClass[i] <- colnames(y)[xx] # Y stores the Class names
        } else {
          predClass[i] <- NaN # There is a tie.
          warning(sprintf("Class for sample %s could not be predicted.", i))
        }
      }
    }

    if (is.matrix(y)) {
      meanY <- colMeans(y)
    } else {
      meanY <- mean(y)
    }

    classy <- colnames(y)
    Q2 <- DQ2 <- specvec <- sensvec <- rep(NA, length(classy))
    confusionMatrix <- list()

    for (i in seq_along(classy)) {
      pos <- which(classVect == classy[i])
      neg <- which(classVect != classy[i])
      tp <- length(which(classVect == classy[i] & predClass == classy[i]))
      tn <- length(which(classVect != classy[i] & predClass != classy[i]))
      fp <- length(which(classVect != classy[i] & predClass == classy[i]))
      fn <- length(which(classVect == classy[i] & predClass != classy[i]))

      if (length(tp) == 0) {
        tp <- 0
      }

      if (length(tn) == 0) {
        tn <- 0
      }

      if (length(fp) == 0) {
        fp <- 0
      }

      if (length(fn) == 0) {
        fn <- 0
      }

      sensvec[i] <- tp / length(pos)
      specvec[i] <- tn / length(neg)

      cm <- matrix(c(tp, fp, fn, tn), ncol = 2, nrow = 2)
      rownames(cm) <- c("predClass1", "predClass0")
      colnames(cm) <- c("trueClass1", "trueClass0")
      confusionMatrix[[classy[i]]] <- cm

      ## DQ2
      E0 <- Ypred[neg, i] - y[neg, i] # Calculate Residuals of Class0 samples
      E1 <- Ypred[pos, i] - y[pos, i] # Calculate Residuals of Class1 samples

      E0count <- which(E0 > 0) # Find predictions for Class0 samples larger than 0
      E1count <- which(E1 < 0) # Find predictions for Class1 samples larger than 1

      SSE0_all <- t(E0) %*% E0
      SSE1_all <- t(E1) %*% E1

      SSE0 <- t(E0[E0count]) %*% E0[E0count]
      SSE1 <- t(E1[E1count]) %*% E1[E1count]

      PRESSD <- SSE0 + SSE1
      PRESS <- SSE0_all + SSE1_all

      Ym <- y[, i] - mean(y[, i])
      TSS <- t(Ym) %*% Ym

      DQ2[i] <- 1 - (PRESSD / TSS)
      Q2[i] <- 1 - PRESS / TSS
    }
    names(Q2)      <- classy
    names(DQ2)     <- classy
    names(sensvec) <- classy
    names(specvec) <- classy
  } else if (method == "OPLS-R") {
    PRESS <- t(Ypred - y) %*% (Ypred - y)
    Ym <- y - mean(y)
    TSS <- t(Ym) %*% Ym
    Q2 <- 1 - PRESS / TSS
  }


  # --- Internal kernel centering helpers ---
  .cv_center_KtrTr <- function(K) scale(t(scale(K, scale = FALSE)), scale = FALSE)

  .cv_center_KteTr <- function(KteTr, KtrTr) {
    ntr <- nrow(KtrTr)
    nte <- nrow(KteTr)
    sc <- matrix(1 / ntr, nte, ntr)
    (KteTr - sc %*% KtrTr) %*% (diag(ntr) - (1 / ntr) * matrix(1, ntr, ntr))
  }

  # --- Internal KOPLS model fit on a centred kernel ---
  # Returns all objects needed by .cv_kopls_predict.
  .cv_kopls_fit <- function(K_cc, Yc, a, nox) {
    n <- nrow(K_cc)
    CSV <- svd(crossprod(Yc, K_cc %*% Yc), nu = a, nv = a)
    Cp <- CSV$u[, seq_len(a), drop = FALSE]
    Sp <- diag(CSV$d[seq_len(a)]^(-0.5), nrow = a)
    Up <- Yc %*% Cp
    UpSp <- Up %*% Sp

    K_right <- K_cc
    K_both <- K_cc
    to_l <- co_l <- so_l <- toNorm_l <- Tp_tr_l <- list()
    K_right_l <- list()
    K_right_l[[1]] <- K_cc

    for (jj in seq_len(nox)) {
      Tp_j <- crossprod(K_right, UpSp)
      K_resid <- K_both - tcrossprod(Tp_j)
      sv <- svd(crossprod(Tp_j, K_resid %*% Tp_j), nu = 1, nv = 1)
      co_j <- sv$u[, 1, drop = FALSE]
      so_j <- sv$d[1]
      to_raw <- K_resid %*% Tp_j %*% co_j / sqrt(max(so_j, .Machine$double.eps))
      tn <- sqrt(sum(to_raw^2))
      if (tn < 1e-14) tn <- 1
      to_j <- to_raw / tn

      Tp_tr_l[[jj]] <- Tp_j
      to_l[[jj]] <- to_j
      co_l[[jj]] <- co_j
      so_l[[jj]] <- so_j
      toNorm_l[[jj]] <- tn

      I_to <- diag(n) - tcrossprod(to_j)
      K_right <- K_right %*% I_to
      K_both <- I_to %*% K_both %*% I_to
      K_right_l[[jj + 1]] <- K_right
    }

    Tp_f <- crossprod(K_right, UpSp)
    Bt_f <- pracma::pinv(crossprod(Tp_f)) %*% crossprod(Tp_f, Up)

    list(
      Cp = Cp, Sp = Sp, Up = Up, UpSp = UpSp, Bt = Bt_f,
      to = to_l, co = co_l, so = so_l, toNorm = toNorm_l,
      Tp_tr = Tp_tr_l, K_right = K_right_l
    )
  }

  # --- Internal KOPLS predict: returns Yhat in centred Y space ---
  .cv_kopls_predict <- function(KteTr_c, m, nox) {
    KteTr_curr <- KteTr_c

    for (i in seq_len(nox)) {
      Tp_te_i <- KteTr_curr %*% m$UpSp
      K_res_tetr <- KteTr_curr - tcrossprod(Tp_te_i, m$Tp_tr[[i]])
      to_te_i <- (K_res_tetr %*% m$Tp_tr[[i]] %*% m$co[[i]]) /
        sqrt(max(m$so[[i]], .Machine$double.eps))
      to_te_i <- to_te_i / m$toNorm[[i]]
      # Deflate: KteTrdeflate[i+1,1] = KteTr_curr - to_te * t(K_right_i * to_tr)
      KteTr_curr <- KteTr_curr -
        tcrossprod(to_te_i, m$K_right[[i]] %*% m$to[[i]])
    }

    Tp_te <- KteTr_curr %*% m$UpSp
    Tp_te %*% m$Bt %*% t(m$Cp)
  }

  if (cv.k >= 2) {
    n_cv   <- nrowMB
    n_lev  <- nort + 1L   # levels: nox = 0, 1, ..., nort

    # One prediction matrix per ort level; level index = nox + 1
    Ypred_cv_all <- vector("list", n_lev)
    for (.lev in seq_len(n_lev))
      Ypred_cv_all[[.lev]] <- matrix(NA_real_, nrow = n_cv, ncol = ncol(y),
                                     dimnames = list(NULL, colnames(y)))
    pressy_cls_lev <- matrix(0, nrow = n_lev, ncol = ncol(y))  # PRESS per (ort-level x class)
    TSS_cv_cls     <- numeric(ncol(y))                           # TSS per class
    cv_test_order  <- integer(0)

    for (fold in seq_len(cv.k)) {
      pred_idx  <- seq.int(fold, n_cv, by = cv.k)
      train_idx <- setdiff(seq_len(n_cv), pred_idx)

      KtrTr_raw <- W_mat[train_idx, train_idx, drop = FALSE]
      KteTr_raw <- W_mat[pred_idx,  train_idx, drop = FALSE]

      KtrTr_c <- .cv_center_KtrTr(KtrTr_raw)
      KteTr_c <- .cv_center_KteTr(KteTr_raw, KtrTr_raw)

      Y_tr_mean <- colMeans(y[train_idx, , drop = FALSE])
      Y_train_c <- sweep(y[train_idx, , drop = FALSE], 2, Y_tr_mean)
      Y_test_c  <- sweep(y[pred_idx,  , drop = FALSE], 2, Y_tr_mean)

      # Fit a separate KOPLS model for each ort level (nox = 0, 1, ..., nort).
      # Using the same model for all levels would give wrong "p" Q2 because
      # Bt is computed from the nort-deflated Tp; fitting separately ensures
      # each level's Bt comes from its own ort-deflated kernel.
      for (nox in 0L:nort) {
        .lev  <- nox + 1L
        m_lev <- .cv_kopls_fit(KtrTr_c, Y_train_c, ndim, nox)
        Yhat  <- .cv_kopls_predict(KteTr_c, m_lev, nox)
        Ypred_cv_all[[.lev]][pred_idx, ] <- sweep(Yhat, 2, Y_tr_mean, "+")
        pressy_cls_lev[.lev, ] <- pressy_cls_lev[.lev, ] + colSums((Y_test_c - Yhat)^2)
      }

      TSS_cv_cls    <- TSS_cv_cls + colSums(Y_test_c^2)
      cv_test_order <- c(cv_test_order, pred_idx)
    }

    # Per-class Q2 from the final (full) model level — named by class, like ComDim_PLS
    Q2_cv <- setNames(1 - pressy_cls_lev[n_lev, ] / TSS_cv_cls, colnames(y))

    # Per ort-level overall Q2 kept in cv$Q2.levels for diagnostics
    q2_names  <- c("p", paste0("po", seq_len(nort)))
    Q2_levels <- setNames(1 - rowSums(pressy_cls_lev) / sum(TSS_cv_cls), q2_names)

    # DQ2 per class from the final model level
    if (method == "OPLS-DA") {
      DQ2_cv      <- setNames(rep(NA_real_, ncol(y)), colnames(y))
      Ypred_final <- Ypred_cv_all[[n_lev]]
      for (i in seq_len(ncol(y))) {
        y_obs  <- y[cv_test_order, i]
        y_pred <- Ypred_final[cv_test_order, i]
        E0     <- y_pred[y_obs == 0]
        E1     <- y_pred[y_obs == 1] - 1
        PRESSD <- sum(E0[E0 > 0]^2) + sum(E1[E1 < 0]^2)
        TSS_i  <- sum((y_obs - mean(y_obs))^2)
        DQ2_cv[i] <- 1 - PRESSD / TSS_i
      }
    } else {
      DQ2_cv <- vector()
    }

    Q2 <- Q2_cv
    if (method == "OPLS-DA") DQ2 <- DQ2_cv

    # Convenience: expose the full-model and pred-only Ypred for compatibility
    Ypred_cv   <- Ypred_cv_all[[n_lev]]   # full model (po{nort})
    Ypred_cv_p <- Ypred_cv_all[[1L]]      # pred-only (p)

    cv_result <- list(
      k          = cv.k,
      Ypred      = Ypred_cv,    # full-model CV predictions (po{nort} level)
      Ypred.p    = Ypred_cv_p,  # pred-only CV predictions (p level)
      Q2         = Q2_cv,
      Q2.levels  = Q2_levels,   # overall Q2 per ort level for diagnostics
      DQ2        = if (method == "OPLS-DA") DQ2_cv else NULL
    )
  } else {
    if (cv.k != 0) {
      warning("'cv.k' must be >= 2 or 0 (skip CV). No cross-validation performed.")
    }
    cv_result <- list()
  }
  ## -------------------------------------------------------------------------

  end_output <- Sys.time()
  running_time <- (end_output - start_ini_time) # Total time of analysis
  if (loquace) {
    print(sprintf("Analysis finished after : %s millisecs", running_time * 1000))
  }

  progress_bar <- utils::txtProgressBar(min = 0, max = 80, style = 3, char = "=")
  utils::setTxtProgressBar(progress_bar, value = 80)

  close(progress_bar)

  # Save data in a ComDim structure

  res_calib <- new("ComDim",
    Method = method, # OPLS-DA or OPLS-R
    ndim = ndim,
    Q.scores = res_calib$Q,
    T.scores = res_calib$T_Loc,
    P.loadings = res_calib$P,
    Saliences = res_calib$saliences,
    Orthogonal = list(
      nort = nort,
      Q.scores = Q_ort,
      T.scores = T_Loc_ort,
      P.loadings.ort = orthoP,
      Saliences.ort = saliences_ort
    ),
    VIP = VIP_global,
    VIP.block = VIP_blocks,
    R2X = res_calib$explained,
    R2Y = r2y,
    Q2 = Q2,
    DQ2 = if (method == "OPLS-DA") {
      DQ2
    } else {
      vector()
    },
    Singular = res_calib$SingVal,
    Mean = list(
      MeanMB = res_calib$MeanMB,
      MeanY = meanY
    ),
    Norm = list(
      NormMB = res_calib$NormMB,
      FrobNorms = frob_norms,
      RVweights = rv_weights
    ),
    PLS.model = list(
      W = W_pred,
      B = B_kopls,
      B0 = B0_kopls,
      Y = y
    ),
    cv = cv_result,
    Prediction = if (method == "OPLS-DA") {
      list(
        Y.pred = Ypred,
        decisionRule = decisionRule,
        trueClass = classVect,
        predClass = data.frame(predClass = predClass, row.names = MB@Samples),
        Sensitivity = sensvec,
        Specificity = specvec,
        confusionMatrix = confusionMatrix
      )
    } else {
      list(Y.pred = Ypred)
    },
    Metadata = res_calib$metadata,
    variable.block = belong_block,
    runtime = as.numeric(running_time)
  )

  return(res_calib)
}


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
