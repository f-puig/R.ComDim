#' ComDim_y - Extending PLS-like supervised methods to multi-block data.
#'
#' Extends any PLS-like method used for regression or discriminant purposes to
#' the multi-block field. The user provides a function (\code{FUN}) that
#' computes one predictive component from the salience-weighted concatenated
#' blocks; global scores, local scores, and loadings are then derived following
#' the traditional ComDim-PLS framework. Optionally, orthogonal components
#' returned by FUN (e.g. from an O-PLS wrapper) are captured. VIP scores and
#' k-fold cross-validation are also supported.
#'
#' @param MB A MultiBlock object.
#' @param y The response: a numeric vector or matrix for regression
#'   (\code{type = 'regression'}), or a class vector or dummy matrix for
#'   discriminant analysis (\code{type = 'discriminant'}).
#' @param ndim Number of predictive Common Dimensions. If \code{NULL},
#'   defaults to the number of blocks.
#' @param FUN The function used as the core of the ComDim analysis. It must
#'   accept \code{W} (salience-weighted concatenated blocks, n x p matrix),
#'   \code{y} (response), and \code{ndim} (number of components to compute)
#'   as its first three named arguments, and return a named list with at least:
#'   \describe{
#'     \item{scores}{X scores vector (length n) for the current component.}
#'     \item{P}{X loadings vector (length p) for the current component.}
#'     \item{W}{X weights vector (length p) for the current component.}
#'     \item{Q}{Y loadings vector for the current component.}
#'     \item{U}{Y scores vector (length n) for the current component.}
#'   }
#'   Optional return fields:
#'   \describe{
#'     \item{y}{y as used internally by FUN (e.g. after centring/scaling),
#'       so that \code{ComDim_y} can detect automatic y-transformation and
#'       back-transform the Q loadings for correct B and B0 computation.}
#'     \item{orthoscores}{A numeric vector (length n) or single-column matrix
#'       of orthogonal X-scores. \strong{Required} when \code{nort = 1}.
#'       Only the first column is used. FUN does not need to change its
#'       behaviour between ort and predictive phases; if it always returns
#'       \code{orthoscores}, the framework uses it only during the ort phase
#'       and ignores it during the predictive phase.}
#'   }
#' @param nort Number of orthogonal Common Dimensions to extract before the
#'   predictive loop. Default \code{0} (no orthogonalization, pure PLS-like).
#'   Only \code{nort = 0} or \code{nort = 1} are accepted; for multiple
#'   orthogonal components use \code{ComDim_OPLS()}.
#'   When \code{nort = 1} the function follows the same two-phase architecture
#'   as \code{ComDim_OPLS}: the orthogonal component is removed from the
#'   (lambda-weighted) concatenated blocks first, then predictive components are
#'   computed from the orthogonally-deflated blocks. \code{FUN} must return
#'   \code{output$orthoscores} (a numeric vector of length n, or a matrix whose
#'   first column is used) when called with \code{ndim = 1}.
#'   Cross-validation is automatically skipped when \code{nort > 0}.
#' @param type \code{'regression'} (default) or \code{'discriminant'}.
#' @param decisionRule Only used when \code{type = 'discriminant'}. If
#'   \code{'fixed'}, samples are assigned to the class whose predicted score
#'   exceeds \code{1/nclasses}; if \code{'max'}, the class with the highest
#'   predicted score is chosen. Default \code{'max'}.
#' @param normalise To apply block normalisation. \code{FALSE} == no
#'   (default), \code{TRUE} == yes.
#' @param scale.y Logical (default \code{FALSE}). When \code{TRUE} and
#'   \code{type = 'regression'}, each column of Y is mean-centred and scaled to
#'   unit variance before being passed to \code{FUN}. The stored B, B0, and
#'   Ypred are back-transformed to the original Y scale, so downstream outputs
#'   (R2Y, Q2, predictions) are always in the original units. Ignored when
#'   \code{type = 'discriminant'} (dummy Y should not be scaled).
#' @param threshold Convergence threshold: iterations stop when the change
#'   in the global score vector falls below this value (default \code{1e-10}).
#' @param loquace Display computation time at each step. \code{TRUE} == yes,
#'   \code{FALSE} == no (default).
#' @param method A string label identifying the method (default: \code{'FUN'}).
#' @param cv.k Number of folds for k-fold cross-validation (default 7). Set
#'   to 0 to skip CV. When \code{cv.k >= 2}, the \code{Q2} and \code{DQ2}
#'   slots in the output reflect cross-validated predictive ability; otherwise
#'   they reflect training-set fit. CV is skipped when \code{nort > 0}.
#' @param ... Additional arguments passed to \code{FUN}.
#' @return A \code{ComDim} object with the following slots:
#' \describe{
#'   \item{\code{Method}}{The label supplied via the \code{method} argument.}
#'   \item{\code{ndim}}{Number of predictive Common Dimensions extracted.}
#'   \item{\code{Q.scores}}{Global consensus scores matrix (\eqn{n \times
#'     ndim}).  Each column \eqn{\mathbf{q}_a} (unit-norm) is derived from
#'     the dominant left direction of FUN applied to the salience-weighted
#'     concatenated blocks.}
#'   \item{\code{T.scores}}{Named list of block-specific local scores
#'     (\eqn{n \times ndim} each).  Local loading
#'     \eqn{\mathbf{p}_{ba} = \tilde{\mathbf{X}}_b'\mathbf{q}_a} (computed
#'     on the ort-deflated block when \code{nort > 0}); local score
#'     \eqn{\mathbf{t}_{ba} =
#'     \tilde{\mathbf{X}}_b\,\mathbf{p}_{ba}(\mathbf{p}_{ba}'\mathbf{p}_{ba})^{-1}}.}
#'   \item{\code{P.loadings}}{Global loadings (\eqn{p_{tot} \times ndim}):
#'     \eqn{\mathbf{P} = \tilde{\mathbf{X}}'\mathbf{Q}}, where
#'     \eqn{\tilde{\mathbf{X}}} is the (optionally ort-deflated) mean-centred
#'     concatenated blocks.}
#'   \item{\code{Saliences}}{Block salience matrix (\eqn{ntable \times ndim}):
#'     \eqn{\lambda_{ba} =
#'     \mathbf{q}_a'\tilde{\mathbf{X}}_b\tilde{\mathbf{X}}_b'\mathbf{q}_a}.}
#'   \item{\code{R2X}}{Proportion of X variance captured by each predictive
#'     component (named vector, length \eqn{ndim}).  Let \eqn{\mathbf{t}_a}
#'     be the X-score vector returned by FUN for component \eqn{a}:
#'     \deqn{R2X_a = \|\mathbf{t}_a\|^4 \big/ \sum_k \|\mathbf{t}_k\|^4.}
#'     When \code{nort > 0}, the denominator also includes the orthogonal
#'     \eqn{\|\mathbf{t}_{ort,k}\|^4} terms, and the orthogonal R2X fractions
#'     are stored separately in \code{Orthogonal$R2X}.}
#'   \item{\code{R2Y}}{Cumulative Y-variance explained (named vector, length
#'     \eqn{ndim}):
#'     \deqn{R2Y_a = 1 - RSS_a / TSS_Y,}
#'     where \eqn{RSS_a} is the residual SS from an OLS regression of
#'     \eqn{\mathbf{Y}} on \eqn{[1, \mathbf{q}_1, \ldots, \mathbf{q}_a]}.
#'     \strong{Note:} \eqn{R2Y_a} is cumulative -- the total Y-variance
#'     explained by the first \eqn{a} components together, not the marginal
#'     contribution of component \eqn{a} alone.}
#'   \item{\code{Q2}}{Predictive Q2 per response column (regression) or per
#'     class (discriminant), named accordingly:
#'     \deqn{Q2 = 1 - PRESS / TSS_Y,}
#'     where \eqn{PRESS = \sum_i (\hat{y}_i - y_i)^2}.
#'     When \code{cv.k >= 2} and \code{nort = 0}: cross-validated (out-of-
#'     sample) predictions are used; otherwise training-set predictions.
#'     CV is automatically skipped when \code{nort > 0}.}
#'   \item{\code{DQ2}}{(Discriminant mode only) Discriminant Q2 per class,
#'     using only penalising residuals:
#'     \deqn{DQ2 = 1 - PRESSD / TSS_Y,}
#'     where \eqn{PRESSD} sums \eqn{\hat{y}_i^2} for class-0 samples with
#'     \eqn{\hat{y}_i > 0}, and \eqn{(\hat{y}_i - 1)^2} for class-1 samples
#'     with \eqn{\hat{y}_i < 1}.  Same cross-validation logic as \code{Q2}.}
#'   \item{\code{Singular}}{Squared L2 norm of the FUN X-score vector per
#'     component (\eqn{\|\mathbf{t}_a\|^2}), used to derive \code{R2X}.}
#'   \item{\code{VIP}}{Global total VIP (named vector, length \eqn{p_{tot}}):
#'     concatenation of \code{VIP.block[[b]]$tot} across blocks.  When
#'     \code{nort = 0}, uses the Wold formula; when \code{nort = 1}, tot
#'     combines predictive and orthogonal VIPs (see \code{VIP.block}).}
#'   \item{\code{VIP.block}}{Named list (one \code{data.frame} per block).
#'     When \code{nort = 0}: columns \code{p} and \code{tot} (= \code{p}),
#'     using the Wold formula:
#'     \deqn{VIPp_j = \sqrt{p_b \cdot
#'       \frac{\sum_a s_a \tilde{w}_{j,a}^2}{\sum_a s_a}},}
#'     where \eqn{s_a = \|\mathbf{t}_a\|^2\|\mathbf{q}_a\|^2} and
#'     \eqn{\tilde{w}_{j,a} = w_{j,a}/\|\mathbf{w}_a\|} is the L2-normalised
#'     \eqn{j}-th element of the \eqn{a}-th weight vector.
#'     When \code{nort = 1}: columns \code{p} (Wold, same as above),
#'     \code{o} (orthogonal VIP, loadings-based:
#'     \eqn{VIPo_j = \sqrt{p_b \cdot \sum_a s_{oa}\tilde{P}_{o,j,a}^2 /
#'     \sum_a s_{oa}}},
#'     where \eqn{s_{oa} = \|\mathbf{q}_{ort}[,a]\|^2} and
#'     \eqn{\tilde{\mathbf{P}}_o} is the column-L2-normalised block-slice of
#'     the ort loadings), and \code{tot}
#'     (\eqn{VIPtot_j = \sqrt{(VIPp_j^2 + VIPo_j^2)/2}}).
#'     Row names are variable names.}
#'   \item{\code{PLS.model}}{List with: \code{W} (X weight matrix collected
#'     from FUN, \eqn{p_{tot} \times ndim}); \code{B} (regression
#'     coefficients,
#'     \eqn{\mathbf{B} = \mathbf{W}(\mathbf{P}'\mathbf{W})^{-1}\mathbf{Q}'},
#'     in original Y units); \code{B0} (intercept,
#'     \eqn{\mathbf{B}_0 = \bar{\mathbf{y}} -
#'     \overline{\tilde{\mathbf{x}}}\mathbf{B}}); \code{Y} (original
#'     response matrix as supplied).
#'     Training-set predictions:
#'     \eqn{\hat{\mathbf{Y}} =
#'     \tilde{\mathbf{X}}\mathbf{B} + \mathbf{B}_0}.}
#'   \item{\code{cv}}{Cross-validation results when \code{cv.k >= 2} and
#'     \code{nort = 0} (empty list otherwise): \code{k}, \code{fold}
#'     (sample-to-fold vector), \code{Ypred} (\eqn{n \times ncol(Y)}
#'     out-of-sample predictions), \code{Q2} (CV Q2 per class/response),
#'     \code{DQ2} (mean CV DQ2, discriminant only),
#'     \code{DQ2.perclass} (CV DQ2 per class, discriminant only).}
#'   \item{\code{Orthogonal}}{When \code{nort > 0}: list with \code{nort},
#'     \code{Q.scores} (global ort scores, \eqn{n \times nort}, unit-norm),
#'     \code{T.scores} (block ort local scores, \eqn{n \times nort} each),
#'     \code{P.loadings.ort} (ort loadings, \eqn{p_{tot} \times nort}),
#'     \code{Saliences.ort} (\eqn{ntable \times nort}), and \code{R2X}
#'     (orthogonal X-variance fractions,
#'     \eqn{R2X_{ort,a} = \|\mathbf{t}_{ort,a}\|^4 / total}).
#'     Empty list when \code{nort = 0}.}
#'   \item{\code{Prediction}}{Training-set predictions: \code{Y.pred}
#'     (\eqn{n \times ncol(Y)}); for discriminant analysis also
#'     \code{decisionRule}, \code{trueClass}, \code{predClass} (data.frame),
#'     \code{Sensitivity} and \code{Specificity} (per class),
#'     \code{confusionMatrix} (named list of 2x2 matrices).}
#'   \item{\code{Mean}}{List with \code{MeanMB} (column means per block),
#'     \code{MeanY} (column means of Y), and \code{ScaleY} (column SDs of Y;
#'     all ones when \code{scale.y = FALSE}).}
#'   \item{\code{Norm}}{List with \code{NormMB}: Frobenius norms for block
#'     normalisation.}
#'   \item{\code{variable.block}}{Character vector (length \eqn{p_{tot}})
#'     mapping each row of \code{P.loadings} and each element of \code{VIP}
#'     to its block.}
#'   \item{\code{runtime}}{Total computation time in seconds.}
#' }
#' @examples
#' b1 <- matrix(rnorm(500), 10, 50) # 10 samples, 50 variables
#' b2 <- matrix(rnorm(800), 10, 80) # 10 samples, 80 variables
#' mb <- MultiBlock(Data = list(b1 = b1, b2 = b2))
#'
#' ## Example 1: ComDim-PLS (regression) ---------------------------------------
#' # Single-step NIPALS PLS wrapper (one predictive component per call).
#' # Note: 'tx' is used instead of 't' to avoid shadowing base::t().
#' fun.PLS <- function(W, y, ndim, ...) {
#'   output <- list()
#'   w <- t(W) %*% y / as.numeric(t(y) %*% y) # X weight (u = y, 1 step)
#'   w <- w / sqrt(sum(w^2)) # L2 normalise
#'   tx <- W %*% w # X score
#'   p <- t(W) %*% tx / as.numeric(t(tx) %*% tx) # X loading
#'   q <- t(y) %*% tx / as.numeric(t(tx) %*% tx) # Y loading
#'   u <- y %*% q / as.numeric(t(q) %*% q) # Y score
#'   output$scores <- as.vector(tx)
#'   output$P <- as.vector(p)
#'   output$W <- as.vector(w)
#'   output$Q <- as.vector(q)
#'   output$U <- as.vector(u)
#'   return(output)
#' }
#'
#' y <- c(1, 1, 1, 1, 1, 5, 5, 5, 10, 10)
#' resultsPLS <- ComDim_y(mb,
#'   y = y, ndim = 2,
#'   type = "regression",
#'   FUN = fun.PLS,
#'   method = "PLS",
#'   cv.k = 0
#' )
#'
#' ## Example 2: ComDim-OPLS-DA (discriminant, nort = 1) ----------------------
#' # Thin wrapper around OPLS_NIPALS_DNR(), the package's NIPALS OPLS engine.
#' # All inputs (W, y, and any extra args such as 'threshold') are forwarded
#' # directly via '...'. Use this pattern when nort > 0; for nort = 0 the
#' # simpler PLS wrapper in Example 1 is sufficient (no orthoscores needed).
#' fun.OPLS <- function(W, y, ndim, ...) {
#'   res <- OPLS_NIPALS_DNR(W = W, y = y, ...)
#'   list(
#'     scores      = as.vector(res$t_pred),
#'     P           = as.vector(res$p),
#'     W           = as.vector(res$w_pred),
#'     Q           = as.vector(res$q),
#'     U           = as.vector(res$u),
#'     orthoscores = matrix(res$t_ort, ncol = 1)
#'   )
#' }
#'
#' groups <- c(rep("A", 5), rep("B", 5))
#' resultsOPLS <- ComDim_y(mb,
#'   y = groups, ndim = 1,
#'   nort = 1,
#'   type = "discriminant",
#'   FUN = fun.OPLS,
#'   method = "OPLS-DA",
#'   cv.k = 0
#' )
#'
#' ## Example 3 (not run): ComDim-OPLS-DA via ropls ---------------------------
#' # Wrapping ropls::opls is also possible. Key points:
#' #   - Use orthoI = 1 (fixed) instead of NA so the output is predictable.
#' #   - Always return output$orthoscores; ComDim_y ignores it in phases
#' #     where ort has already been removed.
#' #   - Expand the single ropls Q loading to match the ncol(y_dummy) width.
#' \donttest{
#' if (requireNamespace("ropls", quietly = TRUE)) {
#' fun.OPLSDA.ropls <- function(W, y, ndim, ...) {
#'   output <- list()
#'   # Convert dummy matrix to ropls-compatible -1/+1 vector
#'   Y <- c(-1, 1)[apply(y, 1, function(x) match(1, x))]
#'   result <- tryCatch(
#'     ropls::opls(
#'       x = W, y = Y, predI = 1, orthoI = 1,
#'       fig.pdfC = "none", info.txtC = "none"
#'     ),
#'     error = function(e) {
#'       ropls::opls(
#'         x = W, y = Y, predI = 1, orthoI = 0,
#'         fig.pdfC = "none", info.txtC = "none"
#'       )
#'     }
#'   )
#'   output$scores <- result@scoreMN[, 1]
#'   output$P <- result@loadingMN[, 1]
#'   output$W <- result@weightMN[, 1]
#'   output$U <- result@uMN[, 1]
#'   # Expand the single ropls Q loading to match the 2-column dummy matrix:
#'   # loadings for class1 and class2 are antisymmetric in binary PLS-DA.
#'   output$Q <- c(-result@cMN[, 1], result@cMN[, 1])
#'   output$y <- result@suppLs$yModelMN # internal y (for scaling detection)
#'   # Orthogonal scores (used during the ort pre-loop when nort > 0)
#'   if (!is.null(result@orthoScoreMN) && ncol(result@orthoScoreMN) > 0) {
#'     output$orthoscores <- result@orthoScoreMN   # n x k matrix; col jj used for jj-th ort
#'   } else {
#'     output$orthoscores <- matrix(0, nrow = nrow(W), ncol = 1)
#'   }
#'   return(output)
#' }
#'
#' b1_r <- matrix(rnorm(8 * 30), 8, 30)
#' b2_r <- matrix(rnorm(8 * 20), 8, 20)
#' mb_r <- MultiBlock(Data = list(b1 = b1_r, b2 = b2_r))
#' resultsOPLSDA <- ComDim_y(mb_r,
#'   y = c(rep("NI", 4), rep("OFF", 4)),
#'   ndim = 1,
#'   nort = 1,
#'   type = "discriminant",
#'   FUN = fun.OPLSDA.ropls,
#'   method = "OPLS-DA(ropls)",
#'   cv.k = 0
#' )
#' }
#' }
#' @export

ComDim_y <- function(MB = MB, y = y,
                     ndim = NULL, FUN = FUN,
                     nort = 0L,
                     type = c("regression", "discriminant")[1],
                     decisionRule = c("fixed", "max")[2],
                     normalise = FALSE, scale.y = FALSE,
                     threshold = 1e-10,
                     loquace = FALSE,
                     method = "FUN",
                     cv.k = 7,
                     ...) {
  # INITIALISATION

  progress_bar <- utils::txtProgressBar(min = 0, max = 80, style = 3, char = "=")

  start_ini_time <- Sys.time() # To start counting calculation time for the initialization


  if (!inherits(MB, "MultiBlock")) {
    stop("'MB' is  not a MultiBlock.")
  }

  if (!(type %in% c("regression", "discriminant"))) {
    stop("'type' must be either 'regression' or 'discriminant'.")
  }

  if (type == "discriminant" && !(decisionRule %in% c("fixed", "max"))) {
    stop("'decisionRule' must be either 'fixed' or 'max'.")
  }

  nort <- as.integer(nort)
  if (nort < 0L) {
    stop("'nort' must be a non-negative integer.")
  }
  if (nort > 1L) {
    stop("'nort > 1' is not supported in ComDim_y(). For multiple orthogonal components, use ComDim_OPLS().")
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
      if (numblock0[[i]] %in% MB@Batch) {
        MB@Batch[[numblock0[i]]] <- NULL
      }
      if (numblock0[[i]] %in% MB@Metadata) {
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


  ## Check y matrix (convert to dummy matrix if it is not already.)
  if (type == "discriminant") {
    tmp <- unique(as.vector(y))
    if (all(tmp %in% c(0, 1))) { # Is it dummy?

      if (!is.matrix(y)) {
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
        stop("Predictive analysis can only be performed if 'y' is a dummy matrix or a class vector.")
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

  if (!is.matrix(y)) { # If it's a vector and the method is PLS-R
    y <- as.matrix(y)
  }

  # Store original y (after validation / dummy conversion) for use in R2Y, Q2,
  # and as the response saved in PLS.model$Y.
  y_raw <- y

  # Scale Y for regression if requested.
  # For discriminant analysis the dummy matrix must not be scaled; emit a
  # warning and skip.
  if (scale.y && type == "discriminant") {
    warning("'scale.y = TRUE' is ignored when type = 'discriminant': dummy Y matrices should not be centred or scaled.")
    scale.y <- FALSE
  }
  if (scale.y) {
    y_center <- colMeans(y)
    y_sd <- apply(y, 2, sd)
    y_sd[y_sd < .Machine$double.eps] <- 1 # guard against constant columns
    y <- sweep(sweep(y, 2, y_center, "-"), 2, y_sd, "/")
  } else {
    y_center <- colMeans(y)
    y_sd <- rep(1, ncol(y))
  }

  if (give_error) {
    stop("The data is not ready for ComDim.")
  } else {
    message("The data can be used for ComDim.")
  }


  if (is.null(ndim)) {
    ndim <- ntable # If the number of components is not defined,
    # the number of components to extract is equal to the number of blocks.
  }

  pieceBar <- 4 + 2 * ntable + ndim + nort # Number of updates in the progress bar.
  pieceBar <- 80 / pieceBar
  total_progress <- pieceBar

  DimLabels <- paste0("CC", 1:ndim) # One label per predictive component.
  OrtLabels <- if (nort > 0) paste0("ort", 1:nort) else character(0)
  TableLabels <- blockNames(MB) # One label per block.

  end_ini_time <- Sys.time() # To end the count of the calculation time.

  if (loquace) {
    message(sprintf("Initialisation finished after : %s millisecs", (end_ini_time - start_ini_time) * 1000))
  }

  utils::setTxtProgressBar(progress_bar, value = total_progress)

  # NORMALISATION

  X_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))
  Xnorm_mat <- matrix(, nrow = nrowMB, ncol = sum(variable_number))

  res_calib <- list()
  temp_tabCalib <- list()
  s_r <- list()
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
      s_r[[i]] <- X_Normed
    } else {
      res_calib$NormMB[[TableLabels[i]]] <- rep(1, length(MB@Variables[[i]]))
      names(res_calib$NormMB[[TableLabels[i]]]) <- MB@Variables[[i]]
      temp_tabCalib[[i]] <- MB@Data[[i]]
      s_r[[i]] <- MB@Data[[i]]
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

    # Compute means (used for B0 and Y-scaling check)
    meanY <- colMeans(y)

    norm_comdim <- Sys.time()
    if (loquace) {
      message(sprintf("Normalization of block %s finished after : %s millisecs", i, (norm_comdim - start_ini_time) * 1000))
    }
    total_progress <- total_progress + pieceBar
    utils::setTxtProgressBar(progress_bar, value = total_progress)
  }

  names(res_calib$NormMB) <- TableLabels

  nR <- nrow(Xnorm_mat)
  nC <- ncol(Xnorm_mat)

  total_progress <- total_progress + pieceBar * ntable
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  # ICs decrements from ndim to 1 across outer loop (consistent with ComDim_PLS)
  ICs <- ndim

  ## PRE-LOOP: orthogonal component extraction
  ##
  ## When nort > 0, extract nort orthogonal components from the
  ## (lambda-weighted) concatenated blocks using FUN. Each block is deflated by
  ## its orthogonal consensus score before the predictive loop runs. This two-
  ## phase structure is identical to ComDim_OPLS. FUN must return
  ## output$orthoscores (an n x 1 matrix or vector) when called with ndim = 1.

  Q_ort <- NULL
  saliences_ort <- NULL
  varexp_ort <- NULL

  if (nort > 0) {
    Q_ort <- matrix(, ncol = nort, nrow = nrowMB)
    saliences_ort <- matrix(, ncol = nort, nrow = ntable)
    varexp_ort <- as.vector(NULL)

    for (jj in 1:nort) {
      lambda_ort <- matrix(rep(1, ntable), nrow = ntable, ncol = 1)
      qini_ort <- as.vector(s_r[[1]][, 1])
      qini_ort <- qini_ort / sqrt(as.vector(qini_ort %*% qini_ort))
      qold_ort <- 100

      iters <- 0
      while (norm(qold_ort - qini_ort, "2") > threshold && iters < 100) {
        iters <- iters + 1
        qold_ort <- qini_ort
        q_ort_sum <- 0

        W <- matrix(, nrow = nrowMB, ncol = 0)
        for (j in 1:ntable) {
          W <- cbind(W, sqrt(lambda_ort[j]) * s_r[[j]])
        }

        output_ort <- FUN(W = W, y = y, ndim = 1, ...)

        if (is.null(output_ort$orthoscores)) {
          stop(
            "'FUN' must return 'orthoscores' (a numeric matrix n x k, k >= 1, ",
            "or a length-n vector) when 'nort' > 0."
          )
        }

        ort_mat <- as.matrix(output_ort$orthoscores)
        if (ncol(ort_mat) < 1L) {
          stop("'FUN$orthoscores' must have at least 1 column when 'nort' > 0.")
        }
        ort_col <- 1L
        ort_vec <- ort_mat[, ort_col]
        sv_ort  <- sqrt(sum(ort_vec^2))
        qtmp_ort <- ort_vec / sv_ort

        for (j in 1:ntable) {
          lambda_ort[j] <- t(qtmp_ort) %*% (s_r[[j]] %*% t(s_r[[j]])) %*% qtmp_ort
          q_ort_sum <- q_ort_sum + lambda_ort[j] * qtmp_ort
        }

        q_ort_sum <- q_ort_sum / sqrt(as.vector(t(q_ort_sum) %*% q_ort_sum))
        if (abs(min(q_ort_sum)) > abs(max(q_ort_sum))) q_ort_sum <- -q_ort_sum
        qini_ort <- q_ort_sum
      } # convergence loop

      saliences_ort[, jj] <- lambda_ort
      Q_ort[, jj] <- qini_ort
      varexp_ort[jj] <- sv_ort^4

      # Deflate each block by the orthogonal consensus score (same as OPLS_MB)
      for (j in 1:ntable) {
        p_j_ort <- t(s_r[[j]]) %*% Q_ort[, jj] /
          as.vector(t(Q_ort[, jj]) %*% Q_ort[, jj])
        s_r[[j]] <- s_r[[j]] - Q_ort[, jj] %*% t(p_j_ort)
      }

      ort_comdim <- Sys.time()
      if (loquace) {
        message(sprintf(
          "Orthogonal component %s determined after : %s millisecs",
          jj, (ort_comdim - start_ini_time) * 1000
        ))
      }
      total_progress <- total_progress + pieceBar
      utils::setTxtProgressBar(progress_bar, value = total_progress)
    } # jj loop

    colnames(Q_ort) <- OrtLabels
    rownames(Q_ort) <- MB@Samples
    colnames(saliences_ort) <- OrtLabels
    rownames(saliences_ort) <- TableLabels
  } # nort > 0

  ## PREDICTIVE EXTRACTION
  ## Run ComDim-PLS on the (orthogonally-deflated) s_r blocks.

  saliences <- matrix(, ncol = ndim, nrow = ntable)
  Q <- matrix(, ncol = ndim, nrow = nrowMB)
  varexp <- as.vector(NULL)
  PLS2_Scores <- matrix(, ncol = ndim, nrow = nrowMB)
  PLS2_P <- matrix(, ncol = ndim, nrow = sum(variable_number))
  PLS2_W <- matrix(, ncol = ndim, nrow = sum(variable_number))
  PLS2_Q <- matrix(, ncol = ndim, nrow = ncol(as.matrix(y)))
  PLS2_U <- matrix(, ncol = ndim, nrow = nrowMB)

  for (dims in 1:ndim) {
    lambda <- matrix(rep(1, ntable), nrow = ntable, ncol = 1)

    qini <- as.vector(s_r[[1]][, 1])
    qini <- qini / sqrt(as.vector(qini %*% qini)) # random initialisation + unit length
    qold <- 100

    iters <- 0
    while (norm(qold - qini, "2") > threshold && iters < 100) {
      iters <- iters + 1
      qold <- qini
      q <- 0

      W <- matrix(, nrow = nrowMB, ncol = 0)
      for (j in 1:ntable) {
        W <- cbind(W, sqrt(lambda[j]) * s_r[[j]])
      }

      # HERE IS THE CORE COMDIM
      # FUN receives ICs (current remaining components), consistent with
      # ComDim_PLS approach.
      output <- FUN(W = W, y = y, ndim = ICs, ...)

      PLS2_Scores[, dims] <- output$scores
      PLS2_P[, dims] <- output$P
      PLS2_W[, dims] <- output$W
      PLS2_Q[, dims] <- output$Q
      PLS2_U[, dims] <- output$U

      sv <- sqrt(t(PLS2_Scores[, dims]) %*% PLS2_Scores[, dims]) # For R2X
      qtmp <- PLS2_Scores[, dims] / as.vector(sv)

      for (j in 1:ntable) {
        lambda[j] <- t(qtmp) %*% (s_r[[j]] %*% t(s_r[[j]])) %*% qtmp
        q <- q + lambda[j] * qtmp
      }

      q <- q / sqrt(as.vector(t(q) %*% q)) # standardise
      if (abs(min(q)) > abs(max(q))) q <- -q
      qini <- q
    } # convergence loop

    saliences[, dims] <- lambda
    Q[, dims] <- q

    res_calib$SingVal[dims] <- sv^2 # Calculated from the scores
    varexp[dims] <- res_calib$SingVal[dims]^2

    # Deflate blocks by predictive consensus score
    aux <- diag(nrowMB) - as.matrix(q) %*% t(as.matrix(q))
    for (j in 1:ntable) {
      s_r[[j]] <- aux %*% s_r[[j]]
    }

    # Decrement ICs for the next component (consistent with ComDim_PLS)
    ICs <- ICs - 1

    iter_comdim <- Sys.time()
    if (loquace) {
      message(sprintf(
        "Component %s determined after : %s millisecs",
        dims, (iter_comdim - start_ini_time) * 1000
      ))
    }
    total_progress <- total_progress + pieceBar
    utils::setTxtProgressBar(progress_bar, value = total_progress)
  } # dims loop

  # Strip dim names from the internally-transformed y (if FUN returned it)
  if (!is.null(output$y)) {
    rownames(output$y) <- NULL
    colnames(output$y) <- NULL
  }

  names(res_calib$SingVal) <- DimLabels

  colnames(Q) <- DimLabels
  rownames(Q) <- MB@Samples
  res_calib$Q <- Q
  rm(Q)

  ## Adding metadata. Metadata extracted from the first block
  res_calib$metadata <- list()
  if (length(MB@Metadata) != 0) {
    res_calib$metadata[[1]] <- MB@Metadata
  }


  end_comdim <- Sys.time()
  if (loquace) {
    message(sprintf("Scores finished after : %s millisecs", (end_comdim - start_ini_time) * 1000))
  }
  total_progress <- total_progress + pieceBar
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  rm(s_r)

  ##
  explained <- varexp / sum(varexp)
  explained.ort <- NULL
  if (nort > 0 && !is.null(varexp_ort)) {
    total.varexp <- sum(varexp) + sum(varexp_ort)
    explained <- varexp / total.varexp
    explained.ort <- varexp_ort / total.varexp
    names(explained.ort) <- OrtLabels
  }
  names(explained) <- DimLabels
  res_calib$explained <- explained
  rm(varexp)
  rm(explained)

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

  ## Calculate Normalised concatenated Xs ('Calib') from col
  # Calculate concatenated CD loadings
  # Reorganise Loadings - 1 matrix / LV

  L_CD_Vec <- NULL
  L_X_Vec <- NULL
  Calib <- NULL
  nCalib <- nrowMB

  # Calculate concatenated CD loadings / Local scores
  L_X <- list()
  T_Loc_ort <- T_Loc <- list()
  for (i in 1:ntable) { # Prepare lists
    T_Loc[[TableLabels[i]]] <- matrix(, ncol = ndim, nrow = nrowMB)
    if (nort > 0) {
      T_Loc_ort[[TableLabels[i]]] <- matrix(, ncol = nort, nrow = nrowMB)
    }
  }
  b <- matrix(, ncol = ndim, nrow = ntable)

  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("The package pracma is needed.")
  }

  # Orthogonal loadings: deflate temp_tabCalib block by block, compute local
  # orthogonal scores and the concatenated orthogonal loading matrix (orthoP).
  # Rebuild Xnorm_mat from the orthogonally-deflated temp_tabCalib so that
  # the global loadings and B0 are computed in the correct subspace.
  orthoP <- NULL
  if (nort > 0) {
    for (j in 1:nort) {
      for (i in 1:ntable) {
        tempo <- t(temp_tabCalib[[i]]) %*% Q_ort[, j] # local ort loadings
        T_Loc_ort[[TableLabels[i]]][, j] <-
          temp_tabCalib[[i]] %*% (tempo %*% pracma::pinv(t(tempo) %*% tempo))
        # Deflate each temp_tabCalib by the orthogonal component
        temp_tabCalib[[i]] <- temp_tabCalib[[i]] - Q_ort[, j] %*% t(tempo)
        if (i == 1) {
          tempo2 <- tempo
        } else {
          tempo2 <- append(tempo2, tempo)
        }
      }
      if (j == 1) {
        orthoP <- tempo2
      } else {
        orthoP <- cbind(orthoP, tempo2)
      }
    }
    orthoP <- as.matrix(orthoP)
    colnames(orthoP) <- OrtLabels
    # Rebuild Xnorm_mat from orthogonally-deflated temp_tabCalib
    for (i in 1:ntable) {
      if (i == 1) {
        Xnorm_mat <- temp_tabCalib[[i]]
      } else {
        Xnorm_mat <- cbind(Xnorm_mat, temp_tabCalib[[i]])
      }
    }
  }

  # Predictive loadings: compute local scores and b-coefficients from the
  # (orthogonally-deflated) temp_tabCalib.
  for (j in 1:ndim) {
    T_mat <- matrix(, nrow = nrowMB, ncol = 0)

    for (i in 1:ntable) {
      temp <- t(temp_tabCalib[[i]]) %*% res_calib$Q[, j] # local loadings

      T_mat <- cbind(
        T_mat,
        temp_tabCalib[[i]] %*% (temp %*% pracma::pinv(t(temp) %*% temp))
      )

      T_Loc[[TableLabels[i]]][, j] <-
        temp_tabCalib[[i]] %*% (temp %*% pracma::pinv(t(temp) %*% temp))

      # Deflate each temp_tabCalib
      temp_tabCalib[[i]] <- temp_tabCalib[[i]] - res_calib$Q[, j] %*% t(temp)
    }

    # MLR b-coefficients between Local and Global Scores
    b[, j] <- pracma::pinv(t(T_mat) %*% T_mat) %*% t(T_mat) %*% res_calib$Q[, j]
  }

  for (i in 1:ntable) {
    rownames(T_Loc[[TableLabels[i]]]) <- rownames(res_calib$Q)
    colnames(T_Loc[[TableLabels[i]]]) <- colnames(res_calib$Q)
    if (nort > 0) {
      rownames(T_Loc_ort[[TableLabels[i]]]) <- MB@Samples
      colnames(T_Loc_ort[[TableLabels[i]]]) <- OrtLabels
    }
  }

  # Calculate Global Loadings
  L_CD_Vec <- t(Xnorm_mat) %*% res_calib$Q # Scaled CD 'global' Loadings

  load_comdim <- Sys.time()
  if (loquace) {
    message(sprintf("Loadings finished after : %s millisecs", (load_comdim - start_ini_time) * 1000))
  }

  total_progress <- total_progress + pieceBar
  utils::setTxtProgressBar(progress_bar, value = total_progress)

  rm(X_mat)
  rm(temp_tabCalib)
  rm(temp)

  ## Output
  end_function <- Sys.time()

  ## Regression coefficients B and intercept B0
  ##
  ## Two back-transformation paths:
  ## (a) scale.y = TRUE: the framework scaled y before passing it to FUN; B and
  ##     B0 are in scaled-y space and must be back-transformed to original units.
  ## (b) scale.y = FALSE but FUN internally scaled y (e.g. ropls): detected by
  ##     comparing output$y to standardised y and back-transforming Q loadings.
  ##     Falls back to un-transformed Q when output$y is NULL or comparison fails.
  if (scale.y) {
    ComDimPLS2_B <- PLS2_W %*% pracma::pinv(t(PLS2_P) %*% PLS2_W) %*% t(PLS2_Q)
  } else {
    y_was_scaled <- FALSE
    if (!is.null(output$y) && ncol(y) == 1) {
      y_was_scaled <- isTRUE(all.equal(
        as.vector(output$y[, 1]),
        as.vector((y[, 1] - mean(y[, 1])) / sd(y[, 1]))
      ))
    }
    if (y_was_scaled) {
      ComDimPLS2_B <- PLS2_W %*% pracma::pinv(t(PLS2_P) %*% PLS2_W) %*% t(PLS2_Q * sd(y) + mean(y))
    } else {
      ComDimPLS2_B <- PLS2_W %*% pracma::pinv(t(PLS2_P) %*% PLS2_W) %*% t(PLS2_Q)
    }
  }
  # B0: use colMeans(Xnorm_mat) here. When nort > 0, Xnorm_mat has been
  # overwritten with the orthogonally-deflated version (lines above), so
  # colMeans(Xnorm_mat) gives the ort-deflated means — consistent with Ypred
  # which is also computed from the deflated Xnorm_mat. When nort == 0,
  # Xnorm_mat is unchanged, so the result equals the original column means.
  ComDimPLS2_B0 <- colMeans(y) - colMeans(Xnorm_mat) %*% ComDimPLS2_B

  # Back-transform B and B0 to the original Y scale when scale.y = TRUE, then
  # restore y to original units for all downstream computations.
  if (scale.y) {
    ComDimPLS2_B <- sweep(ComDimPLS2_B, 2, y_sd, "*")
    ComDimPLS2_B0 <- ComDimPLS2_B0 * y_sd + y_center
    y <- y_raw # restore original y for R2Y, Q2, Ypred, etc.
  }

  ## R2Y — cumulative OLS of Y on [1 | Q_1..j] for j = 1..ndim.
  ## r2y[j] reports the fraction of Y variance explained by the first j global
  ## ComDim scores. An explicit intercept is included so that Y with non-zero
  ## mean is handled correctly.
  Y_c_r2   <- sweep(y, 2, colMeans(y))
  sstot_Y_r2 <- sum(Y_c_r2^2)
  if (sstot_Y_r2 < .Machine$double.eps) {
    r2y <- setNames(rep(0, ndim), paste0("CC", seq_len(ndim)))
  } else {
    r2y <- setNames(numeric(ndim), paste0("CC", seq_len(ndim)))
    for (.j in seq_len(ndim)) {
      T_r2_j     <- cbind(1, res_calib$Q[, seq_len(.j), drop = FALSE])
      B_r2_j     <- pracma::pinv(crossprod(T_r2_j)) %*% crossprod(T_r2_j, y)
      Ypred_r2_j <- T_r2_j %*% B_r2_j
      r2y[.j]    <- 1 - sum((Ypred_r2_j - y)^2) / sstot_Y_r2
    }
  }
  rm(Y_c_r2)

  res_calib$b <- b
  colnames(res_calib$b) <- DimLabels
  rownames(res_calib$b) <- TableLabels

  res_calib$T_Loc <- T_Loc

  # T_Loc is no longer needed
  rm(T_Loc)

  varnames <- vector()
  for (i in 1:ntable) {
    varnames <- append(varnames, MB@Variables[[i]])
  }

  res_calib$P <- L_CD_Vec
  colnames(res_calib$P) <- DimLabels
  rownames(res_calib$P) <- varnames
  rownames(ComDimPLS2_B) <- varnames
  if (!is.null(colnames(y_raw))) colnames(ComDimPLS2_B) <- colnames(y_raw)

  rm(L_CD_Vec)

  ## -------------------------------------------------------------------------
  ## VIP (Variable Importance in Projection)
  ##
  ## Predictive VIP ("p") — Wold formula using PLS2_W, PLS2_Scores, PLS2_Q:
  ##   W_j  = L2-normalised column of the weight matrix (per block)
  ##   s_a  = colSums(T^2) * colSums(Q^2)   [length ndim]
  ##   VIPp_j = sqrt( p_b * sum_a(W_ja^2 * s_a) / sum_a(s_a) )
  ##
  ## Orthogonal VIP ("o", only when nort > 0) — loadings-based formula:
  ##   VIPo_j = sqrt( p_b * sum(ss_to * P_norm_j^2) / sum(ss_to) )
  ## (ort components have no W/Q analogue, so their weight is measured via
  ##  the magnitude of the ort loadings in X-space.)
  ##
  ## Per-block output: data.frame with "p", "o" (if nort>0), "tot" columns.
  ## Global VIP (@VIP): "tot" concatenated across blocks.
  ## -------------------------------------------------------------------------

  s_vip     <- colSums(PLS2_Scores^2) * colSums(PLS2_Q^2)
  sum_s_vip <- sum(s_vip)

  # Loadings-based orthogonal VIP helper (unchanged, nort>0 only)
  .vip_ort_for_loadings_y <- function(ort_scores_mat, ort_loadings_b) {
    ncomp  <- ncol(ort_scores_mat)
    p_b    <- nrow(ort_loadings_b)
    ss_T   <- colSums(ort_scores_mat^2)
    if (any(ss_T < .Machine$double.eps)) return(rep(NA_real_, p_b))
    P_norm <- apply(ort_loadings_b, 2, function(pp) {
      nrm <- sqrt(sum(pp^2)); if (nrm > 1e-14) pp / nrm else pp
    })
    if (is.null(dim(P_norm))) P_norm <- matrix(P_norm, ncol = 1)
    sum_s <- sum(ss_T)
    if (sum_s < .Machine$double.eps) return(rep(NA_real_, p_b))
    as.vector(sqrt(p_b * (P_norm^2 %*% ss_T) / sum_s))
  }

  has_ort <- !is.null(Q_ort) && !is.null(orthoP)

  if (sum_s_vip < .Machine$double.eps) {
    warning("Sum of explained SS is zero; VIP scores cannot be computed.")
    VIP_global <- setNames(rep(NA_real_, length(varnames)), varnames)
    VIP_blocks <- setNames(
      lapply(seq_len(ntable), function(i) {
        v <- rep(NA_real_, variable_number[i])
        if (has_ort) data.frame(p=v, o=v, tot=v, row.names=MB@Variables[[i]])
        else         data.frame(p=v,     tot=v, row.names=MB@Variables[[i]])
      }),
      TableLabels
    )
  } else {
    # Global VIP: L2-normalise W columns across all variables
    W_norm_all <- apply(PLS2_W, 2, function(w) {
      nrm <- sqrt(sum(w^2)); if (nrm > 1e-14) w / nrm else w
    })
    p_total    <- nrow(PLS2_W)
    VIP_global <- setNames(
      as.vector(sqrt(p_total * (W_norm_all^2 %*% s_vip) / sum_s_vip)),
      varnames
    )

    # Per-block VIP
    VIP_blocks <- setNames(vector("list", ntable), TableLabels)
    k_vip     <- 1L
    k_ort_vip <- 1L
    for (i in seq_len(ntable)) {
      p_b <- variable_number[i]
      idx <- k_vip:(k_vip + p_b - 1L)

      # Predictive VIP (Wold): L2-normalise W using block rows only
      W_b_norm <- apply(PLS2_W[idx, , drop = FALSE], 2, function(w) {
        nrm <- sqrt(sum(w^2)); if (nrm > 1e-14) w / nrm else w
      })
      vip_p <- as.vector(sqrt(p_b * (W_b_norm^2 %*% s_vip) / sum_s_vip))

      if (has_ort) {
        ort_P_b <- orthoP[k_ort_vip:(k_ort_vip + p_b - 1L), , drop = FALSE]
        vip_o   <- .vip_ort_for_loadings_y(Q_ort, ort_P_b)
        vip_t   <- sqrt((vip_p^2 + vip_o^2) / 2)
        VIP_blocks[[TableLabels[i]]] <- data.frame(
          p = vip_p, o = vip_o, tot = vip_t, row.names = MB@Variables[[i]]
        )
        k_ort_vip <- k_ort_vip + p_b
      } else {
        VIP_blocks[[TableLabels[i]]] <- data.frame(
          p = vip_p, tot = vip_p, row.names = MB@Variables[[i]]
        )
      }
      k_vip <- k_vip + p_b
    }

    # Global VIP overwritten with "tot" to be consistent with per-block
    VIP_global <- setNames(
      as.vector(unlist(lapply(VIP_blocks, function(b) b$tot))),
      varnames
    )
  }
  ## -------------------------------------------------------------------------

  rm(saliences)

  # Define block membership of each variable in P.loadings
  belong_block <- rep(0, nrow(res_calib$P))
  k <- 1
  for (i in 1:ntable) {
    belong_block[k:(k + variable_number[i] - 1)] <- TableLabels[i]
    k <- k + variable_number[i]
  }

  ## Calculate Y predictions
  Ypred <- Xnorm_mat %*% ComDimPLS2_B

  for (i in seq_len(ncol(y))) {
    Ypred[, i] <- unlist(Ypred[, i]) + ComDimPLS2_B0[i]
  }
  rownames(Ypred) <- MB@Samples
  if (!is.null(colnames(y_raw))) colnames(Ypred) <- colnames(y_raw)

  ## Predict classes / compute Q2
  if (type == "discriminant") {
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
      E0 <- Ypred[neg, i] - y[neg, i] # Residuals of Class0 samples
      E1 <- Ypred[pos, i] - y[pos, i] # Residuals of Class1 samples

      E0count <- which(E0 > 0) # Class0 predictions above 0
      E1count <- which(E1 < 0) # Class1 predictions below 1

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
  } else if (type == "regression") {
    PRESS_r <- colSums((Ypred - y)^2)
    TSS_r <- colSums(sweep(y, 2, colMeans(y))^2)
    Q2 <- setNames(1 - PRESS_r / TSS_r, colnames(y))
  }

  ## -------------------------------------------------------------------------
  ## k-fold cross-validation: Q2 and DQ2
  ## Skipped when nort > 0 (CV with user-defined orthogonal FUNs is not
  ## generically supported) or when cv.k < 2.
  ## -------------------------------------------------------------------------
  if (cv.k >= 2 && nort == 0) {
    # Sequential fold assignment
    fold_ids <- integer(nrowMB)
    for (.f in seq_len(cv.k)) fold_ids[seq.int(.f, nrowMB, by = cv.k)] <- .f

    Ypred_cv <- matrix(NA_real_, nrow = nrowMB, ncol = ncol(y))
    colnames(Ypred_cv) <- colnames(y)

    for (fold in seq_len(cv.k)) {
      test_idx  <- seq.int(fold, nrowMB, by = cv.k)
      train_idx <- setdiff(seq_len(nrowMB), test_idx)

      # Subset MB to training samples
      MB_cv <- MB
      MB_cv@Data <- lapply(MB@Data, function(x) x[train_idx, , drop = FALSE])
      MB_cv@Samples <- MB@Samples[train_idx]

      # Fit inner model on training fold (cv.k = 0 prevents nested CV)
      invisible(utils::capture.output(
        res_cv_fold <- ComDim_y(
          MB = MB_cv,
          y = y_raw[train_idx, , drop = FALSE],
          ndim = ndim,
          FUN = FUN,
          nort = 0L,
          type = type,
          decisionRule = decisionRule,
          normalise = normalise,
          scale.y = scale.y,
          threshold = threshold,
          loquace = FALSE,
          method = method,
          cv.k = 0,
          ...
        )
      ))

      # Build raw test X (concatenated blocks, same order as Xnorm_mat)
      X_test <- matrix(0, nrow = length(test_idx), ncol = sum(variable_number))
      k_idx <- 1
      for (i in seq_len(ntable)) {
        vi <- variable_number[i]
        X_test[, k_idx:(k_idx + vi - 1)] <- MB@Data[[i]][test_idx, , drop = FALSE]
        k_idx <- k_idx + vi
      }

      # Apply training normalisation to test set (only when normalise = TRUE)
      if (normalise) {
        k_idx <- 1
        for (i in seq_len(ntable)) {
          vi <- variable_number[i]
          mn <- res_cv_fold@Mean$MeanMB[[i]] # training col means
          nrm <- res_cv_fold@Norm$NormMB[i] # training Frobenius norm
          X_test[, k_idx:(k_idx + vi - 1)] <-
            (X_test[, k_idx:(k_idx + vi - 1)] -
              matrix(mn, nrow = length(test_idx), ncol = vi, byrow = TRUE)) / nrm
          k_idx <- k_idx + vi
        }
      }

      # Predict Y for the test fold using the training model
      B_cv <- res_cv_fold@PLS.model$B
      B0_cv <- res_cv_fold@PLS.model$B0
      Ypred_cv[test_idx, ] <-
        X_test %*% B_cv +
        matrix(B0_cv, nrow = length(test_idx), ncol = ncol(y), byrow = TRUE)
    }

    # Q2 per response column (or class)
    TSS_cv <- colSums(sweep(y, 2, colMeans(y))^2)
    PRESS_cv <- colSums((Ypred_cv - y)^2)
    Q2_cv <- setNames(1 - PRESS_cv / TSS_cv, colnames(y))

    # DQ2 per class then averaged across classes (discriminant only).
    # Per-class values are stored in cv$DQ2.perclass; @DQ2 holds the mean,
    # consistent with ComDim_OPLS / ConsensusOPLS convention.
    if (type == "discriminant") {
      DQ2_per_class <- setNames(rep(NA_real_, ncol(y)), colnames(y))
      for (i in seq_len(ncol(y))) {
        pos <- which(classVect == colnames(y)[i])
        neg <- which(classVect != colnames(y)[i])
        E0 <- Ypred_cv[neg, i] - y[neg, i]
        E1 <- Ypred_cv[pos, i] - y[pos, i]
        PRESSD <- sum(E0[E0 > 0]^2) + sum(E1[E1 < 0]^2)
        DQ2_per_class[i] <- 1 - PRESSD / TSS_cv[i]
      }
      DQ2_cv <- mean(DQ2_per_class)   # single averaged value (matches ComDim_OPLS)
    } else {
      DQ2_per_class <- vector()
      DQ2_cv <- vector()
    }

    # Override training R2 with CV-based Q2/DQ2
    Q2 <- Q2_cv
    if (type == "discriminant") DQ2 <- DQ2_cv

    cv_result <- list(
      k            = cv.k,
      fold         = fold_ids,
      Ypred        = Ypred_cv,
      Q2           = Q2_cv,
      DQ2          = if (type == "discriminant") DQ2_cv        else NULL,
      DQ2.perclass = if (type == "discriminant") DQ2_per_class else NULL
    )
  } else {
    if (cv.k != 0 && nort == 0) {
      warning("'cv.k' must be >= 2 or 0 (skip CV). No cross-validation performed.")
    }
    cv_result <- list()
  }
  ## -------------------------------------------------------------------------

  end_output <- Sys.time()
  running_time <- (end_output - start_ini_time)
  if (loquace) {
    message(sprintf("Analysis finished after : %s seconds", running_time))
  }
  res_calib$runtime <- running_time # Total time of analysis

  progress_bar <- utils::txtProgressBar(min = 0, max = 80, style = 3, char = "=")
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
    Orthogonal = if (nort > 0) {
      list(
        nort = nort,
        Q.scores = Q_ort,
        T.scores = T_Loc_ort,
        P.loadings.ort = orthoP,
        Saliences.ort = saliences_ort,
        R2X = explained.ort
      )
    } else {
      list()
    },
    VIP = VIP_global,
    VIP.block = VIP_blocks,
    R2X = res_calib$explained,
    R2Y = r2y,
    Q2 = Q2,
    DQ2 = if (type == "discriminant") {
      DQ2
    } else {
      vector()
    },
    Singular = res_calib$SingVal,
    Mean = list(
      MeanMB = res_calib$MeanMB,
      MeanY = y_center, # original Y column means
      ScaleY = y_sd
    ), # original Y column SDs (1 if scale.y=FALSE)
    Norm = list(NormMB = res_calib$NormMB),
    PLS.model = list(
      W = PLS2_W,
      B = ComDimPLS2_B, # in original Y units
      B0 = ComDimPLS2_B0, # in original Y units
      Y = y_raw
    ), # original (unscaled) Y
    cv = cv_result,
    Prediction = if (type == "discriminant") {
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
    runtime = as.numeric(res_calib$runtime)
  )

  return(res_calib)
}
