#' PredictMultiBlock
#'
#' Projects a new MultiBlock dataset into an existing ComDim model.
#' Works with models produced by ComDim_PCA, ComDim_PLS,
#' ComDim_OPLS, ComDim_y, and ComDim_Exploratory. The projection
#' type (PCA-like, PLS-like, OPLS-like) is determined from the model
#' structure rather than the method string, so custom method labels from
#' ComDim_y are handled automatically.
#'
#' @param MB A MultiBlock object containing the new samples to project.
#' @param y Response vector or dummy matrix (optional). When supplied for a
#'   supervised model, Q2, DQ2, and classification statistics are computed for
#'   the new samples.
#' @param model A ComDim object (the calibration model).
#' @param normalise If TRUE, each block is mean-centred using the training
#'   column means and divided by the training Frobenius norm. Must match the
#'   \code{normalise} setting used during calibration. Default FALSE.
#' @param loquace If TRUE, print a message for each set of model elements
#'   that were projected. Default TRUE.
#' @return The \code{model} ComDim object with updated slots:
#'   \itemize{
#'     \item \code{Q.scores} — projected global scores (new samples x ndim).
#'     \item \code{T.scores} — projected local scores (per block).
#'     \item \code{Orthogonal$Q.scores} — projected global ort scores (if
#'       model has orthogonal components).
#'     \item \code{Orthogonal$T.scores} — projected local ort scores.
#'     \item \code{Prediction$Y.pred} — predicted Y (supervised models).
#'     \item \code{Q2}, \code{DQ2}, classification slots — when \code{y}
#'       is supplied for a supervised model.
#'   }
#' @export
PredictMultiBlock <- function(MB = MB, y, model = model,
                              normalise = FALSE, loquace = TRUE) {
  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!inherits(MB, "MultiBlock")) {
    stop("'MB' is not a MultiBlock.")
  }

  if (!inherits(model, "ComDim")) {
    stop("'model' is not class ComDim.")
  }

  if (!(all(names(model@T.scores) %in% blockNames(MB)))) {
    stop("Not all the model blocks are included in MB.")
  }

  # ---------------------------------------------------------------------------
  # Block alignment: remove extra blocks; reorder/subset variables to match model
  # ---------------------------------------------------------------------------
  if (any(!(blockNames(MB) %in% names(model@T.scores)))) {
    dif_names <- setdiff(blockNames(MB), names(model@T.scores))
    for (i in rev(dif_names)) {
      MB@Data[[i]] <- NULL
      MB@Variables[[i]] <- NULL
    }
  }

  for (i in names(model@T.scores)) {
    model_vars_i <- rownames(model@P.loadings)[model@variable.block == i]
    if (all(model_vars_i %in% MB@Variables[[i]])) {
      idx <- match(model_vars_i, MB@Variables[[i]])
      MB@Data[[i]] <- MB@Data[[i]][, idx, drop = FALSE]
      MB@Variables[[i]] <- MB@Variables[[i]][idx]
    } else {
      stop(sprintf("Some variables required by the model are missing in block '%s'.", i))
    }
  }

  # ---------------------------------------------------------------------------
  # Determine model type from structure (works for any Method string)
  # ---------------------------------------------------------------------------
  is_supervised <- length(model@PLS.model) > 0 && !is.null(model@PLS.model$B)
  has_ort <- length(model@Orthogonal) > 0 &&
    !is.null(model@Orthogonal$P.loadings.ort)
  is_discriminant <- is_supervised && !is.null(model@Prediction$decisionRule)

  # ---------------------------------------------------------------------------
  # Process y (when supplied for a supervised model)
  # ---------------------------------------------------------------------------
  if (hasArg(y) && is_supervised) {
    if (is.vector(y) || is.matrix(y) || is.data.frame(y)) {
      n_new <- length(sampleNames(MB))
      if ((is.vector(y) && length(y) != n_new) ||
        (!is.vector(y) && nrow(y) != n_new)) {
        stop("'y' does not have the same number of rows as MB.")
      }
    }

    if (is_discriminant) {
      tmp <- unique(as.vector(y))
      if (all(tmp %in% c(0, 1))) { # already a dummy matrix
        if (!is.matrix(y)) y <- as.matrix(y)
        classVect <- rep(NA_character_, nrow(y))
        if (ncol(y) == 1) {
          classVect <- as.character(as.vector(y))
        } else {
          if (is.null(colnames(y))) {
            stop("Column names for 'y' must be provided when it is a dummy matrix.")
          }
          if (any(rowSums(y) != 1)) {
            if (any(rowSums(y) == 0)) stop("At least one sample was not assigned to a class.")
            if (any(colSums(y) > 1)) stop("At least one sample was assigned to more than one class.")
          }
          for (i in seq_len(ncol(y))) {
            classVect[which(y[, i] == 1)] <- colnames(y)[i]
          }
        }
      } else { # class vector → convert to dummy
        if ((is.matrix(y) || is.data.frame(y)) && ncol(y) > 1) {
          stop("For a discriminant model, 'y' must be a dummy matrix or a class vector.")
        }
        classVect <- as.character(as.vector(y))
        classy <- sort(unique(classVect))
        if (any(!(classy %in% colnames(model@PLS.model$Y)))) {
          stop("Not all classes in 'y' exist in the calibration model.")
        }
        y_dm <- matrix(0, nrow = length(classVect), ncol = length(classy))
        for (i in seq_along(classy)) {
          y_dm[classVect == classy[i], i] <- 1
        }
        colnames(y_dm) <- classy
        y <- y_dm
        rm(y_dm, classy)
      }
    } else { # regression
      if ((is.matrix(y) || is.data.frame(y)) && ncol(y) > 1) {
        stop("For a regression model, 'y' must be a vector or single-column matrix.")
      }
      if (is.vector(y)) y <- as.matrix(y)
    }
  }

  # ---------------------------------------------------------------------------
  # Build normalised Xnorm list and Xnorm_mat
  # ---------------------------------------------------------------------------
  names_table <- names(model@T.scores)
  nrowMB <- length(sampleNames(MB))
  p_total <- nrow(model@P.loadings)

  Xnorm <- list()
  Xnorm_mat <- matrix(0, nrow = nrowMB, ncol = p_total)

  for (i in names_table) {
    blk_idx <- which(model@variable.block == i) # column range in the global matrix

    if (normalise) {
      # Apply training mean-centring and training Frobenius norm
      trn_mean <- model@Mean$MeanMB[[i]]
      # NormMB can be a named list (normalise=FALSE training) or unnamed numeric
      # vector (normalise=TRUE training). Use block position as fallback index.
      blk_pos <- match(i, names_table)
      trn_norm <- tryCatch(
        as.numeric(model@Norm$NormMB[[i]]), # named access
        error = function(e) as.numeric(model@Norm$NormMB[blk_pos]) # positional
      )
      if (length(trn_norm) > 1) trn_norm <- trn_norm[1] # should be scalar
      X_mean <- MB@Data[[i]] -
        matrix(trn_mean, nrow = nrowMB, ncol = length(trn_mean), byrow = TRUE)
      Xnorm[[i]] <- X_mean / trn_norm
    } else {
      Xnorm[[i]] <- MB@Data[[i]]
    }
    Xnorm_mat[, blk_idx] <- Xnorm[[i]]
  }

  # ---------------------------------------------------------------------------
  # STEP 1 — Orthogonal deflation (only when model has ort components)
  # ---------------------------------------------------------------------------
  if (has_ort) {
    nort <- model@Orthogonal$nort
    ort_labs <- colnames(model@Orthogonal$P.loadings.ort)
    if (is.null(ort_labs)) ort_labs <- paste0("ort", seq_len(nort))

    Q_ort <- matrix(NA_real_,
      nrow = nrowMB, ncol = nort,
      dimnames = list(MB@Samples, ort_labs)
    )
    T_Loc_ort <- setNames(
      lapply(names_table, function(i) {
        matrix(NA_real_,
          nrow = nrowMB, ncol = nort,
          dimnames = list(MB@Samples, ort_labs)
        )
      }),
      names_table
    )

    for (jj in seq_len(nort)) {
      p_ort_global <- model@Orthogonal$P.loadings.ort[, jj] # p × 1
      norm2_ort <- sum(p_ort_global^2)

      # Global ort score: project full (current) Xnorm_mat onto ort loading
      Q_ort[, jj] <- as.vector(Xnorm_mat %*% p_ort_global) / norm2_ort

      # Local ort scores + block-level deflation
      for (i in names_table) {
        blk_idx <- which(model@variable.block == i)
        p_ort_loc <- p_ort_global[blk_idx] # p_i × 1
        norm2_loc <- sum(p_ort_loc^2)
        T_Loc_ort[[i]][, jj] <- as.vector(Xnorm[[i]] %*% p_ort_loc) / norm2_loc
        Xnorm[[i]] <- Xnorm[[i]] - Q_ort[, jj] %*% t(p_ort_loc)
      }

      # Rebuild Xnorm_mat from ort-deflated blocks
      Xnorm_mat <- do.call(
        cbind,
        lapply(names_table, function(i) Xnorm[[i]])
      )
    }

    model@Orthogonal$Q.scores <- Q_ort
    model@Orthogonal$T.scores <- T_Loc_ort
    if (loquace) {
      cat("Orthogonal Q.scores and T.scores were projected from the ComDim model.\n")
    }
  }

  # Capture Xnorm_mat here for Y prediction:
  #   - PCA/PLS: original (not ort-deflated) normalised X
  #   - OPLS: ort-deflated X  (B was fitted on this subspace)
  Xnorm_mat_for_ypred <- Xnorm_mat

  # ---------------------------------------------------------------------------
  # STEP 2 — Predictive score projection (sequential deflation)
  # ---------------------------------------------------------------------------
  ndim <- model@ndim
  dim_labs <- colnames(model@Q.scores)
  if (is.null(dim_labs)) dim_labs <- paste0("CC", seq_len(ndim))

  Q <- matrix(NA_real_,
    nrow = nrowMB, ncol = ndim,
    dimnames = list(MB@Samples, dim_labs)
  )
  T_Loc <- setNames(
    lapply(names_table, function(i) {
      matrix(NA_real_,
        nrow = nrowMB, ncol = ndim,
        dimnames = list(MB@Samples, dim_labs)
      )
    }),
    names_table
  )

  for (j in seq_len(ndim)) {
    p_global <- model@P.loadings[, j] # p × 1
    norm2_g <- sum(p_global^2)

    # Global score: project current Xnorm_mat onto the j-th global loading
    Q[, j] <- as.vector(Xnorm_mat %*% p_global) / norm2_g

    # Local scores + block-level deflation
    for (i in names_table) {
      blk_idx <- which(model@variable.block == i)
      p_local <- p_global[blk_idx] # p_i × 1
      norm2_l <- sum(p_local^2)
      T_Loc[[i]][, j] <- as.vector(Xnorm[[i]] %*% p_local) / norm2_l
      Xnorm[[i]] <- Xnorm[[i]] - Q[, j] %*% t(p_local)
    }

    # Deflate Xnorm_mat so the next component operates on the residual
    Xnorm_mat <- Xnorm_mat - Q[, j] %*% t(p_global)
  }

  model@Q.scores <- Q
  model@T.scores <- T_Loc
  if (loquace) {
    cat("Predictive Q.scores and T.scores were projected from the ComDim model.\n")
  }

  # ---------------------------------------------------------------------------
  # STEP 3 — Y prediction (supervised models only)
  # ---------------------------------------------------------------------------
  if (is_supervised) {
    Ypred <- Xnorm_mat_for_ypred %*% model@PLS.model$B
    for (i in seq_len(ncol(Ypred))) {
      Ypred[, i] <- Ypred[, i] + model@PLS.model$B0[i]
    }
    model@Prediction$Y.pred <- Ypred
    if (loquace) cat("Y.pred was predicted from the ComDim model.\n")

    # -----------------------------------------------------------------------
    # STEP 4 — Classification stats (discriminant + y supplied)
    # -----------------------------------------------------------------------
    if (hasArg(y) && is_discriminant) {
      predClass <- rep(NA_character_, nrow(y))
      dr <- model@Prediction$decisionRule

      if (dr == "max") {
        for (i in seq_len(nrow(Ypred))) {
          xx <- which(Ypred[i, ] == max(Ypred[i, ], na.rm = TRUE))
          if (length(xx) == 1) {
            predClass[i] <- colnames(y)[xx]
          } else {
            predClass[i] <- NA_character_
            warning(sprintf("Class for sample %d could not be predicted (tie).", i))
          }
        }
      } else if (dr == "fixed") {
        for (i in seq_len(nrow(Ypred))) {
          xx <- which(Ypred[i, ] > 1 / ncol(y))
          if (length(xx) == 1) {
            predClass[i] <- colnames(y)[xx]
          } else {
            predClass[i] <- NA_character_
            warning(sprintf("Class for sample %d could not be predicted (tie/no class).", i))
          }
        }
      }

      classy <- colnames(y)
      Q2 <- setNames(rep(NA_real_, length(classy)), classy)
      DQ2 <- setNames(rep(NA_real_, length(classy)), classy)
      specvec <- setNames(rep(NA_real_, length(classy)), classy)
      sensvec <- setNames(rep(NA_real_, length(classy)), classy)
      confusionMatrix <- list()

      for (i in seq_along(classy)) {
        pos <- which(classVect == classy[i])
        neg <- which(classVect != classy[i])
        tp <- sum(classVect == classy[i] & predClass == classy[i], na.rm = TRUE)
        tn <- sum(classVect != classy[i] & predClass != classy[i], na.rm = TRUE)
        fp <- sum(classVect != classy[i] & predClass == classy[i], na.rm = TRUE)
        fn <- sum(classVect == classy[i] & predClass != classy[i], na.rm = TRUE)

        sensvec[i] <- if (length(pos) > 0) tp / length(pos) else NA_real_
        specvec[i] <- if (length(neg) > 0) tn / length(neg) else NA_real_

        cm <- matrix(c(tp, fp, fn, tn), ncol = 2, nrow = 2)
        rownames(cm) <- c("predClass1", "predClass0")
        colnames(cm) <- c("trueClass1", "trueClass0")
        confusionMatrix[[classy[i]]] <- cm

        E0 <- Ypred[neg, i] - y[neg, i]
        E1 <- Ypred[pos, i] - y[pos, i]

        SSE0_all <- sum(E0^2)
        SSE1_all <- sum(E1^2)
        SSE0 <- sum(E0[E0 > 0]^2)
        SSE1 <- sum(E1[E1 < 0]^2)

        PRESSD <- SSE0 + SSE1
        PRESS <- SSE0_all + SSE1_all
        Ym <- y[, i] - mean(y[, i])
        TSS <- sum(Ym^2)

        DQ2[i] <- 1 - PRESSD / TSS
        Q2[i] <- 1 - PRESS / TSS
      }

      model@Q2 <- Q2
      model@DQ2 <- DQ2
      model@Prediction$trueClass <- classVect
      model@Prediction$predClass <- data.frame(predClass = predClass, row.names = MB@Samples)
      model@Prediction$Sensitivity <- sensvec
      model@Prediction$Specificity <- specvec
      model@Prediction$confusionMatrix <- confusionMatrix
      if (loquace) {
        cat(
          "Q2, DQ2, predClass, Sensitivity, Specificity, and confusionMatrix",
          "were computed for the new samples.\n"
        )
      }
    } else if (hasArg(y) && !is_discriminant) {
      # -----------------------------------------------------------------------
      # STEP 5 — Q2 for regression (y supplied, not discriminant)
      # -----------------------------------------------------------------------
      PRESS_r <- colSums((Ypred - y)^2)
      TSS_r <- colSums(sweep(y, 2, colMeans(y))^2)
      model@Q2 <- setNames(1 - PRESS_r / TSS_r, colnames(as.matrix(y)))
      if (loquace) cat("Q2 was computed for the new samples.\n")
    }
  }

  return(model)
}
