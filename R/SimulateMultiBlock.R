#' SimulateMultiBlock
#'
#' Generate a synthetic \code{MultiBlock} dataset built from a known number of
#' orthogonal latent sources plus Gaussian noise.  Useful for benchmarking and
#' testing ComDim functions.
#'
#' The dataset is constructed as follows:
#' \enumerate{
#'   \item \code{n_sources} score vectors (\eqn{n \times \text{n\_sources}}) are
#'     drawn from a standard normal distribution and orthonormalised by QR
#'     decomposition.
#'   \item Loading vectors (\eqn{\text{n\_sources} \times p}) are built so that
#'     each source loads primarily (SD = 1) on one equal-sized variable segment,
#'     with small cross-loadings (SD = 0.10) on the remaining variables.
#'   \item The true signal \eqn{X = TP} is computed.
#'   \item Gaussian noise is added such that
#'     \eqn{\text{noise\_var} = \text{signal\_var} \times
#'     \text{noise} / (1 - \text{noise})}.
#'   \item The \eqn{p} variables are split into \code{n_blocks} equal-width
#'     blocks, each assembled as a named element of the returned
#'     \code{MultiBlock}.
#' }
#'
#' @param n Number of samples. Default: \code{500}.
#' @param p Total number of variables (split evenly across blocks). Must be
#'   divisible by \code{n_blocks} and by \code{n_sources}. Default: \code{2000}.
#' @param n_sources Number of orthogonal latent sources. Default: \code{4}.
#' @param noise Fraction of total variance attributed to noise, in (0, 1).
#'   Default: \code{0.05} (5 \% noise).
#' @param n_blocks Number of blocks to split the variables into. Default:
#'   \code{2}.
#' @return A \code{\link{MultiBlock}} object with \code{n_blocks} blocks, each
#'   of size \eqn{n \times (p / \text{n\_blocks})}, named \code{"Block1"},
#'   \code{"Block2"}, etc.
#' @seealso \code{\link{ComDim_PCA}}, \code{\link{MultiBlock}}
#' @examples
#' mb <- SimulateMultiBlock(n = 100, p = 200, n_sources = 4,
#'                          noise = 0.05, n_blocks = 2)
#' mb <- NormalizeMultiBlock(mb, method = 'norm')
#' res <- ComDim_PCA(mb, ndim = 4)
#' @export
SimulateMultiBlock <- function(n         = 500L,
                               p         = 2000L,
                               n_sources = 4L,
                               noise     = 0.05,
                               n_blocks  = 2L) {
  if (p %% n_blocks != 0L)
    stop("'p' must be divisible by 'n_blocks'.")
  if (p %% n_sources != 0L)
    stop("'p' must be divisible by 'n_sources'.")
  if (noise <= 0 || noise >= 1)
    stop("'noise' must be strictly between 0 and 1.")

  # Orthonormal score matrix (n x n_sources)
  T_orth <- qr.Q(qr(matrix(rnorm(n * n_sources), nrow = n, ncol = n_sources)))

  # Loading matrix (n_sources x p): each source loads primarily on one segment
  segment <- p / n_sources
  P <- matrix(0, nrow = n_sources, ncol = p)
  for (k in seq_len(n_sources)) {
    idx   <- ((k - 1L) * segment + 1L):(k * segment)
    other <- setdiff(seq_len(p), idx)
    P[k, idx]   <- rnorm(segment, sd = 1)
    P[k, other] <- rnorm(length(other), sd = 0.10)
  }

  X_signal   <- T_orth %*% P
  signal_var <- var(as.vector(X_signal))
  noise_var  <- signal_var * noise / (1 - noise)
  X          <- X_signal + matrix(rnorm(n * p, sd = sqrt(noise_var)),
                                  nrow = n, ncol = p)

  rownames(X) <- paste0("S", seq_len(n))
  colnames(X) <- paste0("V", seq_len(p))

  # Split into blocks
  block_width <- p / n_blocks
  blocks <- vector("list", n_blocks)
  names(blocks) <- paste0("Block", seq_len(n_blocks))
  for (b in seq_len(n_blocks)) {
    cols <- ((b - 1L) * block_width + 1L):(b * block_width)
    blocks[[b]] <- X[, cols, drop = FALSE]
  }

  MultiBlock(Data = blocks)
}
