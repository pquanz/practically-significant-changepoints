# --- Vectorized helpers for all-coordinates-at-once computation ---
# These return vectors of length p (one value per column).


#' @keywords internal
#' @noRd
.gap_sum_block_all <- function(X, m) {
  u <- nrow(X)
  p <- ncol(X)
  if (u <= 0) return(rep(0, p))
  if (m < 0) {
    s <- colSums(X)
    return(s * s)
  }
  if (u <= m + 1) return(rep(0, p))

  s <- colSums(X)
  offdiag <- s * s - colSums(X * X)

  lag_sum <- rep(0, p)
  if (m > 0) {
    for (h in 1:m) {
      if (u - h >= 1) {
        lag_sum <- lag_sum + colSums(
          X[1:(u - h), , drop = FALSE] * X[(1 + h):u, , drop = FALSE]
        )
      }
    }
  }
  offdiag - 2 * lag_sum
}


#' @keywords internal
#' @noRd
.gap_weighted_sum_all <- function(X, m) {
  u <- nrow(X)
  p <- ncol(X)
  if (u <= 0) return(rep(0, p))
  if (m < 0) return(u * colSums(X))
  if (u <= m + 1) return(rep(0, p))

  idx <- 1:u
  w <- (u - 1) - pmin(m, idx - 1) - pmin(m, u - idx)
  as.numeric(crossprod(X, w))
}


#' Compute t_lambda_statistic for ALL columns simultaneously
#'
#' Returns a numeric vector of length p, where entry ℓ equals
#' `t_lambda_statistic(sample, k, m, lambda, support_set = c(ℓ))`.
#'
#' @param sample Numeric matrix (rows = n, columns = p).
#' @param k Integer change-point location (1-based).
#' @param m1,m2 Trimming parameters (scalars).
#' @param lambda Tuning parameter in (0, 1].
#' @return Numeric vector of length p.
#'
#' @keywords internal
#' @noRd
.t_lambda_all_coords <- function(sample, k, m1, m2, lambda = 1) {
  n <- nrow(sample)
  p <- ncol(sample)

  if (k <= m1 + 1 || k + m2 >= n - 1) return(rep(0, p))

  klam  <- floor(k * lambda)
  nklam <- floor((n - k - m2) * lambda)

  den1 <- nm_factor(k, m1)
  den2 <- nm_factor(n - k - m2, m2)
  if (den1 == 0 || den2 == 0) return(rep(0, p))

  fklam  <- nm_factor(klam, m1)
  fnklam <- nm_factor(nklam, m2)
  if (fklam == 0 || fnklam == 0) return(rep(0, p))

  L <- sample[1:klam, , drop = FALSE]
  R <- sample[(k + m2 + 1):(k + m2 + nklam), , drop = FALSE]

  t1  <- .gap_sum_block_all(L, m1)
  t2  <- .gap_sum_block_all(R, m2)
  t31 <- .gap_weighted_sum_all(L, m1)
  t32 <- .gap_weighted_sum_all(R, m2)

  # p = 1 per coordinate (support_set = c(ℓ)), so no division by p
  num <- fnklam * t1 + fklam * t2 - 2 * t31 * t32
  num / (den2 * den1)
}


#' Fast estimation of the support set (paper: \eqn{\hat S_n})
#'
#' Drop-in replacement for `estimate_support_set` that vectorizes
#' across all p coordinates instead of looping.
#' Speedup is roughly p× in R overhead.
#'
#' @param sample Numeric matrix (rows = n, columns = p).
#' @param k Integer index (1-based) denoting the change-point location.
#' @param m Trimming parameter(s). Either a nonnegative scalar or
#'   a length-2 numeric giving `c(m1 = ..., m2 = ...)`.
#'   Defaults to `c(m1 = 0, m2 = 0)`.
#' @param kappa Numeric exponent of `log(p)` in the set estimator.
#' @param n_points Integer (paper: K). Default is 20.
#'
#' @return Integer vector of indices, or `integer(0)`.
#'
#' @examples
#' X <- matrix(rnorm(200 * 10), 200, 10)
#' s_hat <- estimate_support_set(X, k = 100)
#' @export
estimate_support_set <- function(
  sample, k,
  m = c(m1 = 0, m2 = 0),
  kappa = 3 / 2,
  n_points = 20,
  parallelize = FALSE
) {
  # parse m
  if (all(c("m1", "m2") %in% names(m)) && length(m) == 2) {
    m1 <- m[["m1"]]
    m2 <- m[["m2"]]
  } else if (is.numeric(m)) {
    m1 <- m
    m2 <- m
  } else {
    stop("Please set either m to a number or a list as c(m1 = ..., m2 = ...)")
  }

  p <- ncol(sample)

  # delta_hat^2_ℓ = T_{n,{ℓ}}(k, m; 1) for all ℓ at once
  delta_sq <- .t_lambda_all_coords(sample, k, m1, m2, lambda = 1)

  # v_hat_ℓ = (mean over λ of (T_{n,{ℓ}}(k,m;λ) - λ^4 * delta_sq_ℓ)^2)^{1/2}
  lambdas <- 1:(n_points - 1) / n_points

  sum_sq <- rep(0, p)
  for (lam in lambdas) {
    t_lam <- .t_lambda_all_coords(sample, k, m1, m2, lam)
    sum_sq <- sum_sq + (t_lam - lam^4 * delta_sq)^2
  }
  v_hat <- sqrt(sum_sq / length(lambdas))

  which(delta_sq > v_hat * log(p)^kappa)
}