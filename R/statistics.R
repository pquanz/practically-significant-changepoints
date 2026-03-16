#' @keywords internal
#' @noRd
nm_factor <- function(k, m) {
  if (m >= 0) {
    ifelse(k > m, (k - m) * (k - m - 1), 0)
  } else {
    k^2
  }
}


#' T_n (k, m; λ) statistic (paper)
#'
#' Computes the statistic referred to in the manuscript as T_n (k, m; λ).
#' The argument `k` denotes the (candidate) change-point location in the data.
#' Although the statistic can be evaluated for any valid `k` with respect to
#' `n`, in applications `k` is typically interpreted as the change point.
#'
#' @param sample Numeric matrix or data frame with observations in rows (n)
#'   and variables in columns (p).
#' @param k Integer index (1-based) denoting the change-point location.
#'   Must lie within the valid range for the data size.
#' @param m Trimming parameter(s). Either a nonnegative scalar (applied to both
#'   sides) or a length-2 numeric giving `c(m1 = ..., m2 = ...)`. Defaults to
#'   `c(m1 = 0, m2 = 0)`.
#' @param lambda Tuning parameter (paper: λ). Numeric in (0, 1]; controls the
#'   self-normalization subsample. Smaller values use a more restrictive subset;
#'   `1` uses the full eligible subset. Default `1`.
#' @param support_set Optional subset of coordinates (columns) to use, that is,
#'   a set of integer indices. If `NULL`, all columns are used.
#'
#' @return Numeric scalar: the value of T_n (k, m; λ).
#'
#' @details
#' An empty set `support_set` results in a value of 0 for the statistic.
#'
#' @examples
#' X <- matrix(rnorm(200 * 10), 200, 10)
#'
#' # default trimming
#' t_lambda_statistic(X, k = 100)
#'
#' # Scalar trimming for both sides
#' t_lambda_statistic(X, k = 120, m = 2)
#'
#' # Scalar trimming with unequal trimming parameters
#' t_lambda_statistic(X, k = 120, m = c(m1 = 1, m2 = 3))
#' @export
t_lambda_statistic_old <- function(
  sample, k, m = c(m1 = 0, m2 = 0), lambda = 1, support_set = NULL
) {

  if (all(c("m1", "m2") %in% names(m)) && length(m) == 2) {
    # m is of the form m = c(m1 = 3, m2 = 10)
    m1 <- m[["m1"]]
    m2 <- m[["m2"]]
  } else if (is.numeric(m)) {
    # m is just a number
    m1 <- m
    m2 <- m
  } else {
    stop("Please set either m to a number or a list as c(m1 = ..., m2 = ...)")
  }

  n <- NROW(sample)
  p <- NCOL(sample)

  if (k <= m1 + 1 || k + m2 >= n - 1) return(0)

  klam <- floor(k * lambda)
  nklam <- floor((n - k - m2) * lambda)

  if (identical(support_set, integer(0))) {  # check for empty set S
    return(0)
  } else if (!is.null(support_set)) {
    # check for nonempty given S, else (i.e. S = NULL, it defaults to S = 1:p)
    sample <- sample[, support_set, drop = FALSE]
    p <- length(support_set)
  }


  t1 <- 0
  for (i in 1:klam) {
    # list of indices j = 1, ... lambda k which satisfy |j - i| > m1
    indices_m_apart <- which(abs(i - 1:klam) > m1)

    # prepare matrix of those sample items in above list
    trimmed_sample_t <- t(sample[indices_m_apart, , drop = FALSE])
    t1 <- t1 + sum(sample[i, ] %*% trimmed_sample_t)
  }


  t2 <- 0
  for (j in (k + m2 + 1):(k + m2 + nklam)) {
    # indices i = k + m2 + 1, ... k + m2 + lambda * (n - k - m2), which satisfy |j - i| > m2
    indices_m_apart <- k + m2 + which(abs(j - (k + m2 + 1):(k + m2 + nklam)) > m2)
    trimmed_sample_t <- t(sample[indices_m_apart, , drop = FALSE])
    t2 <- t2 + sum(sample[j, ] %*% trimmed_sample_t)
  }


  t31 <- 0
  t32 <- 0
  if (klam > m1 && nklam > m2 && m1 >= 0 && m2 >= 0) {
    cs31 <- colSums(sample[(1:klam), , drop = FALSE])
    t31 <- (klam - 2 * m1 - 1) * cs31
    if (m1 > 0) {
      for (i in 1:m1) {
        t31 <- t31 + i * sample[i + klam - m1, ] - (i - m1 - 1) * sample[i, ]
      }
    }
    cs32 <- colSums(sample[((k + m2 + 1): (k + m2 + nklam)), , drop = FALSE])
    t32 <- (nklam - 2 * m2 - 1) * cs32
    if (m2 > 0) {
      for (j in 1:m2) {
        t32 <- t32 +
          j * sample[(j + k + nklam), ] -
          (j - m2 - 1) * sample[(j + k + m2), ]
      }
    }
  } else if (klam > m1 && nklam > m2) {
    # here m < 0
    t31 <- klam * colSums(sample[(1:klam), , drop = FALSE])
    t32 <- nklam * colSums(sample[((k + m2 + 1): (k + m2 + nklam)), , drop = FALSE])

  } # else the sum will be zero

  sum <- nm_factor(nklam, m2) * t1 + nm_factor(klam, m1) * t2 -
    2 * (t31 %*% t32)
  sum <- sum / (p * nm_factor(n - k - m2, m2) * nm_factor(k, m1))

  sum[1, 1]
}








# --- fast helpers (internal) ---
.gap_sum_block <- function(X, m) {
  u <- nrow(X)
  if (u <= 0) return(0)
  if (m < 0) {
    # |a-b| > m always true, includes diagonal -> sum_{a,b} Xa^T Xb = ||sum X||^2
    s <- colSums(X)
    return(sum(s * s))
  }
  if (u <= m + 1) return(0)

  s <- colSums(X)
  offdiag <- sum(s * s) - sum(rowSums(X * X))  # sum_{a!=b} Xa^T Xb

  lag_sum <- 0
  if (m > 0) {
    for (h in 1:m) {
      if (u - h >= 1) {
        lag_sum <- lag_sum + sum(rowSums(X[1:(u - h), , drop = FALSE] *
                                        X[(1 + h):u, , drop = FALSE]))
      }
    }
  }
  offdiag - 2 * lag_sum
}

.gap_weighted_sum <- function(X, m) {
  u <- nrow(X)
  if (u <= 0) return(rep(0, ncol(X)))
  if (m < 0) {
    # all pairs allowed incl diagonal -> each row pairs with u rows
    return(u * colSums(X))
  }
  if (u <= m + 1) return(rep(0, ncol(X)))

  idx <- 1:u
  w <- (u - 1) - pmin(m, idx - 1) - pmin(m, u - idx)  # w_u(t)
  colSums(X * w)
}

# --- DROP-IN replacement: same name/signature ---
#' @export
t_lambda_statistic <- function(
  sample, k, m = c(m1 = 0, m2 = 0), lambda = 1, support_set = NULL
) {
  if (all(c("m1", "m2") %in% names(m)) && length(m) == 2) {
    m1 <- m[["m1"]]
    m2 <- m[["m2"]]
  } else if (is.numeric(m)) {
    m1 <- m
    m2 <- m
  } else {
    stop("Please set either m to a number or a list as c(m1 = ..., m2 = ...)")
  }

  n <- NROW(sample)
  p <- NCOL(sample)

  # same early exit as current code
  if (k <= m1 + 1 || k + m2 >= n - 1) return(0)

  klam  <- floor(k * lambda)
  nklam <- floor((n - k - m2) * lambda)

  if (identical(support_set, integer(0))) {
    return(0)
  } else if (!is.null(support_set)) {
    sample <- sample[, support_set, drop = FALSE]
    p <- length(support_set)
  }

  # if the factors are zero, return 0 (avoid 0/0)
  den1 <- nm_factor(k, m1)
  den2 <- nm_factor(n - k - m2, m2)
  if (p == 0 || den1 == 0 || den2 == 0) return(0)
  if (nm_factor(klam, m1) == 0 || nm_factor(nklam, m2) == 0) return(0)

  # blocks for lambda
  L <- sample[1:klam, , drop = FALSE]
  R <- sample[(k + m2 + 1):(k + m2 + nklam), , drop = FALSE]

  t1 <- .gap_sum_block(L, m1)
  t2 <- .gap_sum_block(R, m2)
  t31 <- .gap_weighted_sum(L, m1)
  t32 <- .gap_weighted_sum(R, m2)

  num <- nm_factor(nklam, m2) * t1 + nm_factor(klam, m1) * t2 - 2 * sum(t31 * t32)
  out <- num / (p * den2 * den1)
  as.numeric(out)
}













#' Self-normalizing statistic (paper: \eqn{V_n} and \eqn{\hat v_\ell})
#'
#' Computes the self-normalization quantity used in the manuscript for the
#' statistic T_n. The argument `k` denotes the (candidate) change-point
#' location (1-based). The computation can be restricted to a subset of
#' coordinates (the paper's set S).
#'
#' @param sample Numeric matrix (rows = n, columns = p).
#' @param k Integer index (1-based) denoting the change-point location.
#' @param m Trimming parameter(s). Either a nonnegative scalar (applied to
#'   both sides) or a length-2 numeric giving `c(m1 = ..., m2 = ...)`.
#'   Defaults to `c(m1 = 0, m2 = 0)`.
#' @param support_set Optional subset of coordinates (columns) to use, that is,
#'   a set of integer indices. If `NULL`, all columns are used.
#' @param n_points Integer (paper: K). Number of points drawn from the
#'   measure ν on (0, 1) used in the self-normalization construction.
#'   Default is 20.
#' @param type Integer code selecting the normalization variant. Use `1`
#'   for \eqn{V_n} and `0` for \eqn{\hat v_\ell}. In the estimation of the
#'   set S, type 0 is used.
#'
#' @return Numeric scalar: the value of the self-normalization factor.
#' @details An empty set results in a value of 0 for the normalizer.
#'
#' @keywords internal
normalizer <- function(
  sample, k,
  m = c(m1 = 0, m2 = 0),
  support_set = NULL,
  n_points = 20,
  type = 1
) {

  if (all(c("m1", "m2") %in% names(m)) && length(m) == 2) {
    # m is of the form m = c(m1 = 3, m2 = 10)
    m1 <- m[["m1"]]
    m2 <- m[["m2"]]
  } else if (is.numeric(m)) {
    # m is just a number
    m1 <- m
    m2 <- m
  } else {
    stop("Please set either m to a number or a list as c(m1 = ..., m2 = ...)")
  }


  n <- nrow(sample)
  lambdas <- 1:(n_points - 1) / n_points

  t_at_1 <- tn_statistic(sample, k, m, support_set = support_set)


  v_integrand <- function(lam) {
    if (type == 0) {
      (t_lambda_statistic(sample, k, m, lam, support_set) - lam^4 * t_at_1)^2
    } else if (type == 1) {
      klam <- floor(lam * k)
      nklam <- floor(lam * (n - k - m2))
      factor <- nm_factor(klam, m1) / nm_factor(k, m1) *
        nm_factor(nklam, m2) / nm_factor(n - k - m2, m2)
      (t_lambda_statistic(sample, k, m, lam, support_set) - factor * t_at_1)^2
    } else {
      stop("Type not accepted. Please only use type 0 or 1.")
    }
  }

  tlambdas <- sapply(lambdas, v_integrand, simplify = TRUE)
  sqrt(mean(tlambdas))
}



#' Pivotal set estimator (paper: \eqn{\hat S_n})
#'
#' Computes the estimator for the support set S, which is a subset of 1:p.
#'
#' @param sample Numeric matrix (rows = n, columns = p).
#' @param k Integer index (1-based) denoting the change-point location.
#' @param m Trimming parameter(s). Either a nonnegative scalar (applied to
#'   both sides) or a length-2 numeric giving `c(m1 = ..., m2 = ...)`.
#'   Defaults to `c(m1 = 0, m2 = 0)`.
#' @param kappa Numeric exponent of the logarithm `log(p)` in the definition
#'   of the set estimator.
#' @param n_points Integer (paper: K). Number of points drawn from the
#'   measure ν on (0, 1) used in the self-normalization construction.
#'   Default is 20.
#' @param parallelize Logical; if `TRUE`, evaluates the internal loop in
#'   parallel.
#'
#' @return List of indices or `integer(0)` if the empty set was estimated.
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

  p <- ncol(sample)

  check_coordinate <- function(i) {
    delta_i <- pscp::tn_statistic(sample, k, m, c(i))
    v_i <- pscp::vl_statistic(sample, k, m, c(i), n_points)
    delta_i > v_i * log(p)^kappa
  }

  indices <- .map_apply(1:p, check_coordinate, parallelize = parallelize)
  which(indices)
}



# --- Aliases ---

#' @rdname normalizer
#' @export
vn_statistic <- function(
  sample, k,
  m = c(m1 = 0, m2 = 0),
  support_set = NULL,
  n_points = 20
) {
  normalizer(sample, k, m, support_set, n_points)
}

#' @rdname normalizer
#' @export
vl_statistic <- function(
  sample, k,
  m = c(m1 = 0, m2 = 0),
  support_set = NULL,
  n_points = 20
) {
  normalizer(sample, k, m, support_set, n_points, 0)
}

#' @rdname t_lambda_statistic
#' @export
tn_statistic <- function(sample, k, m = c(m1 = 0, m2 = 0), support_set = NULL) {
  t_lambda_statistic(sample, k, m, 1, support_set)
}