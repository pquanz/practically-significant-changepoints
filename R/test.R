#' Change-point test for relevant hypotheses (wrt. Δ)
#'
#' Runs the paper's test by estimating the change point, the trimming
#' parameter `m` and (if `adaptive = TRUE`) also the support set S.
#'
#' @name changepoint_test
#' @param sample Numeric matrix or numeric data frame (rows = n, cols = p).
#' @param delta Numeric scalar (paper: \eqn{\Delta}). Size size of relevant
#'   change that will be tested for.
#' @param m Trimming parameter. Either a nonnegative scalar applied to both
#'   sides or a length-2 numeric `c(m1 = ..., m2 = ...)`.
#' @param adaptive Logical; if `TRUE` and `support_set = NULL`, estimates
#'   the support set via `estimate_support_set`.
#'   If `FALSE`, uses `1:p` as support_set.
#' @param support_set Optional subset of column indices (paper: set \eqn{S}).
#'   In this case, `adaptive` must be set to `TRUE`.
#' @param alpha Test level in (0, 1). Default `0.05`.
#' @param method_m Character, either `"max"` or `"individual"`.
#'   - `"individual"`: returns `c(m1 = ..., m2 = ...)`, where m1 is the
#'     estimated trimming parameter before the change point and m2 the
#'     estimated trimming parameter after the change point.
#'   - `"max"`: select a single trimming index \eqn{\hat m} by taking the
#'     maximum of the estimators obtained as described in individual.
#'   Default is `"max"`.
#' @param n_points Integer (paper: \eqn{K}); number of grid points used by
#'   the self-normalizing statistic \eqn{V_n}. Default `20`.
#' @param parallelize Logical; if `TRUE`, parallelizes internal loops when a
#'   backend is registered; otherwise runs sequentially. Default `TRUE`.
#'
#' @return A data.frame object with entries
#'   - `tn`: value of either \eqn{T_n} or \eqn{T_{n, \hat S_n}}
#'           (depending on the value of `adaptive`)
#'   - `vn`: value of either \eqn{V_n} or \eqn{V_{n, \hat S_n}}
#'           (depending on the value of `adaptive`)
#'   - `khat`: named scalar `critical_value`
#'   - `reject`: Logical; decision of the test
#'   - `max_delta_reject`: the maximum value such that any relevant
#'     hypothesis with `delta <= max_delta_reject` will reject the
#'     null hypotheses.
#'
#' @template parallelization
#'
#' @seealso [tn_statistic()], [quantile_g()]
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200 * 10), 200, 10)
#' # Minimal call (adaptive, no trimming)
#' # out <- cp_test(X, delta = 2, alpha = 0.05)
#' # Fixed trimming example, nonadaptive
#' # out <- cp_test(X, delta = 1, m = 2, adaptive = FALSE)
#' @export
cp_test <- function(
  sample,
  delta,
  m = c(m1 = 0, m2 = 0),
  adaptive = TRUE,
  support_set = NULL,
  alpha = 0.05,
  n_points = 20,
  parallelize = FALSE
) {
  khat <- estimate_changepoint(sample)

  if (adaptive) {
    if (is.null(support_set)) {
      support_set <- estimate_support_set(
        sample, khat, m, parallelize = parallelize
      )
    }
    tn <- tn_statistic(sample, khat, m, support_set)
    vn <- vn_statistic(sample, khat, m, support_set, n_points)
  } else {
    tn <- tn_statistic(sample, khat, m)
    vn <- vn_statistic(sample, khat, m, n_points = n_points)
  }

  quantile <- quantile_g(alpha = alpha, n_points = n_points)
  max_delta_rej <- tn - quantile * vn
  df <- data.frame(
    tn = round(tn, digits = 3),
    vn = round(vn, digits = 3),
    khat = khat,
    reject = max_delta_rej >= delta,
    max_delta_reject = round(max_delta_rej, digits = 3)
  )

  if (adaptive) df$shat <- length(support_set)

  df
}



#' @rdname changepoint_test
#' @export
cp_test_sparsity_adj <- function(
  sample,
  delta,
  alpha = 0.05,
  method_m = "individual",
  n_points = 20,
  parallelize = FALSE
) {
  # everything will be estimated here
  khat <- estimate_changepoint(sample)
  m <- estimate_m(sample, khat, method = method_m, parallelize = parallelize)

  df <- cp_test(
    sample, delta, m,
    alpha = alpha,
    n_points = n_points,
    parallelize = parallelize
  )

  if (method_m == "individual") {
    df$m1 <- m[["m1"]]
    df$m2 <- m[["m2"]]
  } else {
    df$m <- m
  }
  df
}


#' @rdname changepoint_test
#' @export
cp_test_normalized <- function(
  sample,
  delta,
  method_m = "individual",
  alpha = 0.05,
  n_points = 20,
  parallelize = FALSE
) {
  khat <- estimate_changepoint(sample)
  m <- estimate_m(sample, khat, method = method_m, parallelize = parallelize)

  df <- cp_test(
    sample, delta, m,
    adaptive = FALSE,
    alpha = alpha,
    n_points = n_points,
    parallelize = parallelize
  )

  if (method_m == "individual") {
    df$m1 <- m[["m1"]]
    df$m2 <- m[["m2"]]
  } else {
    df$m <- m
  }
  df
}