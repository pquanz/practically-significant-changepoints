#' Change-point test for relevant hypotheses (wrt. Δ)
#'
#' Runs the paper's test by estimating the change point, the trimming
#' parameter `m` and (if `adaptive = TRUE`) also the support set S.
#'
#' @name changepoint_test
#' @param sample Numeric matrix or numeric data frame (rows = n, cols = p).
#' @param delta Numeric scalar (paper: \eqn{\Delta}). Size size of relevant
#'   change that will be tested for. If set to `NULL`, it only returns
#'   `max_delta_reject` and no test decision. Default `NULL`
#' @param m Trimming parameter. Either a nonnegative scalar applied to both
#'   sides or a length-2 numeric `c(m1 = ..., m2 = ...)`. Default `NULL`
#' @param adaptive Logical; if `TRUE` and `support_set = NULL`, estimates
#'   the support set via `estimate_support_set`.
#'   If `FALSE`, uses `1:p` as support_set.
#' @param support_set Optional subset of column indices (paper: set \eqn{S}).
#'   In this case, `adaptive` must be set to `TRUE`.
#' @param alpha Test level in (0, 1). Default `0.05`.
#' @param method_m Character, either `"max"` or `"individual"`. Only relevant
#'   if `m = NULL`. Options are:
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
#' @return A data.frame object with entries (not all are nec. given)
#'   - `tn`: value of either \eqn{T_n} or \eqn{T_{n, \hat S_n}}
#'           (depending on the value of `adaptive`)
#'   - `vn`: value of either \eqn{V_n} or \eqn{V_{n, \hat S_n}}
#'           (depending on the value of `adaptive`)
#'   - `khat`: named scalar `critical_value`
#'   - `reject`: Logical; decision of the test (only given if `delta`
#'     is not equal to `NULL`)
#'   - `delta_max`: the maximum value such that any relevant
#'     hypothesis with `delta <= delta_max` will reject the
#'     null hypotheses.
#'
#' @template parallelization
#'
#' @seealso [tn_statistic()], [quantile_g()]
cp_test <- function(
  sample,
  delta = NULL,
  m = NULL,
  adaptive = FALSE,
  support_set = NULL,
  alpha = 0.05,
  method_m = "individual",
  n_points = 20,
  parallelize = FALSE
) {
  khat <- estimate_changepoint(sample)

  # estimate trimming param m, if it is NULL
  estimate_m <- is.null(m)
  if (estimate_m) {
    m <- estimate_m(
      sample, khat, method_m = method_m, parallelize = parallelize
    )
  }

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
    tn        = round(tn, digits = 3),
    vn        = round(vn, digits = 3),
    khat      = khat,
    delta_max = round(max_delta_rej, digits = 3),
    n_points  = n_points
  )

  # add m estimators to data frame
  if (is.numeric(m) && all(c("m1", "m2") %in% names(m))) {
    df$m1 <- m[["m1"]]
    df$m2 <- m[["m2"]]
  } else {
    df$m <- m
  }

  if (!is.null(delta)) df$reject <- (max_delta_rej >= delta)
  if (adaptive) df$shat <- length(support_set)

  df
}



#' @rdname changepoint_test
#' @export
cp_test_sparsity_adj <- function(
  sample,
  delta = NULL,
  m = NULL,
  alpha = 0.05,
  method_m = "individual",
  n_points = 20,
  parallelize = FALSE
) {

  cp_test(
    sample, delta, m, TRUE,
    alpha = alpha,
    method_m = method_m,
    n_points = n_points,
    parallelize = parallelize
  )
}


#' @rdname changepoint_test
#' @export
cp_test_normalized <- function(
  sample,
  delta = NULL,
  m = NULL,
  method_m = "individual",
  alpha = 0.05,
  n_points = 20,
  parallelize = FALSE
) {

  cp_test(
    sample, delta, m, FALSE,
    method_m = method_m,
    alpha = alpha,
    n_points = n_points,
    parallelize = parallelize
  )
}