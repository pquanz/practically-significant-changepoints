#' @describeIn estimate_m \eqn{F_1} function for the estimation of m (paper)
#'
#' Computes the quantity denoted \eqn{F_1} in the manuscript for a given
#' candidate change-point location `k`. The input is a numeric matrix.
#' Trimming is capped by a proportion of the available segment length
#' (see `max_proportion_m`).
#'
#' @return A data frame with two columns:
#'   - `m` (integer): the trimming value/index evaluated (sorted ascending).
#'   - `f1` (numeric): the corresponding f1 value at that `m`.
#'
#' @examples
#' X <- matrix(rnorm(200 * 10), 200, 10)
#' get_f1(X, k = 100)
#' get_f1(X, k = 120, max_proportion_m = 1 / 4, parallelize = TRUE)
#' @export
get_f1 <- function(sample, k, max_proportion_m = 1 / 3, parallelize = FALSE) {
  p <- ncol(sample)

  max_m1 <- floor(k * max_proportion_m)

  compute_for_m1 <- function(m) {
    t1 <- 0
    cmeans <- matrix(colMeans(sample[1:k, ]), nrow = k, ncol = p, byrow = TRUE)
    centered_sample_1 <- sample[1:k, ] - cmeans
    for (i in 1:k) {
      # get a list of indices j = 1, ... lambda k which satisfy |j - i| > m
      indices_m_apart <- which(abs(i - 1:k) > m)

      # prepare matrix of those sample items in above list
      trimmed_sample_t <- t(centered_sample_1[indices_m_apart, ])
      t1 <- t1 + sum(centered_sample_1[i, ] %*% trimmed_sample_t)
    }
    t1 / ((k - m) * (k - m - 1))
  }

  f1 <- .map_apply(seq_len(max_m1), compute_for_m1, parallelize = parallelize)
  f1 / compute_for_m1(0)
}


#' @describeIn estimate_m Difference Function Delta F1
#' @examples
#' get_delta_f1(X, k = 100)
#' @export
get_delta_f1 <- function(
  sample, k, max_proportion_m = 1 / 3, parallelize = FALSE
) {
  f1 <- get_f1(sample, k, max_proportion_m, parallelize)
  del_f1 <- abs(c(f1, 0) - c(0, f1))[1:(length(f1) - 1)]
  data.frame(m = seq_len(length(del_f1)), F = del_f1)
}



#' @describeIn estimate_m F_2 function for the estimation of m (paper)
#'
#' Computes the quantity denoted F_2 in the manuscript for a given
#' candidate change-point location `k`. The input is a numeric matrix.
#' Trimming is capped by a proportion of the available segment length
#' (see `max_proportion_m`).
#'
#' @return A data frame with two columns:
#'   - `m` (integer): the trimming value/index evaluated (sorted ascending).
#'   - `f2` (numeric): the corresponding f2 value at that `m`.
#'
#' @examples
#' X <- matrix(rnorm(200 * 10), 200, 10)
#' get_f2(X, k = 100)
#' get_f2(X, k = 120, max_proportion_m = 1 / 4, parallelize = TRUE)
#' @export
get_f2 <- function(sample, k, max_proportion_m = 1 / 3, parallelize = FALSE) {
  n <- nrow(sample)
  p <- ncol(sample)

  max_m2 <- floor((n - k) * max_proportion_m)

  compute_for_m2 <- function(m) {
    t2 <- 0
    cmeans <- matrix(
      colMeans(sample[(k + 1): n, ]), nrow = n - k, ncol = p, byrow = TRUE
    )
    centered_sample_2 <- sample[(k + 1):n, ] - cmeans
    for (j in 1:(n - k)) {
      # get a list of indices i = k + 1, ... k + lambda * (n - k)
      # which satisfy |j - i| > m
      indices_m_apart <- which(abs(j - 1:(n - k)) > m)
      # prepare matrix of those sample items in above list
      trimmed_sample_t <- t(centered_sample_2[indices_m_apart, ])
      t2 <- t2 + sum(centered_sample_2[j, ] %*% trimmed_sample_t)
    }
    t2 / ((n - k - m) * (n - k - m - 1))
  }

  f2 <- .map_apply(seq_len(max_m2), compute_for_m2, parallelize = parallelize)
  f2 / compute_for_m2(0)
}


#' @describeIn estimate_m Difference Function Delta F2
#'
#' @examples
#' get_delta_f2(X, k = 100)
#' @export
get_delta_f2 <- function(
  sample, k, max_proportion_m = 1 / 3, parallelize = FALSE
) {
  f2 <- get_f2(sample, k, max_proportion_m, parallelize)
  del_f2 <- abs(c(f2, 0) - c(0, f2))[1:(length(f2) - 1)]
  data.frame(m = seq_len(length(del_f2)), F = del_f2)
}


#' Estimation of trimming parameter \eqn{\hat m_n}
#'
#' Estimates the trimming parameter for a given change point
#' candidate `k` using the difference function Delta F1. The maximum
#' value can only be a proportion of the whole respective segmet
#' (see `max_proportion_m`).
#'
#' @name estimate_m
#' @param sample Numeric matrix or data frame (rows = n, columns = p).
#'   If a data frame is supplied, all columns must be numeric.
#' @param k Integer index (1-based) denoting the candidate change-point
#'   location.
#' @param threshold Numeric scalar used by threshold-based rules (see `method`).
#'   Default `0.01`.
#' @param max_proportion_m Numeric in (0, 1]; caps the maximum `m` considered
#'   as `floor(max_proportion_m * segment_length)`, where segment_length denotes
#'   either the sample before or after the change point. Default `1/3`.
#' @param method Character selecting the return value. Supported:
#'   - `"max"`: pick the maximum of m1 and m2 and return numeric m.
#'   - `"individual"`: keep both values m1 and m2 and return both as list.
#'   Default `"max"`.
#' @param parallelize Logical; if `TRUE`, evaluates the internal loop in
#'   parallel.
#'
#' @return Depending on method:
#'   - if `method` is set to `max`, it returns an integer scalar,
#'   - if `method` is set to `individual`, it returns a list
#'     c(m1 = ..., m2 = ...)
#'
#' @template parallelization
#'
#' @seealso [get_f1], [get_delta_f1()], [get_f2()], [get_delta_f2()]
#'
#' @examples
#' X <- matrix(rnorm(200 * 10), 200, 10)
#' estimate_m(X, k = 120)
#' estimate_m(X, k = 120, method = "individual", threshold = 0.02)
#' @export
estimate_m <- function(
  sample, k,
  threshold = 0.01,
  max_proportion_m = 1 / 3,
  method = "max",
  parallelize = FALSE
) {
  n <- nrow(sample)

  estimate_m_from_diffs <- function(diffs) {
    for (i in seq_len(length(diffs))) {
      if (diffs[i] < threshold) {
        return(i - 1)
      }
    }
    i
    # which.min(diffs)  # TODO: add option
  }

  if (floor(k * max_proportion_m) > 0) {
    m1 <- estimate_m_from_diffs(
      get_f1(sample, k, max_proportion_m, parallelize)
    )
  } else {
    m1 <- 0
  }
  if (floor((n - k) * max_proportion_m) > 0) {
    m2 <- estimate_m_from_diffs(
      get_f2(sample, k, max_proportion_m, parallelize)
    )
  } else {
    m2 <- 0
  }

  if (method == "max") {
    max(m1, m2)
  } else if (method == "individual") {
    c("m1" = m1, "m2" = m2)
  }
}
