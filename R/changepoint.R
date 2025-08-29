#' Helper function for the change point estimation
#'
#' This function will be used to calculate and/or plot the
#' change point function.
#'
#' @param sample Numeric matrix or numeric data frame (rows = n, cols = p).
#' @param dist_to_margin Integer index (1-based). Bounds the domain
#'   of this function away from `1` and `nrow(sample)`, where uncontrollable
#'   spikes of the function values might occur.
#'
#' @seealso [estimate_changepoint()], [plot_changepoint_function()]
changepoint_function <- function(sample, dist_to_margin = 1) {
  n <- nrow(sample)
  cp_at_k <- function(k) {
    mean1 <- colMeans(sample[1:k, , drop = FALSE])
    mean2 <- colMeans(sample[(k + 1):n, , drop = FALSE])
    meandiff <- mean1 - mean2
    k * (n - k) / n^2 * (meandiff %*% meandiff)
  }
  lapply((dist_to_margin + 1):(n - dist_to_margin), cp_at_k)
}

#' Estimation of k hat (Change point location)
#'
#' Computes the estimator for the change point as described
#' in the manuscript.
#'
#' @param sample Numeric matrix or numeric data frame (rows = n, cols = p).
#' @param dist_to_margin Integer index (1-based). Bounds the possible
#'   values for khat away from `1` and `nrow(sample)`, where uncontrollable
#'   spikes might occur.
#'
#' @return Numeric khat: a value between `1 + dist_to_margin` and
#'   `nrow(sample) - dist_to_margin`.
#'
#' @seealso [plot_changepoint_function()], [changepoint_function()]
#'
#' @examples
#' X <- matrix(rnorm(200 * 20), 200, 20)
#' khat <- estimate_changepoint(X)
#' khat_safe <- estimate_changepoint(X, 5)
#'
#' @export
estimate_changepoint <- function(sample, dist_to_margin = 1) {
  vals <- changepoint_function(sample, dist_to_margin)
  return(which.max(vals) + dist_to_margin)
}

#' Estimation of k hat (Change point location)
#'
#' Computes the estimator for the change point as described
#' in the manuscript.
#'
#' @param sample Numeric matrix or numeric data frame (rows = n, cols = p).
#'
#' @return Invisibly returns `NULL`.
#'
#' @seealso [changepoint_function()]
#'
#' @examples
#' X <- matrix(rnorm(200 * 20), 200, 20)
#' plot_changepoint_function(X)
#'
#' @export
plot_changepoint_function <- function(sample) {
  plot(unlist(changepoint_function(sample)), type = "l")
  invisible(NULL)
}