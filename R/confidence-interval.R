#' (Upper) one-sided (1 - `alpha`) confidence interval for
#' the (squared) normalized norm
#'
#' @param test_result data frame that can be obtained by applying
#'   [cp_test_normalized()] or [cp_test_sparsity_adj()]
#' @param alpha Confidence level in (0,1). Default 0.05.
#' @return Named numeric vector: c(lower, upper)
#' @export
ci_sq_norm <- function(test_result, alpha = 0.05) {
  q <- quantile_g(1 - alpha,test_result$n_points)
  c(lower = 0, upper = round(test_result$tn + q * test_result$vn, digits = 2))
}

#' Two-sided (1 - `alpha`) confidence interval for
#' the (squared) normalized norm
#'
#' @param test_result data frame that can be obtained by applying
#'   [cp_test_normalized()] or [cp_test_sparsity_adj()]
#' @param alpha Confidence level in (0,1). Default 0.05.
#' @return Named numeric vector: c(lower, upper)
#' @export
ci_sq_norm_twoside <- function(test_result, alpha = 0.05) {
  q <- quantile_g(1 - alpha / 2, test_result$n_points)
  c(
    lower = round(max(0, test_result$tn - q * test_result$vn), digits = 2),
    upper = round(test_result$tn + q * test_result$vn, digits = 2)
  )
}