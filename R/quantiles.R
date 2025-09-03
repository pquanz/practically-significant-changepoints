#' Critical value for self-normalization (G/V)
#'
#' Lookup by distribution, level `alpha`, and parameter `s`. If `s` is
#' off-grid and `interpolate = TRUE`, linear interpolation in `s` is used.
#'
#' @param dist "G", "V", "H", "W".
#' @param alpha Numeric in (0, 1); must match a tabulated level.
#' @param n_points Numeric; interpolated in `n_points` if off-grid when
#'    `interpolate = TRUE`.
#' @param interpolate Logical; default TRUE.
#' @param extrapolate Logical; allow linear extrapolation in `n_points`;
#'    default FALSE.
#' @param ... For wrappers [q_g()], [quantile_g()], ..., [quantile_w()] to pass
#'    the arguments `interpolate` and `extrapolate`.
#' @return Numeric scalar critical value (or `NA_real_`).
quantile_function <- function(dist = c("G", "V", "H", "W"),
                              alpha = 0.95,
                              n_points = 20,
                              interpolate = TRUE,
                              extrapolate = FALSE) {
  dist <- match.arg(dist)
  tab <- quantiles

  alphas <- sort(unique(tab$alpha))
  if (!any(abs(alphas - alpha) <= 1e-12)) {
    stop("`alpha` must be one of: ", paste(format(alphas), collapse = ", "))
  }

  subs <- tab[tab$dist == dist & tab$alpha == alpha,
              c("n_points", "quantile"), drop = FALSE]
  if (!nrow(subs)) return(NA_real_)

  hit <- which(abs(subs$n_points - n_points) <= 1e-12)
  if (length(hit)) return(subs$quantile[hit[1L]])

  if (!interpolate) return(NA_real_)

  o <- order(subs$n_points)
  subs <- subs[o, ]
  rule <- if (extrapolate) 2L else 1L
  as.numeric(
    stats::approx(subs$n_points, subs$quantile, xout = n_points, rule = rule)$y
  )
}


#' @rdname quantile_function
#' @export
q_g <- function(alpha = 0.95, n_points = 20, ...) {
  quantile_function("G", alpha = alpha, n_points = n_points, ...)
}
#' @rdname quantile_function
#' @export
quantile_g <- function(alpha = 0.95, n_points = 20, ...) {
  quantile_function("G", alpha = alpha, n_points = n_points, ...)
}

#' @rdname quantile_function
#' @export
q_v <- function(alpha = 0.95, n_points = 20, ...) {
  quantile_function("V", alpha = alpha, n_points = n_points, ...)
}
#' @rdname quantile_function
#' @export
quantile_v <- function(alpha = 0.95, n_points = 20, ...) {
  quantile_function("V", alpha = alpha, n_points = n_points, ...)
}

#' @rdname quantile_function
#' @export
q_h <- function(alpha = 0.95, n_points = 20, ...) {
  quantile_function("H", alpha = alpha, n_points = n_points, ...)
}
#' @rdname quantile_function
#' @export
quantile_h <- function(alpha = 0.95, n_points = 20, ...) {
  quantile_function("H", alpha = alpha, n_points = n_points, ...)
}

#' @rdname quantile_function
#' @export
q_w <- function(alpha = 0.95, n_points = 20, ...) {
  quantile_function("W", alpha = alpha, n_points = n_points, ...)
}
#' @rdname quantile_function
#' @export
quantile_w <- function(alpha = 0.95, n_points = 20, ...) {
  quantile_function("W", alpha = alpha, n_points = n_points, ...)
}
