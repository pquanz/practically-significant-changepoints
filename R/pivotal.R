#' Simulate a Brownian motion path
#'
#' Generates one path of a standard Brownian motion returned as a list with
#' n_points + 1 elements.
#'
#' @param n_points Integer, number of support points. The path has
#'  `n_points + 1` points including time 0.
#'
#' @return A numeric vector of length `n_points + 1`.
#'
draw_brownian_motion <- function(n_points = 500) {
  c(0, cumsum(stats::rnorm(n_points - 1, 0, sqrt(1 / n_points))))
}

#' Draw i.i.d. samples from the distribution \eqn{G_\alpha}
#'
#' @param n Integer number of draws (>= 1).
#' @param n_points Number of support points of the probability
#'   measure nu.
#' @param alpha Exponent in the integral, see manuscript.
#' @return Numeric vector of length `n`.
#'
#' @examples
#' draw_g(5, alpha = 4)
#'
#' @export
draw_g <- function(n = 1, n_points = 20, alpha = 6) {
  g <- replicate(n, {
    bm <- draw_brownian_motion()
    lambdas <- 1:(n_points - 1) / n_points

    integrand <- function(lam) {
      lam^alpha * (bm[lam * length(bm)] - lam * bm[length(bm)])^2
    }
    v <- sqrt(mean(sapply(lambdas, integrand, simplify = TRUE)))
    bm[length(bm)] / v
  })
  g
}

#' Draw i.i.d. samples from the distribution \eqn{V_\alpha}
#'
#' @param n Integer number of draws (>= 1).
#' @param n_points Number of support points of the probability
#'   measure nu.
#' @param alpha Exponent in the integral, see manuscript.
#' @return Numeric vector of length `n`.
#'
#' @examples
#' draw_v(5, alpha = 4)
#'
#' @export
draw_v <- function(n = 1, n_points = 20, alpha = 6) {
  v <- replicate(n, {
    bm <- draw_brownian_motion()
    lambdas <- 1:(n_points - 1) / n_points

    integrand <- function(lam) {
      lam^alpha * (bm[lam * length(bm)] - lam * bm[length(bm)])^2
    }
    sqrt(mean(sapply(lambdas, integrand, simplify = TRUE)))
  })
  v
}


#' Draw i.i.d. samples from the distribution \eqn{H_\alpha}
#'
#' @param n Integer number of draws (>= 1).
#' @param n_points Number of support points of the probability
#'   measure nu.
#' @param alpha Exponent in the integral, see manuscript.
#' @return Numeric vector of length `n`.
#'
#' @examples
#' draw_h(5, alpha = 4)
#'
#' @export
draw_h <- function(n = 1, n_points = 20, alpha = 4) {
  h <- replicate(n, {
    bm <- draw_brownian_motion()
    lambdas <- 1:(n_points - 1) / n_points

    integrand <- function(lam) {
      lam^alpha * (bm[lam * length(bm)] - lam^2 * bm[length(bm)])^2
    }
    w <- sqrt(mean(sapply(lambdas, integrand, simplify = TRUE)))
    (bm[length(bm)]^2 - 1) / w
  })
  h
}


#' Draw i.i.d. samples from the distribution \eqn{W_\alpha}
#'
#' @param n Integer number of draws (>= 1).
#' @param n_points Number of support points of the probability
#'   measure nu.
#' @param alpha Exponent in the integral, see manuscript.
#' @return Numeric vector of length `n`.
#'
#' @examples
#' draw_w(5, alpha = 4)
#'
#' @export
draw_w <- function(n = 1, n_points = 20, alpha = 4) {
  w <- replicate(n, {
    bm <- draw_brownian_motion()
    lambdas <- 1:(n_points - 1) / n_points

    integrand <- function(lam) {
      lam^alpha * (
        (bm[lam * length(bm)]^2 - lam) - lam^2 * (bm[length(bm)]^2 - 1)
      )^2
    }
    sqrt(mean(sapply(lambdas, integrand, simplify = TRUE)))
  })
  w
}


# # into other file
# get_quantiles <- function(
#   sample_fnc = draw_g,
#   quantiles = c(0.8, 0.90, 0.95, 0.975, 0.99),
#   n_runs = 50000,
#   n_points = 20
# ) {

#   vals <- c()
#   for (i in 1:n_runs) {
#     vals <- c(vals, sample_fnc(n_points = n_points))
#   }

#   mean(vals)
#   quantile(vals, c(0.8, 0.90, 0.95, 0.975, 0.99, 0.995))
# }


# qs <- get_quantiles(random.g, n_points = 10)
# print(qs)

# qs <- get_quantiles(random.g, n_points = 15)
# print(qs)

# qs <- get_quantiles(random.g, n_points = 20)
# print(qs)

# qs <- get_quantiles(random.g, n_points = 25)
# print(qs)


# get_quantiles(random.v)
# get_quantiles(random.h)
# get_quantiles(random.w)