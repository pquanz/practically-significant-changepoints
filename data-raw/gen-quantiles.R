# remotes::install_github("youruser/pscp@v0.1.0")
library(pscp)


# --- grids (edit as needed) ---
alpha_grid    <- c(
  0.005, 0.01, 0.025, 0.05, 0.10, 0.20,
  0.80, 0.90, 0.95, 0.975, 0.99, 0.995
)
n_points_grid <- c(5:50)
dists         <- c("G", "V", "H", "W")

n_samples <- 10000L
set.seed(101)

# pick sampler by distribution
r_by_dist <- function(n, n_points = 20, dist = "G") {
  switch(dist,
    G = draw_g(n, n_points, 6),
    V = draw_v(n, n_points, 6),
    H = draw_h(n, n_points, 4),
    W = draw_w(n, n_points, 4)
  )
}

# Monte Carlo function for quantile
quant_at <- function(alpha, n_points, dist) {
  x <- r_by_dist(n_samples, n_points, dist)
  stats::quantile(x, probs = 1 - alpha, names = FALSE, type = 8)
}

# build table
quantiles <- expand.grid(
  dist = dists,
  alpha = alpha_grid,
  n_points = n_points_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

quantiles$quantile <- mapply(
  quant_at, quantiles$alpha, quantiles$n_points, quantiles$dist
)

attr(quantiles, "meta") <- list(
  version = "0.0.1",
  n_samples = n_samples,
  seed = 101L,
  date = as.character(Sys.Date())
)

# save as internal data in R/sysdata.rda
usethis::use_data(
  quantiles, internal = TRUE, compress = "xz", overwrite = TRUE
)
