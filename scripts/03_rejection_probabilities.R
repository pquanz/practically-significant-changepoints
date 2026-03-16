#!/usr/bin/env Rscript
#
# 03_rejection_probabilities.R
#
# Reproduces Table 5 from the paper:
#   Rejection probabilities of tests (2.14) and (3.6) at the boundary
#   of the hypotheses and under an alternative.
#
# Usage:
#   Rscript 03_rejection_probabilities.R
#
#   For parallel execution (e.g. 10 cores):
#   PSCP_WORKERS=10 Rscript 03_rejection_probabilities.R
#
# Dependencies:
#   - pscp

library(pscp)

# Suppress per-call progress bars from pscp internals.
# We parallelize at the replication level instead.
if (requireNamespace("progressr", quietly = TRUE)) {
  progressr::handlers("void")
}

N_WORKERS <- as.integer(Sys.getenv("PSCP_WORKERS", "1"))
if (N_WORKERS > 1 && requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)) {
  future::plan(future::multisession, workers = N_WORKERS)
  cat(sprintf("Using %d workers for replication-level parallelism.\n", N_WORKERS))
  USE_FUTURE <- TRUE
} else {
  USE_FUTURE <- FALSE
}


# ============================================================================
# Data generating processes (Section 5.1)
# ============================================================================

#' Generate data from model (2.1)
#'
#' @param n Sample size.
#' @param p Dimension.
#' @param k0 True change point location.
#' @param delta Vector of length p (change in mean after k0).
#' @param mu Common mean vector of length p.
#' @param model One of "IND", "MA2", "MA4", "MA6", "AR0.5", "AR0.6".
#' @return n x p matrix.
generate_data <- function(n, p, k0, delta, mu, model = "IND") {

  ma_coeffs <- c(0.5, 0.25, 0.2, 0.1, 0.05, 0.025)

  # Covariance matrix Sigma for epsilon_tilde
  # Sigma_ij = 0.5 * 0.9^|i-j|
  make_sigma_chol <- function(p) {
    idx <- 0:(p - 1)
    Sigma <- 0.5 * outer(idx, idx, function(i, j) 0.9^abs(i - j))
    chol(Sigma)  # upper Cholesky: Sigma = t(L) %*% L
  }

  if (model == "IND") {
    # {epsilon_j} iid ~ N(0, diag_p(0.5))
    errors <- matrix(rnorm(n * p, sd = sqrt(0.5)), n, p)

  } else if (model %in% c("MA2", "MA4", "MA6")) {
    q <- as.integer(sub("MA", "", model))
    L <- make_sigma_chol(p)

    # generate (n + q) iid innovations, then form MA(q)
    n_total <- n + q
    eps <- matrix(rnorm(n_total * p), n_total, p) %*% L
    errors <- matrix(0, n, p)
    for (t in seq_len(n)) {
      errors[t, ] <- eps[t + q, ]
      for (lag in seq_len(q)) {
        errors[t, ] <- errors[t, ] + ma_coeffs[lag] * eps[t + q - lag, ]
      }
    }

  } else if (model %in% c("AR0.5", "AR0.6")) {
    c_ar <- as.numeric(sub("AR", "", model))
    L <- make_sigma_chol(p)

    burn <- 200
    n_total <- n + burn
    eps <- matrix(rnorm(n_total * p), n_total, p) %*% L
    errors <- matrix(0, n_total, p)
    errors[1, ] <- eps[1, ]
    for (t in 2:n_total) {
      errors[t, ] <- c_ar * errors[t - 1, ] + eps[t, ]
    }
    errors <- errors[(burn + 1):n_total, , drop = FALSE]

  } else {
    stop("Unknown model: ", model)
  }

  # X_j = mu + eta_j  (j <= k0)  or  mu + rho_j + delta  (j > k0)
  # since f1 = f2, we use the same error process for both parts
  X <- sweep(errors, 2, mu, "+")
  if (k0 < n) {
    X[(k0 + 1):n, ] <- sweep(X[(k0 + 1):n, , drop = FALSE], 2, delta, "+")
  }
  X
}


# ============================================================================
# Build the signal vector delta
# ============================================================================

#' Construct delta = (d, ..., d, 0, ..., 0) with s = round(s_ratio * p)
#' nonzero entries.
#'
#' For type = "sparse":  ||delta||^2_{2,0} = d^2 = norm_val
#' For type = "l2":      ||delta||^2 = s * d^2 / p = norm_val
make_delta <- function(p, s_ratio, norm_val, type = "sparse") {
  s <- max(1, round(s_ratio * p))
  if (type == "sparse") {
    d <- sqrt(norm_val)
  } else if (type == "l2") {
    d <- sqrt(norm_val * p / s)
  }
  c(rep(d, s), rep(0, p - s))
}


# ============================================================================
# Single-replication test wrappers (using cp_test_* from pscp)
# ============================================================================

#' Run test (2.14) on one dataset
#'
#' @return Logical: TRUE if H0 is rejected.
run_test_l2 <- function(X, Delta, alpha, K, method_m) {
  out <- cp_test_normalized(
    X, delta = Delta, alpha = alpha, method_m = method_m, n_points = K
  )
  out$reject
}


#' Run test (3.6) on one dataset
#'
#' @return Logical: TRUE if H0 is rejected.
run_test_sparse <- function(X, Delta, alpha, K, method_m) {
  out <- cp_test_sparsity_adj(
    X, delta = Delta, alpha = alpha, method_m = method_m, n_points = K
  )
  out$reject
}


# ============================================================================
# Main simulation (Table 5)
# ============================================================================

run_simulation <- function(
  n_rep    = 1000,
  models   = c("IND", "MA2", "MA6", "AR0.5", "AR0.6"),
  ns       = c(200, 400),
  ps       = c(200, 400, 800),
  s_ratios = c(0.25, 0.50, 0.75, 1.00),
  Delta    = 2.0,
  H0_norm  = 2.0,      # boundary of H0
  H1_norm  = 2.25,     # alternative
  alpha    = 0.05,
  K        = 20,
  method_m = "individual",
  theta_0  = 0.6,
  seed     = 42
) {

  set.seed(seed)
  mu_val <- 10

  results <- list()

  for (model in models) {
    for (n in ns) {
      for (p in ps) {

        k0 <- floor(n * theta_0)
        mu <- rep(mu_val, p)

        cat(sprintf("\n=== Model: %s, n=%d, p=%d ===\n", model, n, p))

        # --- Test (2.14): normalized l2 norm (does not depend on s/p) ---
        cat("  Test (2.14)...\n")

        delta_l2_H0 <- make_delta(p, 1.0, H0_norm, type = "l2")
        delta_l2_H1 <- make_delta(p, 1.0, H1_norm, type = "l2")

        rej_H0 <- 0L
        rej_H1 <- 0L

        run_one_l2 <- function(r) {
          X0 <- generate_data(n, p, k0, delta_l2_H0, mu, model)
          r0 <- run_test_l2(X0, Delta, alpha, K, method_m)
          X1 <- generate_data(n, p, k0, delta_l2_H1, mu, model)
          r1 <- run_test_l2(X1, Delta, alpha, K, method_m)
          c(H0 = r0, H1 = r1)
        }

        if (USE_FUTURE) {
          res_l2 <- do.call(rbind, future.apply::future_lapply(
            seq_len(n_rep), run_one_l2, future.seed = TRUE
          ))
        } else {
          res_l2 <- do.call(rbind, lapply(seq_len(n_rep), run_one_l2))
        }
        rej_H0 <- sum(res_l2[, "H0"])
        rej_H1 <- sum(res_l2[, "H1"])
        cat(sprintf("    done. H0=%.3f, H1=%.3f\n", rej_H0 / n_rep, rej_H1 / n_rep))

        results[[length(results) + 1]] <- data.frame(
          model = model, n = n, p = p, test = "(2.14)",
          s_ratio = 1.0,
          H0_rej = rej_H0 / n_rep,
          H1_rej = rej_H1 / n_rep
        )

        # --- Test (3.6): sparsity-adjusted norm (for each s/p) ---
        for (s_ratio in s_ratios) {
          cat(sprintf("  Test (3.6), s/p=%.2f...\n", s_ratio))

          delta_sp_H0 <- make_delta(p, s_ratio, H0_norm, type = "sparse")
          delta_sp_H1 <- make_delta(p, s_ratio, H1_norm, type = "sparse")

          rej_H0 <- 0L
          rej_H1 <- 0L

          run_one_sp <- function(r) {
            X0 <- generate_data(n, p, k0, delta_sp_H0, mu, model)
            r0 <- run_test_sparse(X0, Delta, alpha, K, method_m)
            X1 <- generate_data(n, p, k0, delta_sp_H1, mu, model)
            r1 <- run_test_sparse(X1, Delta, alpha, K, method_m)
            c(H0 = r0, H1 = r1)
          }

          if (USE_FUTURE) {
            res_sp <- do.call(rbind, future.apply::future_lapply(
              seq_len(n_rep), run_one_sp, future.seed = TRUE
            ))
          } else {
            res_sp <- do.call(rbind, lapply(seq_len(n_rep), run_one_sp))
          }
          rej_H0 <- sum(res_sp[, "H0"])
          rej_H1 <- sum(res_sp[, "H1"])
          cat(sprintf("    done. H0=%.3f, H1=%.3f\n", rej_H0 / n_rep, rej_H1 / n_rep))

          results[[length(results) + 1]] <- data.frame(
            model = model, n = n, p = p, test = "(3.6)",
            s_ratio = s_ratio,
            H0_rej = rej_H0 / n_rep,
            H1_rej = rej_H1 / n_rep
          )
        }
      }
    }
  }

  do.call(rbind, results)
}


# ============================================================================
# Run and save
# ============================================================================

if (!dir.exists("results")) dir.create("results")

cat("Starting simulation (Table 5)...\n")

res <- run_simulation(n_rep = 1000)

saveRDS(res, file = "results/table5_rejection_probs.rds")
write.csv(res, file = "results/table5_rejection_probs.csv", row.names = FALSE)

cat("\nDone. Results saved to results/\n")
print(res)