# pscp 0.2.1


## Changed

- Removed alias `s_hat()`
- Changed order of parameter `method` in `estimate_m()`.


## Fixes

- Fixed `estimate_m()`, where `get_f1()` was used, instead of `get_delta_f1()`.



# pscp 0.2.0


## Added
- confidence intervals: `ci_sq_norm()` and `ci_sq_norm_twoside()`. 


## Changed
- testing withouth setting `delta`.
- minor fix in `trimming.R`.



# pscp 0.1.0

First public, complete release.

## Added
- Core statistic: `tn_statistic()` and `vn_statistic()`.
- Estimators for parameters: `estimate_support_set()`, `estimate_changepoint()`,
  `estimate_m()`, 
- For estimation of the trimming parameter $m$, we provide `get_delta_f1()`
  and `get_delta_f2()` for visual inspection.
- Tests for relevant hypotheses: `cp_test_normalized()` and `cp_test_sparsity_adj()`
  as described in .
- Finite-sample quantiles (as internal dataset) with accessors:
  E.g. `q_g()`, `quantile_g()`, as well as those for distributions
  $\mathbb{V}, \mathbb{H}$ and $\mathbb{W}$.
- Optional parallelization via **future.apply** (`parallelize = TRUE`).
- Script in data-raw to regenerate quantiles.
- Basic documentation.

## Notes
- No datasets are bundled in the package. It focuses on methods for the paper.
