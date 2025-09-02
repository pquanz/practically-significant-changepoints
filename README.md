# practically-significant-changepoints (pscp)


This repository (including the R package pscp) accompanies the publication
"Practically significant change points in high dimension - measuring signal
strength pro active component" (see [arXiv](https://arxiv.org/abs/2508.21520)).



## Install


#### Requirements (and optional packages)

- R (≥ 4.1 recommended)

#### For parallel computing (see usage below)

- Optional: `future`, `future.apply`, `progressr`


```r
# install.packages("remotes")
remotes::install_github("pquanz/practically-significant-changepoints")
```



## Quick Start & Example

```r
library(pscp)

n <- 200
p <- 400
theta <- 0.5    # relative change point

set.seed(101)

# before change point, mean vector (0, ..., 0)
X_prior <- matrix(rnorm(floor(n * theta) * p), floor(n * theta), p)

# after change point, mean vector (1, ..., 1)
X_post <- matrix(rnorm(floor(n * theta) * p, 1), floor(n * theta), p)

# full sample
X <- rbind(X_prior, X_post)

khat <- estimate_changepoint(X)
print(khat)     # 100

# statistics
tn <- tn_statistic(X, khat)
vn <- vn_statistic(X, khat)

q <- quantile_g()   # default alpha = 0.95 and n_points = 20
Δ <- 25


# Test decision (at boundary of H0: 5% rejections on average)
tn > Δ + q * vn   # FALSE
```

The statistics above use `m = c(m1 = 0, m2 = 0)` as trimming parameters as default values.
This will be enough for independent data. 

Instead of estimating `khat`, `tn` and `vn` and checking the inequality `tn > Δ + q * vn`,
we can run the command `cp_test_normalized()`.
We can pass `m = 0` (or `m = c(m1 = 0, m2 = 0)`) as well, otherwise it will estimate
these parameters. Here they turn out to be `m = c(m1 = 1, m2 = 1)` as can be seen in
the output.

It will also give the maximum value for $\Delta$ (`delta_max`) such that the test still
rejects.


```r
# all of the above in short: (delta = Δ)
out <- cp_test_normalized(X, delta = 1, alpha = 0.05)
print(out)
#      tn    vn khat delta_max n_points m1 m2 reject
# 1 0.989 0.001  100      0.97       20  1  1  FALSE
```




Another example demonstrates, how to use parallelization with the `future` package.

```r
library(pscp)

n <- 100
p <- 200
theta <- 0.5

set.seed(33)

# use 10 cores
future::plan(future::multisession, workers = 10)
out <- future.apply::future_lapply(
  seq_len(100),
  function(i) {
    X <- rbind(     # full sample as above
      matrix(rnorm(floor(n * theta) * p), floor(n * theta), p),
      matrix(rnorm(floor(n * theta) * p, 1), floor(n * theta), p)
    )
    cp_test_normalized(X, delta = 1, m = 0, alpha = 0.05)
  },
  future.seed = TRUE
)

# do.call(rbind, out) looks as follows:
#        tn    vn khat delta_max n_points m reject
# 1   0.962 0.006   50     0.844       20 0  FALSE
# 2   0.989 0.002   50     0.951       20 0  FALSE
# 3   1.026 0.004   50     0.941       20 0  FALSE
# ...
# 96  1.044 0.002   50     1.002       20 0   TRUE
# 97  1.005 0.003   50     0.951       20 0  FALSE
# 98  1.023 0.005   50     0.931       20 0  FALSE
# 99  0.977 0.003   50     0.923       20 0  FALSE
# 100 0.981 0.003   50     0.919       20 0  FALSE


# take the average of all rejections
print(mean(do.call(rbind, out)$reject))    # 0.04
```
