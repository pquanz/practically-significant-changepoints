# practically-significant-changepoints (pscp)


This repository (including the R package pscp) accompanies the publication
"Practically significant change points in high dimension - measuring signal
strength pro active component."



## Install


#### Requirements (and optional packages)

- R (≥ 4.1 recommended)

#### For parallel computing (see usage below)

- Optional: `future`, `future.apply`, `progressr`
- For building from source on Windows: **Rtools**; on macOS: **Xcode Command Line Tools**


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

# after change point, mean vector (5, ..., 5)
X_post <- matrix(rnorm(floor(n * theta) * p, 5), floor(n * theta), p)

# full sample
X <- rbind(X_prior, X_post)

khat <- estimate_changepoint(X)
print(khat)     # 100

# statistics
tn <- tn_statistic(X, khat)
vn <- vn_statistic(X, khat)

q <- quantile_g(0.05, 20)
Δ <- 25

# Test decision (at boundary of H0: 5% rejections on average)
print(tn > Δ + q * vn)  # FALSE
```

All of the above can be executed in short by using the command `cp_test()`.
It will also give the maximum value for $\Delta$ such that the test still rejects.

```r
# all of the above in short: (delta = Δ)
out <- cp_test(X, delta = 25, adaptive = FALSE, alpha = 0.05)
print(out)
#       tn    vn khat reject max_delta_reject
# 1 24.951 0.005  100  FALSE           24.863
```




Another example demonstrates, how to use parallelization with the `future` package.

```r
library(pscp)

n <- 100
p <- 200
theta <- 0.5

set.seed(101)

# use 10 cores
future::plan(future::multisession, workers = 10)
out <- future.apply::future_lapply(
  seq_len(100),
  function(i) {
    X <- rbind(     # full sample as above
      matrix(rnorm(floor(n * theta) * p), floor(n * theta), p),
      matrix(rnorm(floor(n * theta) * p, 1), floor(n * theta), p)
    )
    cp_test(X, delta = 1, adaptive = FALSE, alpha = 0.05)
  },
  future.seed = TRUE
)

# do.call(rbind, out) looks as follows:
#        tn    vn khat reject max_delta_reject
# 1   0.986 0.002   50  FALSE            0.938
# 2   0.995 0.002   50  FALSE            0.957
# 3   0.976 0.001   50  FALSE            0.954
# 4   1.004 0.002   50  FALSE            0.966
# 5   1.018 0.003   50  FALSE            0.967
# 6   0.991 0.002   50  FALSE            0.944
# 7   0.956 0.003   50  FALSE            0.898
# 8   1.008 0.002   50  FALSE            0.979
# 9   1.032 0.001   50   TRUE            1.006
# 10  0.997 0.002   50  FALSE            0.953
# ...
# 100 1.040 0.004   50  FALSE            0.955


# take the average of all rejections
print(mean(do.call(rbind, out)$reject))    # 0.05
```
