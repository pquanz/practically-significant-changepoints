# Overview 

This directory contains the files `statistic.R`, `changepoint.R` and `trimming.R`. The former contains the scripts of the computation of the statistics $T_n$, $V_n$ and $\hat S_n$ used to perform the tests. The second file contains all functions related to the computation of $\hat k_n$ corresponding to the quantity $k_0$. The latter file contains all estimation procedures for the trimming parameter $m$.



## Contents of the file statistic.R

In the following, we remark on some core computational steps.

### Computation of $T_n$

- The function `tn_statistic()` will compute the statistic $T_n$.
  <br>
  Expanding (ignoring the prefactor) yields
  $$
  \begin{aligned}
  T_n (k, m; \lambda) &= \sum_{\substack{i_1, i_2 = 1 \\ |i_1 - i_2| > m}}^{\lfloor \lambda k \rfloor} \sum_{\substack{j_1, j_2 = k + 1 \\ |j_1 - j_2| > m}}^{k + \lfloor \lambda (n-k) \rfloor} (X_{i_1} - X_{j_1})^\top (X_{i_2} - X_{j_2})\\
  &= N_m (\lfloor \lambda (n-k) \rfloor) \cdot T_1 (\lambda) + N_m (\lfloor \lambda k \rfloor) \cdot T_2 (\lambda) - 2 \cdot T_3 (\lambda),
  \end{aligned}
  $$
  where
  $$
  \begin{aligned}
  T_1 (\lambda) &= \sum_{\substack{i_1, i_2 = 1 \\ |i_1 - i_2| > m}}^{\lfloor \lambda k \rfloor} X_{i_1}^\top X_{i_2},\\
  T_2 (\lambda) &= \sum_{\substack{i_1, i_2 = k + 1 \\ |i_1 - i_2| > m}}^{k + \lfloor \lambda (n-k) \rfloor} X_{j_1}^\top X_{j_2},\\
  T_3 (\lambda) &= \sum_{\substack{i_1, i_2 = 1 \\ |i_1 - i_2| > m}}^{\lfloor \lambda k \rfloor} \sum_{\substack{j_1, j_2 = k + 1 \\ |j_1 - j_2| > m}}^{k + \lfloor \lambda (n-k) \rfloor} X_{i_1}^\top X_{j_2}.
  \end{aligned}
  $$
  These terms will be calculated as follows:
    - $T_1 (\lambda)$: Rewrite it as
      $$
      \begin{aligned}
      T_1 (\lambda)
      &= \bigg( X_1^\top (0, ..., 0, X_{m+2}, ..., X_{\lfloor \lambda k \rfloor}) + \cdots + X_{\lfloor \lambda k \rfloor}^\top (X_1, ..., X_{\lfloor \lambda k \rfloor - m - 1}, 0, ..., 0) \bigg) \begin{pmatrix}
      1 \\ \vdots \\ 1
      \end{pmatrix}\\
      &= X_1^\top (X_{j})_{|1-j| > m}  \begin{pmatrix}
      1 \\ \vdots \\ 1
      \end{pmatrix} + \cdots +  X_{\lfloor \lambda k \rfloor}^\top (X_{j})_{|\lfloor \lambda k \rfloor-j| > m}  \begin{pmatrix}
      1 \\ \vdots \\ 1
      \end{pmatrix}
      \end{aligned}
      $$
      The indices in the respective matrices on the right hand sides are stored in the list `indices_m_apart`. Then, the variable `trimmed_sample_t` is the respective matrix on the right hand side in each term. The vector containing ones is used to sum up all entries of the vector of inner products that is left. In the code, this is done by the line `sum(sample[i,] %*% trimmed_sample_t)`. Finally, the outer for loop `for (i in 1:klam) {...}` represents the sum that is calculated over all these terms. The variable `t1` accumulates the running sum and will take the value of $T_1 (\lambda)$ in the end.
    - $T_2 (\lambda)$: Similar to $T_1 (\lambda)$. Its value will be stored in `t2`.

    - $T_3 (\lambda)$: This term will only be nonzero if $\lfloor \lambda k \rfloor > m$ and $\lfloor \lambda (n-k) \rfloor > m$. By Lemma S2.1
      $$
      \begin{aligned}
      \sum_{\substack{i_1, i_2 = 1 \\ |i_1 - i_2| > m}}^{\lfloor \lambda k \rfloor} X_{i_1} = \sum_{i=1}^{\lfloor \lambda k \rfloor} (\lfloor \lambda k \rfloor - 2m - 1) X_i + \sum_{i=1}^{m} i X_{i + \lfloor \lambda k \rfloor - m} - \sum_{i=1}^{m} (i-m-1) X_i,
      \end{aligned}
      $$
      which will be stored in `t31`. The function `colSums(sample[(1:klam), , drop = FALSE])` will sum up all values $X_1, ..., X_{\lfloor \lambda k \rfloor}$. A similar representation can be found for
      $$
      \begin{aligned}
      \sum_{\substack{j_1, j_2 = k+1 \\ |j_1 - j_2| > m}}^{k + \lfloor \lambda (n-k) \rfloor} X_{j_1},
      \end{aligned}
      $$
      whose value will be stored in `t32`. Close to the end, the inner product of those to terms will be calculated by `t31 %*% t32`.

### Choice of normalizers $V_n$

There are many possibilities to compute a self-normalizer, which will work as well, but might yield different limiting distributions. The one used in the paper is type 1 for $T_n$ and type 0 for $\hat \delta_j$.


Let
$$
\Lambda_n (\lambda) = \frac{N_m (\lfloor \lambda \hat k_n \rfloor)}{N_m (\hat k_n)} \frac{N_m (\lfloor \lambda (n - \hat k_n ) \rfloor )}{N_m (n - \hat k_n)}
$$
and by using the parameter `type` in the `normalizer()` function, different self-normalizers can be used. The supported types are

- **type 0** (exported under the alias `vn_statistic()` and used as the normalizing statistic)
  $$
  \bigg( \int_0^1 \big( T_n (\lambda) - \lambda^4 T_n (1)  \big)^2 \mathrm{d} \lambda \bigg)^{1/2}
  $$


- **type 1** (exported under the alias `vl_statistic()` and used in the definition of $\hat S_n$)
  $$
  \bigg( \int_0^1 \big( T_n (\lambda) - \Lambda_n (\lambda) T_n (1)  \big)^2 \mathrm{d} \lambda \bigg)^{1/2}
  $$




### Estimation of $S$

In the estimator $\hat S_n$, you can set the exponent $\kappa$, which results in less conservative estimate of $S$, in the definition
$$
\hat S_n = \{ \ell = 1, \ldots, p \mid \hat \delta_\ell^2 > \hat v_\ell \cdot \log^\kappa (p) \}.
$$
In the manuscript, this quantity is set to $\kappa = 3 / 2$. The computation can be long and it is recommended to call this function with `parallelize = TRUE`.




## Contents of the file changepoint.

This file contains all functions related to the estimation of the change point. The function `estimate_changepoint()` is the main function of this file. Sometimes it is useful to plot the change point function (of which the maximum is taken in the definition of the estimator), which can be done with the function `plot_changepoint_function()`.



## Contents of the file trimming.R

The implementation of the estimator for the parameter $m$ allows for two or a single trimming parameter. Both are handled accordingly in the functions `tn_statistic()`, `vn_statistic()`, `vl_statistic()` and `estimate_support_set()`. The methodology described in the manuscript is used, where the parameter `threshold` refers to the value $T$. Choosing `method = "max"` returns the estimator described in the paper, whereas `method = "individual"` handles the more general case.

For a visiual approach to choosing the estimators, the functions `get_delta_f1()` and `get_delta_f2()` help finding the appropriate value for $\hat m_n$ (or $\hat m_n^{(1)}$ and $\hat m_n^{(2)}$). As these functions are computationally intense, it is recommended to run these in parallel with the parameter `parallelize = TRUE`.
