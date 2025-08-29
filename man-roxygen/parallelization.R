#' @section Parallelization:
#' If `parallelize = TRUE`, this function uses the **future.apply** backend.
#' It runs **sequentially** unless you set a future plan with more than one
#' worker. In your script/session, do for example:
#'
#' \preformatted{
#'   # choose a plan (Windows/macOS/Linux)
#'   future::plan(future::multisession, workers = 4)
#'
#'   # run your code here ...
#'
#'   # return to sequential mode
#'   future::plan(future::sequential)
#' }
#'
#' Reproducibility: set a seed **before** calling (e.g., `set.seed(1)`);
#' internally we use `future.seed = TRUE`. Package code and dependencies are
#' shipped to workers automatically; just ensure any packages used internally
#' are listed in `Imports:`/`Suggests:` in `DESCRIPTION`.
