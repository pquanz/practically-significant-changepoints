#'
#' @keywords internal
#' @noRd
.map_apply <- function(index, fn, parallelize = TRUE, combine = "c") {
  parallelize   <- isTRUE(parallelize)
  progress      <- interactive()

  if (isTRUE(progress)) progressr::handlers(
    Sys.getenv("PSCP_PROGRESS_TYPE", "txtprogressbar")
  )

  use_future <- parallelize &&
    requireNamespace("future.apply", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE) &&
    isTRUE(tryCatch(future::nbrOfWorkers() > 1L, error = function(e) FALSE))

  use_progress <- progress &&
    requireNamespace("progressr", quietly = TRUE)

  if (use_future) {
    if (use_progress) {
      progressr::with_progress({
        p <- progressr::progressor(along = index)
        res <- future.apply::future_lapply(
          index,
          function(i) { on.exit(p(), add = TRUE); fn(i) },
          future.seed = TRUE
        )
      })
    } else {
      res <- future.apply::future_lapply(index, fn, future.seed = TRUE)
    }
  } else {
    if (use_progress) {
      progressr::with_progress({
        p <- progressr::progressor(along = index)
        res <- lapply(index, function(i) { on.exit(p(), add = TRUE); fn(i) })
      })
    } else {
      res <- lapply(index, fn)
    }
  }


  combine_fun <- switch(
    combine,
    c     = base::c,
    rbind = base::rbind,
    cbind = base::cbind
  )
  do.call(combine_fun, res)
}


#' @keywords internal
#' @noRd
.progress_bar_default <- function() {
  val <- Sys.getenv("PSCP_PROGRESS", "")
  if (!nzchar(val)) return(FALSE)
  v <- tolower(trimws(val))
  truthy <- c("1", "true", "t", "yes", "y", "on")
  falsy  <- c("0", "false", "f", "no", "n", "off")
  if (v %in% truthy) return(TRUE)
  if (v %in% falsy) return(FALSE)
  if (v == "auto") return(interactive())
  FALSE
}


.onLoad <- function(libname, pkgname) {
  h <- Sys.getenv("PSCP_PROGRESS_TYPE", "cli")
  if (nzchar(h) && requireNamespace("progressr", quietly = TRUE)) {
    # try to set a handler (e.g., "txtprogressbar", "progress", "cli")
    try({
      # progressr::handlers(reset = TRUE)
      # progressr::handlers("cli")
      # progressr::handlers(progressr::handler_cli(interval = 0.3))
      progressr::handlers(progressr::handler_txtprogressbar())
    }, silent = FALSE)
  }
}
