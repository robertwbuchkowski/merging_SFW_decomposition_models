# ============================================================
# PARALLEL BACKEND - run independent scenario spin-ups on several cores.
# ------------------------------------------------------------
# Scenario spin-ups are embarrassingly parallel: each scenario is an
# independent ODE integration that writes its own .rds file, so there is
# nothing to synchronise. This file provides ONE portable helper, par_lapply(),
# that behaves like lapply() but spreads the work over cores:
#
#   * macOS / Linux -> parallel::mclapply() (fork: workers inherit everything,
#                      no setup cost, lowest overhead)
#   * Windows       -> a PSOCK cluster (fork is unavailable, so each worker is a
#                      fresh R session: it must setwd(), load packages and
#                      source() the project files itself -- handled here)
#   * ncores = 1    -> plain lapply(), so everything stays debuggable
#
# Only base R's `parallel` package is needed (ships with R).
# ============================================================
suppressPackageStartupMessages(library(parallel))

# ------------------------------------------------------------
# detect_cores(): how many workers to use by default. Leaves `reserve` cores
# free so the machine stays usable, and never returns less than 1.
# ------------------------------------------------------------
detect_cores <- function(reserve = 1, max_cores = Inf) {
  n <- parallel::detectCores(logical = FALSE)
  if (is.na(n) || n < 1) n <- 1
  max(1, min(n - reserve, max_cores))
}

# project files every worker needs (order matters: setup.R defines model_table)
.default_srcs <- c("R/climate_forcing.R", "R/spinup.R", "R/plot_ode_output.R",
                   "R/setup.R", "R/compare_functions.R", "R/fit_animals.R",
                   "R/dynamic_spinup.R")
.default_pkgs <- c("deSolve", "rootSolve", "yaml")

# ------------------------------------------------------------
# par_lapply(): lapply() over `X`, in parallel.
#   FUN      must be SELF-CONTAINED: take everything it needs as arguments (or
#            via ...), not from the calling environment. On Windows the workers
#            are fresh sessions and will NOT see your global variables.
#   ncores   NULL -> detect_cores(). 1 -> sequential lapply (good for debugging).
#   packages / source_files  loaded on each PSOCK worker.
# Each element is wrapped in tryCatch, so one failing scenario returns an error
# object instead of killing the whole run.
# ------------------------------------------------------------
par_lapply <- function(X, FUN, ...,
                       ncores       = NULL,
                       packages     = .default_pkgs,
                       source_files = .default_srcs,
                       wd           = getwd(),
                       verbose      = TRUE) {

  if (is.null(ncores)) ncores <- detect_cores()
  ncores <- max(1L, min(as.integer(ncores), length(X)))

  safe <- function(x, ...) tryCatch(FUN(x, ...),
                                    error = function(e)
                                      structure(list(item = x, message = conditionMessage(e)),
                                                class = "par_error"))

  if (ncores == 1L) {
    if (verbose) message("par_lapply: 1 core (sequential)")
    return(lapply(X, safe, ...))
  }

  if (.Platform$OS.type != "windows") {
    if (verbose) message("par_lapply: ", ncores, " cores (fork/mclapply)")
    res <- parallel::mclapply(X, safe, ..., mc.cores = ncores,
                              mc.preschedule = FALSE)
  } else {
    if (verbose) message("par_lapply: ", ncores, " cores (PSOCK cluster)")
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # each fresh worker: move to the project root, load packages, source the code
    parallel::clusterCall(cl, function(wd, pkgs, srcs) {
      setwd(wd)
      for (p in pkgs) suppressPackageStartupMessages(
        library(p, character.only = TRUE))
      for (s in srcs) source(s)
      invisible(NULL)
    }, wd, packages, source_files)
    res <- parallel::parLapply(cl, X, safe, ...)
  }

  failed <- vapply(res, inherits, logical(1), "par_error")
  if (any(failed) && verbose)
    for (f in res[failed])
      message("  FAILED: ", f$item, " -- ", f$message)
  res
}

# drop failed elements from a par_lapply() result (and report them)
par_ok <- function(res, verbose = TRUE) {
  failed <- vapply(res, inherits, logical(1), "par_error")
  if (any(failed) && verbose)
    message(sum(failed), " of ", length(res), " task(s) failed; keeping the rest.")
  res[!failed]
}
