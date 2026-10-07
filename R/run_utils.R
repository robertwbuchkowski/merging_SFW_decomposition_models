#------------------------------------------------------------------------#
# Run helpers: step runner and step-through debugging ----
#------------------------------------------------------------------------#
# Used in: every script in Scripts/ (debugonce_functions(), called right after
# each script's source() block) and Scripts/0_run_all.R (run_one()).

#------------------------------------------------------------------------#
# debugonce_functions() ----
#------------------------------------------------------------------------#
# Used in: Scripts/1_ ... 6_*.R, right after the source() calls at the top.
# Does nothing unless the option sfw.debugonce is TRUE. When it is, every
# function currently defined in `envir` (by default the global environment,
# i.e. everything the script has just sourced from R/) is flagged with
# debugonce(), so running the script top to bottom opens the browser the first
# time each function is called. Turn it on from Scripts/0_run_all.R
# (debug_functions <- TRUE) or, for a single script, run
#   options(sfw.debugonce = TRUE)
# before sourcing it. `exclude` lists function names to skip.
# Note: model evaluations sent to parallel workers are not debugged there, so
# scripts that use workers fall back to serial execution while this is on.
debugonce_functions <- function(envir = globalenv(),
                                exclude = c("debugonce_functions", "run_one")) {
  if (!isTRUE(getOption("sfw.debugonce", FALSE))) return(invisible(character(0)))
  nms <- ls(envir)
  nms <- nms[vapply(nms, function(n) is.function(get(n, envir = envir)), logical(1))]
  nms <- setdiff(nms, exclude)
  for (n in nms) debugonce(get(n, envir = envir))
  message("debugonce set on ", length(nms), " functions: ", paste(nms, collapse = ", "))
  invisible(nms)
}

#------------------------------------------------------------------------#
# run_one() ----
#------------------------------------------------------------------------#
# Used in: Scripts/0_run_all.R (runs each analysis step).
# Sources one script into the GLOBAL environment, as if you had opened and run
# it yourself, so all of its functions and objects can be inspected afterwards.
# If clean = TRUE, the global environment is first cleared (except the objects
# named in `keep`), so each step starts fresh, like a new R session.
run_one <- function(path, clean = TRUE, keep = character(0)) {
  if (clean) {
    drop <- setdiff(ls(globalenv(), all.names = TRUE), keep)
    rm(list = drop, envir = globalenv())
  }
  message("\n", strrep("=", 60), "\n== RUN: ", path, "\n", strrep("=", 60))
  t0 <- Sys.time()
  source(path, local = globalenv())
  message(sprintf("== DONE: %s  (%.1f min)", path,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  invisible(TRUE)
}
