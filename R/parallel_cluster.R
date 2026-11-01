#------------------------------------------------------------------------#
# Parallel worker pool ----
#------------------------------------------------------------------------#
# Used in: Scripts/4_uncertainty_animal_effect.R, Scripts/5_parameter_identifiability.R,
#   Scripts/6_sobol_variance_partition.R

#------------------------------------------------------------------------#
# model_cluster() ----
#------------------------------------------------------------------------#
# model_cluster(): start a parallel cluster whose workers can run the model.
# Works on Windows, macOS and Linux (PSOCK sockets). Each worker sets the
# project root as its working directory, loads `packages`, sources
# `source_files`, then receives the objects named in `export` from `envir`.
# Returns NULL when n_cores <= 1 (serial). Stop it with parallel::stopCluster().
model_cluster <- function(n_cores, source_files, packages = c("deSolve", "rootSolve", "yaml"),
                          export = character(0), envir = globalenv()) {
  if (is.null(n_cores) || n_cores <= 1) return(NULL)
  cl <- parallel::makeCluster(n_cores)
  wd <- getwd()
  parallel::clusterCall(cl, function(wd, pkgs, files) {
    setwd(wd)
    for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
    for (f in files) source(f)
    NULL
  }, wd, packages, source_files)
  if (length(export)) parallel::clusterExport(cl, export, envir = envir)
  cl
}
