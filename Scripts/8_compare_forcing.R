#------------------------------------------------------------------------#
# COMPARE RESULTS UNDER SEASONAL VS CONSTANT FORCING ----
#------------------------------------------------------------------------#
# Reads the saved outputs of Scripts/1, 4, 5 and 7 from both results folders
# and puts them side by side. Nothing is re-simulated.
#   Results/                   seasonal forcing (annual limit cycle)
#   Results/constant_forcing/  constant mean forcing
# (results_dir_by_forcing in R/sensitivity_settings.R.)
#
# HOW TO PRODUCE THE TWO RUNS (Scripts/0_run_all.R):
#   1. forcing_run <- "seasonal": steps 1, 4, 5, 7 (your existing run)
#   2. forcing_run <- "constant": steps 1, 4, 5, 7 (steps 2, 3, 6 are skipped)
#   3. run_step8 <- TRUE (this script; forcing_run does not matter)
#
# Sections whose files are missing in either folder are skipped with a message.
# Outputs (Results/forcing_comparison/):
#   forcing_comparison.pdf                 all figures, one per page
#   overview.csv                           one row per scenario, key numbers
#   calibration.csv          + calibration.png       fitted animal parameters
#   headline_effect.csv      + headline_effect.png   effect at defaults and 95% interval
#   headline_acceptance.csv  + headline_samples.png  accepted-sample distributions
#   morris[_standardized].csv, morris[_standardized]_parameters.csv,
#     morris[_standardized].png            ranking agreement (Spearman, top-k overlap)
#   identifiability.csv, identifiability_animal_status.csv,
#     identifiability_counts.png, identifiability_animal.png
#
# HELPER FUNCTIONS live in R/compare_forcing.R (sourced below).
library(pacman); p_load(tidyverse)
source("R/compare_functions.R")       # pretty_scenario(), scenario_labeller()
source("R/sensitivity_settings.R")    # results_dir_by_forcing, param_labels
source("R/sensitivity_helpers.R")     # pretty_param()
source("R/compare_forcing.R")         # cmp_*() comparison functions
source("R/run_utils.R")
debugonce_functions()   # only active when options(sfw.debugonce = TRUE); see Scripts/0_run_all.R

out_dir <- "Results/forcing_comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
morris_n_top <- 10        # parameters per scenario in the Morris figure
morris_top_k <- 5         # size of the top set compared between forcings

cat("Comparing", paste(sprintf("%s (%s)", names(results_dir_by_forcing),
                                results_dir_by_forcing), collapse = " vs "), "\n")

#------------------------------------------------------------------------#
# Run every comparison ----
#------------------------------------------------------------------------#
res <- list()
cat("Calibration (Scripts/1):\n");       res$calibration <- cmp_calibration()
cat("Headline effect (Scripts/5):\n");   res$headline    <- cmp_headline()
                                         res$samples     <- cmp_headline_samples()
cat("Morris (Scripts/4):\n");            res$morris      <- cmp_morris("", morris_n_top, morris_top_k)
                                         res$morris_std  <- cmp_morris("_standardized", morris_n_top, morris_top_k)
cat("Identifiability (Scripts/7):\n");   res$ident       <- cmp_identifiability()

#------------------------------------------------------------------------#
# Tables ----
#------------------------------------------------------------------------#
save_tab <- function(x, file) if (!is.null(x)) write_csv(x, file.path(out_dir, file))
save_tab(res$calibration$table, "calibration.csv")
save_tab(res$headline$table,    "headline_effect.csv")
save_tab(res$samples$table,     "headline_acceptance.csv")
save_tab(res$morris$table,      "morris.csv")
save_tab(res$morris$params,     "morris_parameters.csv")
save_tab(res$morris_std$table,  "morris_standardized.csv")
save_tab(res$morris_std$params, "morris_standardized_parameters.csv")
save_tab(res$ident$table,       "identifiability.csv")
save_tab(res$ident$params,      "identifiability_animal_status.csv")
overview <- cmp_overview(res)
save_tab(overview, "overview.csv")

#------------------------------------------------------------------------#
# Figures ----
#------------------------------------------------------------------------#
figs <- list(calibration           = list(res$calibration$plot, 10, 5),
             headline_effect       = list(res$headline$plot,    10, 5),
             headline_samples      = list(res$samples$plot,     10, 7),
             morris                = list(res$morris$plot,      12, 9),
             morris_standardized   = list(res$morris_std$plot,  12, 9),
             identifiability_counts = list(res$ident$plot_counts, 12, 6),
             identifiability_animal = list(res$ident$plot_animal, 12, 9))
figs <- figs[!vapply(figs, function(f) is.null(f[[1]]), logical(1))]
for (nm in names(figs))
  ggsave(file.path(out_dir, paste0(nm, ".png")), figs[[nm]][[1]],
         width = figs[[nm]][[2]], height = figs[[nm]][[3]], dpi = 200, bg = "white")
if (length(figs)) {
  pdf(file.path(out_dir, "forcing_comparison.pdf"), width = 12, height = 8.5)
  for (f in figs) print(f[[1]])
  dev.off()
}

if (!is.null(overview)) {
  cat("\nOverview (seasonal vs constant):\n")
  print(as.data.frame(overview), row.names = FALSE, digits = 3)
}
cat("\nWrote", length(figs), "figures and the tables to", out_dir, "\n")
