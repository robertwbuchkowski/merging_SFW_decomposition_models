# ============================================================
# RENDER SCENARIO REPORTS  -  one self-contained HTML per scenario.
#
# Reads the saved spin-up + follow-up outputs and renders
# Reports/scenario_report.Rmd for each scenario, writing the HTML to Reports/.
# The report recomputes the equilibrium + animal effect from the saved
# parameters, and reads the seasonal spin-up and follow-up trajectories from
# Data/spinup and Data/followup.
#
# HOW TO USE
#   * First produce the saved outputs:
#       Scripts/spinup_dynamic.R      -> Data/spinup/<model>_<scenario>_*.rds
#       Scripts/followup_analysis.R   -> Data/followup/<model>_<scenario>_*.rds
#   * Then run this script from the project root. Needs the `rmarkdown`
#     package (and pandoc, which RStudio bundles).
#   * Edit `scenarios` / `model` below to choose what to render; leaving
#     scenarios = NULL auto-detects every scenario with a saved baseline.
# ============================================================
if (!requireNamespace("rmarkdown", quietly = TRUE))
  stop("Install the 'rmarkdown' package to render reports: install.packages('rmarkdown')")

model     <- "millennial"
scenarios <- NULL                 # NULL = auto-detect from Data/spinup
out_dir   <- "Reports"
rmd       <- file.path(out_dir, "scenario_report.Rmd")

if (is.null(scenarios)) {
  fs <- list.files("Data/spinup", pattern = sprintf("^%s_.*_baseline\\.rds$", model))
  scenarios <- sub(sprintf("^%s_(.*)_baseline\\.rds$", model), "\\1", fs)
}
if (!length(scenarios)) stop("No saved baseline spin-ups found in Data/spinup.")

message("Rendering reports for: ", paste(scenarios, collapse = ", "))
root <- normalizePath(".")

for (sc in scenarios) {
  out_html <- sprintf("%s_%s_report.pdf", model, sc)
  message("  -> ", file.path(out_dir, out_html))
  rmarkdown::render(
    input       = rmd,
    output_file = out_html,
    output_dir  = out_dir,
    params      = list(scenario = sc, model = model, project_root = root),
    envir       = new.env(parent = globalenv()),
    quiet       = TRUE)
}
message("Done. Reports written to ", normalizePath(out_dir))
