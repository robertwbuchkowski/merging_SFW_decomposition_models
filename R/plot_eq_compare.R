# ============================================================
# plot_eq_compare.R - compare equilibrium pools across models with ggplot
# ------------------------------------------------------------
# Pools differ across models AND span very different magnitudes (C_wood_tree
# ~ thousands vs MIC ~ tens), so the default view FACETS BY POOL with FREE
# axes: every pool gets its own panel and scale, and each panel shows one bar
# per model. Pools that exist in only some models simply appear with fewer
# bars. Needs only ggplot2.
# ============================================================

plot_eq_compare <- function(compare_list,
                            metric      = c("percent_change", "difference",
                                            "treatment", "baseline"),
                            facet_by    = c("pool", "none"),
                            free_scales = TRUE,
                            ncol        = NULL,
                            drop_na     = TRUE) {

  metric   <- match.arg(metric)
  facet_by <- match.arg(facet_by)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_eq_compare() needs the ggplot2 package.")
  library(ggplot2)

  # ---- accept a named list of eq_compare data frames, or a combined df ----
  if (is.data.frame(compare_list)) {
    if (!all(c("model", "name", metric) %in% names(compare_list)))
      stop("combined data frame needs columns: model, name, ", metric)
    df <- data.frame(model = compare_list$model,
                     name  = as.character(compare_list$name),
                     value = suppressWarnings(as.numeric(compare_list[[metric]])),
                     stringsAsFactors = FALSE)
  } else {
    if (is.null(names(compare_list)) || any(names(compare_list) == ""))
      stop("compare_list must be a NAMED list, e.g. list(Century = df1, MIMICS = df2).")
    parts <- lapply(names(compare_list), function(m) {
      d <- compare_list[[m]]
      if (!all(c("name", metric) %in% names(d)))
        stop("Model '", m, "' is missing column(s): ",
             paste(setdiff(c("name", metric), names(d)), collapse = ", "))
      data.frame(model = m, name = as.character(d$name),
                 value = suppressWarnings(as.numeric(d[[metric]])),
                 stringsAsFactors = FALSE)
    })
    df <- do.call(rbind, parts)
  }

  if (drop_na) {
    nd <- sum(is.na(df$value))
    if (nd) message(nd, " row(s) dropped where ", metric,
                    " is NA (e.g. pools absent from the baseline).")
    df <- df[!is.na(df$value), , drop = FALSE]
  }
  if (!nrow(df)) stop("Nothing to plot after filtering.")

  # order pools by how many models share them, then name
  n_models    <- tapply(df$model, df$name, function(x) length(unique(x)))
  ord         <- names(sort(n_models[unique(df$name)], decreasing = TRUE))
  df$name     <- factor(df$name, levels = ord)
  df$model    <- factor(df$model, levels = unique(df$model))

  ylab <- switch(metric,
                 percent_change = "Percent change vs baseline (%)",
                 difference     = "Difference (baseline - treatment)",
                 treatment      = "Treatment equilibrium pool (g C m-2)",
                 baseline       = "Baseline equilibrium pool (g C m-2)")

  p <- ggplot(df, aes(x = model, y = value, fill = model)) +
    geom_col(width = 0.75) +
    labs(x = NULL, y = ylab, fill = "Model",
         title = "Equilibrium comparison across models",
         subtitle = if (facet_by == "pool")
           "One panel per pool; free axes so different magnitudes stay readable" else NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 30, hjust = 1))

  if (metric %in% c("percent_change", "difference"))
    p <- p + geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40")

  if (requireNamespace("RColorBrewer", quietly = TRUE))
    p <- p + scale_fill_brewer(palette = "Dark2")

  if (facet_by == "pool")
    p <- p + facet_wrap(~ name,
                        scales = if (free_scales) "free_y" else "fixed",
                        ncol = ncol)

  p
}
