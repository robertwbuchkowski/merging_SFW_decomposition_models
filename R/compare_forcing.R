#------------------------------------------------------------------------#
# Compare results under seasonal vs constant forcing ----
#------------------------------------------------------------------------#
# Helper functions for Scripts/8_compare_forcing.R. Each reads the saved
# outputs of Scripts/1, 4, 5 and 7 from the two results folders
# (results_dir_by_forcing in R/sensitivity_settings.R: Results/ and
# Results/constant_forcing/) and returns a tidy table and a figure. A section
# whose files are missing in either folder returns NULL with a message, so the
# comparison runs with whatever has been produced so far.
#
# Colours: seasonal = blue, constant = orange (validated categorical pair);
# point shape repeats the forcing so the figures also read in grey scale.

forcing_cols   <- c(seasonal = "#2a78d6", constant = "#eb6834")
forcing_shapes <- c(seasonal = 16, constant = 17)
forcing_labs   <- c(seasonal = "Seasonal forcing (limit cycle)",
                    constant = "Constant mean forcing")

#------------------------------------------------------------------------#
# read_both() ----
#------------------------------------------------------------------------#
# Used in: every cmp_*() function below
# Read the same file from both results folders and stack them with a
# `forcing` column. NULL (with a message) if it is missing in either.
read_both <- function(file, dirs = results_dir_by_forcing) {
  paths <- file.path(dirs, file); names(paths) <- names(dirs)
  miss  <- paths[!file.exists(paths)]
  if (length(miss)) {
    message("  skipped (", file, " missing in: ", paste(dirname(miss), collapse = ", "), ")")
    return(NULL)
  }
  bind_rows(lapply(names(paths), function(f)
    mutate(read_csv(paths[[f]], show_col_types = FALSE, guess_max = 1e5), forcing = f)))
}

# shared look of the comparison figures
theme_forcing <- function(base_size = 10)
  theme_minimal(base_size = base_size) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

scale_forcing <- function()
  list(scale_colour_manual(values = forcing_cols, labels = forcing_labs, name = NULL),
       scale_shape_manual(values = forcing_shapes, labels = forcing_labs, name = NULL))

#------------------------------------------------------------------------#
# cmp_calibration() ----
#------------------------------------------------------------------------#
# Used in: Scripts/8_compare_forcing.R (section: Calibration)
# Fitted animal parameters (Scripts/1) under both forcings.
# Returns list(table, plot); table has one row per scenario x parameter with
# both values and the ratio seasonal / constant.
cmp_calibration <- function() {
  d <- read_both("fitted_animal_params.csv")
  if (is.null(d)) return(NULL)
  tab <- d %>% select(model, scenario, animal, param, role, forcing, value) %>%
    pivot_wider(names_from = forcing, values_from = value) %>%
    mutate(ratio_seasonal_to_constant = seasonal / constant)
  p <- ggplot(d, aes(value, pretty_scenario(scenario), colour = forcing, shape = forcing)) +
    geom_line(aes(group = scenario), colour = "grey75", linewidth = 0.6) +
    geom_point(size = 3) +
    facet_wrap(~param, scales = "free_x") +
    scale_x_log10() + scale_forcing() +
    labs(x = "Fitted value (log scale)", y = NULL,
         title = "Calibrated animal parameters",
         subtitle = "Each animal fitted to the same target biomass under each forcing") +
    theme_forcing()
  list(table = tab, plot = p)
}

#------------------------------------------------------------------------#
# cmp_headline() ----
#------------------------------------------------------------------------#
# Used in: Scripts/8_compare_forcing.R (section: Headline effect)
# Headline animal effect on total soil + root C (Scripts/5): the effect at the
# default parameters and the 95% interval of the accepted Latin-hypercube
# samples, total and direct, under both forcings.
cmp_headline <- function(subset_used = "accepted") {
  s <- read_both("animal_effect_uncertainty_summary.csv")
  if (is.null(s)) return(NULL)
  s <- s %>% filter(subset == subset_used)
  tab <- s %>% select(scenario, forcing, n, total_default, total_median, total_q025,
                      total_q975, p_total_positive, direct_default, direct_median,
                      direct_q025, direct_q975, pct_direct_default, pct_direct_median)
  long <- bind_rows(
    s %>% transmute(scenario, forcing, effect = "Total effect", default = total_default,
                    median = total_median, lo = total_q025, hi = total_q975),
    s %>% transmute(scenario, forcing, effect = "Direct effect (feeding only)",
                    default = direct_default, median = direct_median,
                    lo = direct_q025, hi = direct_q975)) %>%
    mutate(effect = factor(effect, levels = c("Total effect", "Direct effect (feeding only)")))
  dodge <- position_dodge(width = 0.5)
  p <- ggplot(long, aes(y = pretty_scenario(scenario), colour = forcing, shape = forcing)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y",
                  position = dodge) +
    geom_point(aes(x = default), size = 3, position = dodge) +
    facet_wrap(~effect, scales = "free_x") + scale_forcing() +
    labs(x = expression("Animal effect on total soil + root C (g C m"^-2*")"), y = NULL,
         title = "Headline animal effect",
         subtitle = paste0("Point = default parameters; bar = 95% interval of the ",
                           subset_used, " Latin-hypercube samples")) +
    theme_forcing()
  list(table = tab, plot = p)
}

#------------------------------------------------------------------------#
# cmp_headline_samples() ----
#------------------------------------------------------------------------#
# Used in: Scripts/8_compare_forcing.R (section: Headline effect)
# Distribution of the TOTAL effect over the accepted Latin-hypercube samples,
# both forcings side by side (violin + median), with the acceptance rate of
# the baseline realism filter in the facet label.
cmp_headline_samples <- function() {
  d  <- read_both("animal_effect_uncertainty_samples.csv")
  dg <- read_both("animal_effect_uncertainty_diagnostics.csv")
  if (is.null(d) || is.null(dg)) return(NULL)
  d <- d %>% filter(accepted %in% c(TRUE, "TRUE", 1))
  if (!nrow(d)) { message("  no accepted samples in either run"); return(NULL) }
  acc <- dg %>% select(scenario, forcing, n_accepted, acceptance_rate)
  p <- ggplot(d, aes(forcing, total, fill = forcing)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_violin(colour = NA, alpha = 0.55, scale = "width") +
    stat_summary(fun = median, geom = "point", shape = 21, size = 2.5,
                 fill = "white", colour = "grey20") +
    facet_wrap(~scenario, scales = "free_y", labeller = scenario_labeller) +
    scale_fill_manual(values = forcing_cols, labels = forcing_labs, name = NULL) +
    scale_x_discrete(labels = c(seasonal = "Seasonal", constant = "Constant")) +
    labs(x = NULL, y = expression("Total animal effect (g C m"^-2*")"),
         title = "Headline effect across accepted parameter samples",
         subtitle = "White point = median") +
    theme_forcing()
  list(table = acc, plot = p)
}

#------------------------------------------------------------------------#
# cmp_morris() ----
#------------------------------------------------------------------------#
# Used in: Scripts/8_compare_forcing.R (section: Morris)
# Morris screening (Scripts/4) under both forcings, for one range definition
# (tag "" = main, "_standardized"). mu* is scaled to the most influential
# parameter within each scenario x forcing (raw mu* is not comparable when the
# effect size differs). Per scenario: Spearman rank correlation of mu* over
# all parameters, overlap of the top-k sets, and the convergence status of
# each run (if the *_convergence.csv files exist).
# Returns list(table, params, plot).
cmp_morris <- function(tag = "", n_top = 10, top_k = 5) {
  d <- read_both(paste0("animal_effect_morris", tag, ".csv"))
  if (is.null(d)) return(NULL)
  d <- d %>% group_by(scenario, forcing) %>%
    mutate(s_max = max(mu_star, na.rm = TRUE), rel = mu_star / s_max,
           rel_lo = mu_star_lo / s_max, rel_hi = mu_star_hi / s_max,
           rank = rank(-mu_star, ties.method = "first", na.last = "keep")) %>% ungroup()
  wide <- d %>% select(scenario, parameter, parameter_label, forcing, mu_star, rank) %>%
    pivot_wider(names_from = forcing, values_from = c(mu_star, rank))
  stats <- wide %>% group_by(scenario) %>%
    summarise(n_parameters = n(),
              spearman_mu_star = suppressWarnings(cor(mu_star_seasonal, mu_star_constant,
                                                      method = "spearman", use = "complete.obs")),
              top_k = top_k,
              top_k_shared = sum(rank_seasonal <= top_k & rank_constant <= top_k, na.rm = TRUE),
              top_k_seasonal = paste(parameter[order(rank_seasonal)][seq_len(top_k)], collapse = ", "),
              top_k_constant = paste(parameter[order(rank_constant)][seq_len(top_k)], collapse = ", "),
              .groups = "drop")
  cv <- read_both(paste0("animal_effect_morris", tag, "_convergence.csv"))
  if (!is.null(cv)) {
    cv <- cv %>% filter(effects == "all") %>% select(scenario, forcing, status) %>%
      pivot_wider(names_from = forcing, values_from = status, names_prefix = "convergence_")
    stats <- left_join(stats, cv, by = "scenario")
  }
  keep <- d %>% filter(rank <= n_top) %>% distinct(scenario, parameter)
  pd <- d %>% semi_join(keep, by = c("scenario", "parameter")) %>%
    group_by(scenario, parameter) %>% mutate(rel_s = rel[forcing == "seasonal"][1]) %>% ungroup() %>%
    mutate(key = paste(parameter_label, scenario, sep = "___"))
  pd$key <- factor(pd$key, levels = unique(pd$key[order(pd$scenario, pd$rel_s)]))
  nudge <- c(seasonal = 0.18, constant = -0.18)
  p <- ggplot(pd, aes(y = key, colour = forcing, shape = forcing)) +
    lapply(names(nudge), function(f) list(
      geom_errorbar(data = filter(pd, forcing == f), aes(xmin = rel_lo, xmax = rel_hi),
                    width = 0.25, orientation = "y", position = position_nudge(y = nudge[[f]])),
      geom_point(data = filter(pd, forcing == f), aes(x = rel), size = 2,
                 position = position_nudge(y = nudge[[f]])))) +
    facet_wrap(~scenario, scales = "free_y", labeller = scenario_labeller) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) + scale_forcing() +
    labs(x = expression(mu*"* relative to the most influential parameter (95% bootstrap CI)"),
         y = NULL,
         title = paste0("Morris screening", if (nzchar(tag)) paste0(" (", sub("^_", "", tag), ")") else ""),
         subtitle = paste0("Top ", n_top, " parameters under either forcing, ordered by the seasonal result"),
         caption = paste0("Spearman rho (all parameters): ",
                          paste(sprintf("%s %.2f", pretty_scenario(stats$scenario),
                                        stats$spearman_mu_star), collapse = "; "))) +
    theme_forcing(9)
  list(table = stats, params = wide, plot = p)
}

#------------------------------------------------------------------------#
# cmp_identifiability() ----
#------------------------------------------------------------------------#
# Used in: Scripts/8_compare_forcing.R (section: Identifiability)
# Identifiability (Scripts/7) under both forcings: number of identifiable
# parameters per data design and parameter set, and the status of every
# animal parameter (identifiable / collinear / no information) side by side.
# Returns list(table, params, plot_counts, plot_animal).
cmp_identifiability <- function() {
  s <- read_both("identifiability_summary.csv")
  p <- read_both("identifiability_parameters.csv")
  if (is.null(s) || is.null(p)) return(NULL)
  designs <- unique(s$design)
  tab <- s %>% select(scenario, param_set, design, forcing, n_params, n_informative,
                      n_identifiable_greedy, identifiable_greedy) %>%
    pivot_wider(names_from = forcing,
                values_from = c(n_informative, n_identifiable_greedy, identifiable_greedy))
  cd <- s %>% mutate(design = factor(design, levels = rev(designs)),
                     param_set = factor(ifelse(param_set == "all", "All parameters",
                                               "Animal parameters only"),
                                        levels = c("All parameters", "Animal parameters only")))
  dodge <- position_dodge(width = 0.6)
  p1 <- ggplot(cd, aes(n_identifiable_greedy, design, colour = forcing, shape = forcing)) +
    geom_linerange(aes(xmin = 0, xmax = n_identifiable_greedy), position = dodge,
                   linewidth = 0.8, alpha = 0.6) +
    geom_point(size = 2.8, position = dodge) +
    facet_grid(param_set ~ scenario, scales = "free_x",
               labeller = labeller(scenario = scenario_labeller)) +
    scale_forcing() +
    labs(x = "Identifiable parameters (greedy subset)", y = "Artificial data used",
         title = "Parameter identifiability") +
    theme_forcing()
  st <- p %>% filter(param_set == "animal") %>%
    mutate(design = factor(design, levels = designs),
           forcing = factor(forcing, levels = c("seasonal", "constant")),
           status = factor(status, levels = c("identifiable", "collinear", "no information")),
           key = paste(pretty_param(parameter), scenario, sep = "___"))
  p2 <- ggplot(st, aes(forcing, key, fill = status)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    facet_grid(scenario ~ design, scales = "free_y", space = "free_y",
               labeller = labeller(scenario = scenario_labeller,
                                   design = label_wrap_gen(22))) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_x_discrete(labels = c(seasonal = "Seasonal", constant = "Constant")) +
    scale_fill_manual(values = c(identifiable = "#2166ac", collinear = "#f4a582",
                                 `no information` = "grey85"), name = NULL, drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Animal parameters: identifiability status by forcing") +
    theme_forcing(9) + theme(panel.grid = element_blank(),
                             strip.text.y = element_text(angle = 0))
  changed <- st %>% select(scenario, design, parameter, forcing, status) %>%
    pivot_wider(names_from = forcing, values_from = status) %>%
    mutate(changed = as.character(seasonal) != as.character(constant))
  list(table = tab, params = changed, plot_counts = p1, plot_animal = p2)
}

#------------------------------------------------------------------------#
# cmp_overview() ----
#------------------------------------------------------------------------#
# Used in: Scripts/8_compare_forcing.R (section: Overview)
# One row per scenario with the headline numbers of every section that ran.
cmp_overview <- function(res) {
  out <- NULL
  add <- function(x) out <<- if (is.null(out)) x else full_join(out, x, by = "scenario")
  if (!is.null(res$calibration))
    add(res$calibration$table %>% filter(grepl("^biomass", role)) %>%
          transmute(scenario, feeding_param = param, feeding_seasonal = seasonal,
                    feeding_constant = constant, feeding_ratio = ratio_seasonal_to_constant))
  if (!is.null(res$headline))
    add(res$headline$table %>% select(scenario, forcing, total_default, total_q025, total_q975) %>%
          pivot_wider(names_from = forcing, values_from = c(total_default, total_q025, total_q975)))
  if (!is.null(res$morris))
    add(res$morris$table %>% select(scenario, spearman_mu_star, top_k_shared,
                                    any_of(c("convergence_seasonal", "convergence_constant"))))
  if (!is.null(res$ident))
    add(res$ident$table %>% filter(param_set == "animal", design == "All data") %>%
          transmute(scenario, animal_ident_all_data_seasonal = n_identifiable_greedy_seasonal,
                    animal_ident_all_data_constant = n_identifiable_greedy_constant))
  out
}
