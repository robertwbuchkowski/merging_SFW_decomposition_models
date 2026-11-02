#------------------------------------------------------------------------#
# Carbon flux network at equilibrium ----
#------------------------------------------------------------------------#
# Draws the soil food-web / Millennial carbon network as a diagram: one box
# per pool and one arrow per flux, with the flux value written on the arrow.
# Three panels (choose with `panels`):
#   "baseline"  fluxes at the no-animal equilibrium
#   "animal"    fluxes at the equilibrium with the scenario's animals
#   "delta"     with - without animals (blue = flux increased, red = decreased)
# Arrow width scales with sqrt(flux); every flux value is written on its
# arrow, and labels are slid along their arrows to avoid boxes, other labels
# and other arrows. All pools are always drawn (empty ones faded).
# Inputs (NPP) and losses (respiration, leaching) are drawn as short stub
# arrows on the pool they enter or leave, so the diagram has no long arrows
# running to a shared "CO2" box.
#
# All fluxes are at the CONSTANT-FORCING equilibrium (the same steady state
# the animals are fitted at, see spinup_equilibrium()), in g C m-2 yr-1.
# Stocks are in g C m-2.
#
# Main functions
#   millennial_fluxes(state, parms)       -> edge table for one state
#   flux_network_scenario(model, scenario, scen, fitted)
#                                         -> both equilibria + flux table
#   plot_flux_network(net)                -> ggplot (one or more panels)
#   save_flux_network(net, file)          -> PNG sized for the panels (+ .csv)
#   plot_flux_network_deltas(nets) / save_flux_network_deltas(nets, file)
#                                         -> change panels of several scenarios,
#                                            one panel per scenario
#   check_flux_layout()                   -> arrows passing under a box (should
#                                            be empty after editing the layout)
#
# Example (from the project root, after Scripts/1_fit_all_animals.R):
#   source("R/flux_network.R")             # plus the usual R/ helpers
#   scen <- read_scenarios("Data/scenarios.xlsx")
#   fit  <- load_fitted_params("Results/fitted_animal_params.csv")
#   net  <- flux_network_scenario("millennial", "RootHerbivore", scen, fit)
#   save_flux_network(net, "Results/figures/flux_network_RootHerbivore.png")
#   plot_flux_network(net, panels = "delta")          # change only

#------------------------------------------------------------------------#
# Flux definitions ----
#------------------------------------------------------------------------#
# One row per flux. `expr` is evaluated in the model's local environment, so
# it can use any variable defined inside millennial_model_wplant(). from/to
# "IN" = external input, "CO2" = respiration, "LEACH" = leaching loss.
# Every term of every derivative in the model appears exactly once here; the
# mass-balance check in millennial_fluxes() confirms this.
.flux_defs <- data.frame(
  from = c(
    "IN", "IN", "IN", "IN",
    "C_root_herb", "C_root_herb", "C_root_herb",
    "C_root_tree", "C_root_tree",
    "Litter", "Litter", "Litter", "Litter",
    "CWD", "CWD",
    "P", "P", "P", "P",
    "A", "A", "A",
    "M", "M", "M",
    "L", "L", "L", "L",
    "B", "B", "B", "B",
    "Earthworm", "Earthworm", "Earthworm", "Earthworm",
    "Detritivore", "Detritivore",
    "RootHerb", "RootHerb"),
  to = c(
    "Litter", "CWD", "C_root_herb", "C_root_tree",
    "P", "L", "RootHerb",
    "P", "L",
    "P", "L", "Earthworm", "Detritivore",
    "P", "L",
    "A", "L", "Earthworm", "Detritivore",
    "P", "M", "Earthworm",
    "A", "L", "Earthworm",
    "M", "B", "LEACH", "Earthworm",
    "M", "L", "CO2", "Detritivore",
    "L", "A", "P", "CO2",
    "P", "CO2",
    "P", "CO2"),
  expr = c(
    "leaf_litter_input", "cwd_input", "root_growth_herb", "root_growth_tree",
    "root_mortality_herb", "exudates_herb", "Fed_rootherb_herb",
    "root_mortality_tree", "exudates_tree",
    "fragmentation_litter", "F_Litter_L", "Fed_earthworm_litter", "Fed_det_lit",
    "fragmentation_CWD", "F_CWD_L",
    "F_pa", "F_pl", "Fed_earthworm_P", "Fed_det_P",
    "p_a * F_a", "(1 - p_a) * F_a", "Fed_earthworm_A",
    "F_ma", "F_ld", "Fed_earthworm_M",
    "F_lm", "F_lb", "F_l", "Fed_earthworm_L",
    "p_b * F_bm", "(1 - p_b) * F_bm", "F_mr", "Fed_det_mic",
    "Waste_earthworm_L", "Waste_earthworm_A", "Carcass_earthworm_P", "Respiration_earthworm",
    "Waste_det_P + Carcass_det_P", "Respiration_detritivore",
    "Waste_rootherb_P + Carcass_rootherb_P", "Respiration_rootherb"),
  stringsAsFactors = FALSE)

#------------------------------------------------------------------------#
# Node layout and labels ----
#------------------------------------------------------------------------#
# Fixed positions chosen so that no arrow passes under a box it does not
# start or end at (check_flux_layout() tests this). POM sits in the middle
# because almost every pool exchanges C with it. Rows from top to bottom:
# plant inputs (roots, litter, CWD); fauna, POM and LMWC; earthworms and
# microbes; aggregates and MAOM.
.flux_nodes <- data.frame(
  pool  = c("C_root_herb", "C_root_tree", "Litter", "CWD",
            "RootHerb", "Detritivore", "P", "L",
            "Earthworm", "B", "A", "M"),
  label = c("Herb roots", "Tree roots", "Litter", "CWD",
            "Root herbivores", "Detritivores", "POM", "LMWC",
            "Earthworms", "Microbes", "Aggregates", "MAOM"),
  group = c("plant", "plant", "detritus", "detritus",
            "animal", "animal", "soil", "soil",
            "animal", "microbe", "soil", "soil"),
  x     = c(-1.0,  3.5,  8.0, 12.5,
            -4.2,  0.0,  5.5, 12.0,
             7.7, 15.5,  1.5,  9.8),
  y     = c(10.0, 10.0, 10.0, 10.0,
             7.8,  6.1,  6.5,  6.5,
             4.0,  3.0,  1.0,  0.0),
  stringsAsFactors = FALSE)

# stub direction for inputs (pointing into the box) and losses (out of it)
.flux_stubs <- data.frame(
  pool = c("Litter", "CWD", "C_root_herb", "C_root_tree",
           "B", "Earthworm", "Detritivore", "RootHerb", "L"),
  kind = c("IN", "IN", "IN", "IN",
           "CO2", "CO2", "CO2", "CO2", "LEACH"),
  dx   = c(0, 0, 0, 0,  1, 0,  0, -1,  1),
  dy   = c(1, 1, 1, 1,  0, -1, -1,  0,  0),
  stringsAsFactors = FALSE)

# the node table with its vertical positions multiplied by y_scale (< 1
# packs the rows closer, e.g. for the page-width multi-scenario figure). In
# a packed layout the earthworm box moves down a little more, to keep the
# microbe -> detritivore arrow clear of both it and POM.
.nodes_at <- function(y_scale = 1) {
  n <- .flux_nodes; n$y <- n$y * y_scale
  if (y_scale < 1) n$y[n$pool == "Earthworm"] <- n$y[n$pool == "Earthworm"] - 0.3
  n
}

.flux_node_fill <- c(plant = "#d9ead3", detritus = "#f4e3c1", animal = "#fbd5b5",
                     soil = "#e6dccf", microbe = "#d6e4f0")

#------------------------------------------------------------------------#
# millennial_fluxes() ----
#------------------------------------------------------------------------#
# Every flux of millennial_model_wplant() for one state, as an edge table
# (from, to, flux) in g C m-2 yr-1. The model is evaluated once with its
# last line replaced so it returns all of its local variables; the fluxes are
# then read from those. forcing defaults to the constant mean climate (the
# same forcing spinup_equilibrium() uses), with daily weights of 1/365, so
# daily fluxes are multiplied by 365.
# Stops if the fluxes in the table do not reproduce every derivative (i.e.
# the model gained a flux that .flux_defs does not list).
millennial_fluxes <- function(state, parms, time = 0,
                              forcing = make_climate_forcing_equilibrium(parms),
                              per_year = 365, tol = 1e-8) {
  fl <- millennial_model_wplant
  body(fl)[[length(body(fl))]] <- quote(return(environment()))
  parms$climate_forcing <- forcing
  env <- fl(time, state, parms)

  flux <- vapply(.flux_defs$expr, function(e) {
    v <- eval(parse(text = e), envir = env)
    as.numeric(unname(v))[1]
  }, numeric(1))
  edges <- data.frame(.flux_defs[, c("from", "to")], flux = per_year * unname(flux),
                      stringsAsFactors = FALSE)

  # mass balance: inflow - outflow of each pool must equal the derivative
  dvec <- per_year * unlist(lapply(.POOLS_MILLENNIAL, function(p)
    as.numeric(unname(get(paste0("d", p), envir = env)))[1]))
  names(dvec) <- .POOLS_MILLENNIAL
  net <- vapply(.POOLS_MILLENNIAL, function(p)
    sum(edges$flux[edges$to == p]) - sum(edges$flux[edges$from == p]), numeric(1))
  err <- abs(net - dvec)
  if (any(err > tol * pmax(1, abs(dvec), sum(abs(edges$flux)))))
    stop("flux table does not close the budget of: ",
         paste(names(err)[err > tol * pmax(1, abs(dvec))], collapse = ", "),
         " -- update .flux_defs in R/flux_network.R to match the model.")

  stocks <- setNames(rep(0, length(.POOLS_MILLENNIAL)), .POOLS_MILLENNIAL)
  shared <- intersect(names(state), .POOLS_MILLENNIAL)
  stocks[shared] <- state[shared]
  list(edges = edges, stocks = stocks, ddt = dvec)
}

#------------------------------------------------------------------------#
# flux_network_scenario() ----
#------------------------------------------------------------------------#
# Equilibria with and without the scenario's animals (same steps as
# Scripts/2_spinup_dynamic.R: baseline equilibrium, apply fitted feeding
# rates, treatment equilibrium warm-started from the baseline), then the flux
# table for both and their difference.
# Returns list(edges, stocks, model, scenario) where
#   edges : from, to, baseline, animal, delta   (g C m-2 yr-1)
#   stocks: pool, baseline, animal, delta        (g C m-2)
# `pair` can be supplied to skip the setup (a setup_scenario_pair() list,
# with equilibria already in init_state_spin, or not).
flux_network_scenario <- function(model, scenario, scen, fitted = NULL,
                                  pair = NULL, verbose = FALSE) {
  if (model != "millennial") stop("flux_network_scenario() supports the millennial model only")
  if (is.null(pair)) pair <- setup_scenario_pair(model, scen, scenario)
  if (is.null(pair$baseline$init_state_spin))
    pair$baseline <- spinup_equilibrium(pair$baseline, verbose = verbose)
  if (!is.null(fitted))
    pair$treatment <- apply_fitted_params(pair$treatment, fitted, model, scenario,
                                          verbose = verbose)
  if (is.null(pair$treatment$init_state_spin))
    pair$treatment <- spinup_equilibrium(pair$treatment,
                                         warm_start = pair$baseline$init_state_spin,
                                         verbose = verbose)
  if (!isTRUE(pair$baseline$spin_info$converged) ||
      !isTRUE(pair$treatment$spin_info$converged))
    warning("equilibrium not converged for ", scenario, " (fluxes are not balanced)")

  fb <- millennial_fluxes(pair$baseline$init_state_spin,  pair$baseline$parms)
  ft <- millennial_fluxes(pair$treatment$init_state_spin, pair$treatment$parms)

  edges <- data.frame(fb$edges[, c("from", "to")],
                      baseline = fb$edges$flux, animal = ft$edges$flux)
  edges$delta <- edges$animal - edges$baseline
  stocks <- data.frame(pool = .POOLS_MILLENNIAL,
                       baseline = unname(fb$stocks), animal = unname(ft$stocks))
  stocks$delta <- stocks$animal - stocks$baseline
  list(edges = edges, stocks = stocks, model = model, scenario = scenario)
}

#------------------------------------------------------------------------#
# Drawing helpers ----
#------------------------------------------------------------------------#
# number format: 3 significant figures, thousands separator, no sci notation
.fmt_flux <- function(x, signed = FALSE) {
  s <- formatC(signif(x, 3), format = "fg", digits = 3, big.mark = ",",
               flag = if (signed) "+" else "")
  s <- sub("\\.$", "", trimws(s))
  s[abs(x) < 5e-4] <- "0"
  s
}

# quadratic Bezier from (x0,y0) to (x1,y1) bowed sideways by `curv` x length
# (positive = to the left of the direction of travel), trimmed where it
# enters the start and end boxes (half sizes hw, hh) plus a gap.
.edge_path <- function(x0, y0, x1, y1, curv, hw, hh, gap = 0.12, n = 160) {
  dx <- x1 - x0; dy <- y1 - y0
  cx <- (x0 + x1) / 2 - curv * dy
  cy <- (y0 + y1) / 2 + curv * dx
  t  <- seq(0, 1, length.out = n)
  px <- (1 - t)^2 * x0 + 2 * (1 - t) * t * cx + t^2 * x1
  py <- (1 - t)^2 * y0 + 2 * (1 - t) * t * cy + t^2 * y1
  inside <- function(xc, yc) abs(px - xc) <= hw + gap & abs(py - yc) <= hh + gap
  keep <- which(!inside(x0, y0) & !inside(x1, y1))
  if (length(keep) < 2) keep <- c(1, n)
  idx <- min(keep):max(keep)
  list(x = px[idx], y = py[idx])
}

# Curvature of every pool-to-pool arrow. Two-way pairs bow apart so both
# arrows and both labels are visible; a few single arrows are bent to keep
# them clear of boxes they would otherwise run under.
.edge_curvature <- function(from, to, y_scale = 1) {
  key  <- paste(from, to)
  curv <- ifelse(paste(to, from) %in% key, 0.20, 0)
  bend <- c("B Detritivore" = -0.07, "C_root_herb L" = -0.06,
            "Litter Earthworm" = 0, "RootHerb P" = 0.05, "B M" = -0.12)
  hit <- key %in% names(bend)
  curv[hit] <- bend[key[hit]]
  if (y_scale < 1) curv[key == "B Detritivore"] <- -0.03   # packed layout
  curv
}

# all pool-to-pool arrows of the model with their drawn paths
.network_paths <- function(hw, hh, y_scale = 1) {
  nodes <- .nodes_at(y_scale)
  e <- unique(.flux_defs[.flux_defs$from %in% nodes$pool &
                         .flux_defs$to %in% nodes$pool, c("from", "to")])
  e$curv <- .edge_curvature(e$from, e$to, y_scale)
  a <- match(e$from, nodes$pool); b <- match(e$to, nodes$pool)
  e$path <- lapply(seq_len(nrow(e)), function(i)
    .edge_path(nodes$x[a[i]], nodes$y[a[i]], nodes$x[b[i]], nodes$y[b[i]],
               e$curv[i], hw, hh))
  e
}

#------------------------------------------------------------------------#
# check_flux_layout() ----
#------------------------------------------------------------------------#
# Checks the fixed layout: lists every arrow that passes under a box other
# than its own two ends (with `margin` of clearance). An empty result means
# the layout is clean. Run after editing .flux_nodes or .edge_curvature().
check_flux_layout <- function(hw = 1.45, hh = 0.5, margin = 0.1, y_scale = 1) {
  nodes <- .nodes_at(y_scale)
  e <- .network_paths(hw, hh, y_scale)
  out <- list()
  for (i in seq_len(nrow(e))) {
    p <- e$path[[i]]
    for (k in seq_len(nrow(nodes))) {
      if (nodes$pool[k] %in% c(e$from[i], e$to[i])) next
      if (any(abs(p$x - nodes$x[k]) < hw + margin & abs(p$y - nodes$y[k]) < hh + margin))
        out[[length(out) + 1]] <- data.frame(from = e$from[i], to = e$to[i],
                                             passes_under = nodes$pool[k])
    }
  }
  if (length(out)) do.call(rbind, out)
  else data.frame(from = character(0), to = character(0), passes_under = character(0))
}

# rectangle overlap (rectangles as c(xmin, xmax, ymin, ymax))
.rect_overlap <- function(r, R) {
  if (!length(R)) return(0)
  R <- matrix(R, ncol = 4, byrow = FALSE)
  sum(r[1] < R[, 2] & r[2] > R[, 1] & r[3] < R[, 4] & r[4] > R[, 3])
}

# Greedy label placement: each label slides along its own arrow (12-88% of
# its length) to the spot with the fewest clashes with boxes, other labels
# and other arrows. Arrows with the least room are placed first.
.place_labels <- function(lab, paths, box_rects, fixed_rects, char_w, line_h) {
  placed <- fixed_rects
  lab$x <- NA_real_; lab$y <- NA_real_
  ord <- order(vapply(paths, function(p) length(p$x), 1))
  other_pts <- function(i) {
    o <- paths[-i]
    list(x = unlist(lapply(o, `[[`, "x")), y = unlist(lapply(o, `[[`, "y")),
         id = rep(seq_along(o), vapply(o, function(p) length(p$x), 1)))
  }
  for (i in ord) {
    p  <- paths[[i]]; n <- length(p$x)
    w  <- nchar(lab$label[i]) * char_w + 0.18
    h  <- line_h
    op <- other_pts(i)
    best <- Inf
    for (f in seq(0.12, 0.88, by = 0.02)) {
      j <- max(1, min(n, round(f * n)))
      r <- c(p$x[j] - w / 2, p$x[j] + w / 2, p$y[j] - h / 2, p$y[j] + h / 2)
      hits_box <- .rect_overlap(r, box_rects)
      hits_lab <- .rect_overlap(r, placed)
      inb <- op$x > r[1] & op$x < r[2] & op$y > r[3] & op$y < r[4]
      n_cross <- length(unique(op$id[inb]))
      score <- 100 * hits_box + 100 * hits_lab + 4 * n_cross + 3 * abs(f - 0.5)
      if (score < best) { best <- score; lab$x[i] <- p$x[j]; lab$y[i] <- p$y[j]; rb <- r }
    }
    placed <- rbind(placed, rb)
  }
  lab
}

#------------------------------------------------------------------------#
# plot_flux_network() ----
#------------------------------------------------------------------------#
# net     : output of flux_network_scenario()
# panels  : any of "baseline", "animal", "delta"
# ncol    : panels per row (default 1: panels stacked, the most readable)
# width   : width of ONE panel in inches as it will be saved (used to size
#           the label boxes so labels do not overlap; save_flux_network()
#           passes it through)
# min_width / max_width : arrow line widths (mm) for the smallest / largest
#           flux; width scales with sqrt(flux) so small fluxes stay visible.
#           "No animals" and "With animals" share one scale; "Change" has its own.
# label_size, node_size : text sizes (mm) of flux values and box text
# Arrows with zero flux in a panel (feeding by absent animals) are not drawn;
# pools that are empty in a panel are drawn faded with a dashed border.
plot_flux_network <- function(net, panels = c("baseline", "animal", "delta"),
                              ncol = 1, width = 13,
                              min_width = 0.35, max_width = 4.5,
                              label_size = 3.4, node_size = 3.6,
                              hw = 1.45, hh = 0.5, stub_len = 1.1,
                              title = NULL) {
  L <- .flux_network_layers(net, panels, width, min_width, max_width,
                            label_size, node_size, hw, hh, stub_len)
  if (is.null(title))
    title <- sprintf("Carbon fluxes at equilibrium: %s",
                     if (exists("pretty_scenario")) pretty_scenario(net$scenario) else net$scenario)
  subtitle <- paste0("Arrows: fluxes (g C m⁻² yr⁻¹). Boxes: stocks ",
                     "(g C m⁻²). Equilibrium under constant mean climate.",
                     if ("delta" %in% panels)
                       "\nChange panel: blue = larger with animals, red = smaller."
                     else "")
  .flux_network_ggplot(L, ncol, hw, hh, label_size, title, subtitle)
}

# .flux_network_layers(): the drawing data (arrow paths, value labels, boxes,
# stub texts) for the requested panels of one scenario. panel_names renames
# panels, e.g. c(delta = "Earthworm") for the multi-scenario figure.
.flux_network_layers <- function(net, panels, width, min_width, max_width,
                                 label_size, node_size, hw, hh, stub_len,
                                 panel_names = NULL, drop_zero = FALSE,
                                 drop_inputs = FALSE, y_scale = 1,
                                 uniform_text = FALSE) {
  panels <- intersect(c("baseline", "animal", "delta"), panels)
  panel_lab <- c(baseline = "Without animals", animal = "With animals",
                 delta = "Change (with − without animals)")
  if (!is.null(panel_names)) panel_lab[names(panel_names)] <- panel_names
  nodes <- .nodes_at(y_scale)
  lim <- .flux_limits(hw, hh, stub_len, drop_inputs, y_scale)
  units_per_mm <- diff(lim$x) / (width * 25.4 * 0.94)
  char_w <- 0.66 * label_size * units_per_mm      # bold digits ~0.6 em + margin
  line_h <- 1.55 * label_size * units_per_mm

  geo <- .network_paths(hw, hh, y_scale)
  key_geo <- paste(geo$from, geo$to)
  is_pool <- net$edges$from %in% nodes$pool & net$edges$to %in% nodes$pool
  pool_edges <- net$edges[is_pool, ]
  stub_edges <- net$edges[!is_pool, ]
  pool_edges$g <- match(paste(pool_edges$from, pool_edges$to), key_geo)

  ref_max <- max(abs(unlist(net$edges[, c("baseline", "animal")])), 1e-12)
  del_max <- max(abs(net$edges$delta), 1e-12)
  width_of <- function(v, pn)
    min_width + (max_width - min_width) * sqrt(abs(v) / if (pn == "delta") del_max else ref_max)
  col_of <- function(v, pn) {
    if (pn != "delta") return(rep("#3d3d3d", length(v)))
    ifelse(abs(v) < 5e-4, "#a0a0a0", ifelse(v > 0, "#2166ac", "#b2182b"))
  }
  box_rects <- cbind(nodes$x - hw, nodes$x + hw, nodes$y - hh, nodes$y + hh)

  paths <- list(); labs <- list(); boxes <- list(); stub_txt <- list()
  for (pn in panels) {
    plab <- panel_lab[[pn]]
    signed <- pn == "delta"

    # stubs first (their text is fixed; labels on arrows avoid it)
    fixed <- matrix(numeric(0), ncol = 4)
    for (j in seq_len(nrow(.flux_stubs))) {
      s <- .flux_stubs[j, ]
      rows <- if (s$kind == "IN") stub_edges$from == "IN" & stub_edges$to == s$pool
              else stub_edges$from == s$pool & stub_edges$to == s$kind
      v   <- sum(stub_edges[[pn]][rows])
      ref <- sum(pmax(abs(stub_edges$baseline[rows]), abs(stub_edges$animal[rows])))
      if (ref == 0 || (pn != "delta" && v == 0)) next
      if (drop_inputs && s$kind == "IN") next                   # NPP stubs left out
      if (drop_zero && abs(v) < 5e-4) next                      # no change: not drawn
      k  <- match(s$pool, nodes$pool)
      ex <- nodes$x[k] + s$dx * (hw + 0.12); ey <- nodes$y[k] + s$dy * (hh + 0.12)
      fx <- ex + s$dx * stub_len;            fy <- ey + s$dy * stub_len
      xs <- if (s$kind == "IN") c(fx, ex) else c(ex, fx)
      ys <- if (s$kind == "IN") c(fy, ey) else c(ey, fy)
      id <- paste(pn, "stub", j)
      paths[[id]] <- data.frame(panel = plab, id = id, x = xs, y = ys,
                                lw = width_of(v, pn), col = col_of(v, pn))
      word <- c(IN = "NPP", CO2 = "CO₂", LEACH = "Leached")[[s$kind]]
      txt  <- paste(word, .fmt_flux(v, signed))
      tx <- fx + s$dx * 0.15; ty <- fy + s$dy * 0.15
      hj <- if (s$dx > 0) 0 else if (s$dx < 0) 1 else 0.5
      vj <- if (s$dy > 0) 0 else if (s$dy < 0) 1 else 0.5
      stub_txt[[id]] <- data.frame(panel = plab, x = tx, y = ty, label = txt,
                                   hjust = hj, vjust = vj, col = col_of(v, pn))
      tw <- nchar(txt) * char_w
      fixed <- rbind(fixed, c(tx - hj * tw, tx + (1 - hj) * tw,
                              ty - vj * line_h, ty + (1 - vj) * line_h))
    }

    # arrows between pools
    v_all <- pool_edges[[pn]]
    ref   <- pmax(abs(pool_edges$baseline), abs(pool_edges$animal))
    draw  <- which(if (pn == "delta") ref > 0 & (!drop_zero | abs(v_all) >= 5e-4)
                   else abs(v_all) > 0)
    plist <- lapply(draw, function(i) geo$path[[pool_edges$g[i]]])
    lab <- data.frame(label = .fmt_flux(v_all[draw], signed), col = col_of(v_all[draw], pn))
    lab <- .place_labels(lab, plist, box_rects, fixed, char_w, line_h)
    for (m in seq_along(draw)) {
      i <- draw[m]; id <- paste(pn, i)
      paths[[id]] <- data.frame(panel = plab, id = id, x = plist[[m]]$x, y = plist[[m]]$y,
                                lw = width_of(v_all[i], pn), col = col_of(v_all[i], pn))
    }
    labs[[pn]] <- data.frame(panel = plab, lab)

    # boxes: pool name + stock (or change in stock)
    st <- net$stocks[match(nodes$pool, net$stocks$pool), ]
    active <- if (pn == "delta") st$baseline > 0 | st$animal > 0 else st[[pn]] > 0
    stock_txt <- if (pn == "delta") {
      pct <- ifelse(st$baseline > 0, 100 * st$delta / st$baseline, NA)
      pct_txt <- ifelse(abs(pct) >= 0.1, sprintf("%+.1f%%", pct),
                        ifelse(abs(pct) < 0.001, "<0.001%", sprintf("%+.1g%%", pct)))
      paste0(.fmt_flux(st$delta, TRUE),
             ifelse(is.na(pct), " (new)",
                    ifelse(st$delta == 0, "", paste0(" (", pct_txt, ")"))))
    } else .fmt_flux(st[[pn]])
    text <- ifelse(active, paste0(nodes$label, "\n", stock_txt),
                   paste0(nodes$label, "\n(absent)"))
    # shrink the text of any box whose longest line would not fit
    longest <- vapply(strsplit(text, "\n"), function(z) max(nchar(z)), 1)
    fits    <- (2 * hw - 0.2) / (0.62 * node_size * units_per_mm)
    boxes[[pn]] <- data.frame(panel = plab, nodes, active = active, text = text,
                              tsize = if (uniform_text) node_size * min(1, fits / longest)
                                      else node_size * pmin(1, fits / longest))
  }
  paths    <- do.call(rbind, paths)
  labs     <- do.call(rbind, labs)
  boxes    <- do.call(rbind, boxes)
  stub_txt <- do.call(rbind, stub_txt)
  lev <- unname(panel_lab[panels])
  paths$panel <- factor(paths$panel, lev); labs$panel <- factor(labs$panel, lev)
  boxes$panel <- factor(boxes$panel, lev); stub_txt$panel <- factor(stub_txt$panel, lev)
  boxes$fill  <- ifelse(boxes$active, .flux_node_fill[boxes$group], "#f4f4f4")
  boxes$lty   <- ifelse(boxes$active, "solid", "22")
  boxes$tcol  <- ifelse(boxes$active, "black", "grey60")

  list(paths = paths, labs = labs, boxes = boxes, stub_txt = stub_txt, lim = lim)
}

# .flux_network_ggplot(): draw layers from .flux_network_layers() (one facet
# per panel).
.flux_network_ggplot <- function(L, ncol, hw, hh, label_size, title, subtitle, scale = 1,
                                 panel_tag = FALSE) {
  paths <- L$paths; labs <- L$labs; boxes <- L$boxes; stub_txt <- L$stub_txt; lim <- L$lim

  ggplot2::ggplot() +
    ggplot2::geom_path(data = paths,
      ggplot2::aes(x, y, group = id, linewidth = lw, colour = col),
      lineend = "butt", linejoin = "mitre",
      arrow = grid::arrow(length = grid::unit(0.22 * max(scale, 0.6), "cm"), type = "closed")) +
    ggplot2::geom_rect(data = boxes,
      ggplot2::aes(xmin = x - hw, xmax = x + hw, ymin = y - hh, ymax = y + hh,
                   fill = fill, linetype = lty),
      colour = "grey35", linewidth = 0.35 * max(scale, 0.6)) +
    ggplot2::geom_text(data = boxes,
      ggplot2::aes(x, y, label = text, colour = tcol, size = tsize),
      lineheight = 0.9) +
    ggplot2::scale_size_identity() +
    ggplot2::geom_label(data = labs,
      ggplot2::aes(x, y, label = label, colour = col),
      size = label_size, fill = "white", label.size = 0,
      label.padding = grid::unit(0.06, "lines"), fontface = "bold") +
    ggplot2::geom_label(data = stub_txt,
      ggplot2::aes(x, y, label = label, hjust = hjust, vjust = vjust, colour = col),
      size = label_size, fontface = "bold", fill = "white", label.size = 0,
      label.padding = grid::unit(0.04, "lines")) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_linetype_identity() +
    (if (panel_tag)       # scenario name inside the panel, top left (no strip)
      ggplot2::geom_text(data = data.frame(panel = factor(levels(boxes$panel), levels(boxes$panel))),
        ggplot2::aes(x = lim$x[1], y = lim$y[2], label = panel), hjust = 0, vjust = 1,
        fontface = "bold", size = 4.2 * scale / 1.35 * 1.6, colour = "black")) +
    ggplot2::facet_wrap(~ panel, ncol = ncol) +
    ggplot2::coord_fixed(xlim = lim$x, ylim = lim$y, expand = FALSE, clip = "off") +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void(base_size = 12 * scale) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   strip.text = ggplot2::element_text(face = "bold", size = 14 * scale,
                                                      hjust = if (scale < 1) 0 else 0.5,
                                                      margin = ggplot2::margin(8 * scale, 0, 4 * scale, 0)),
                   plot.title = ggplot2::element_text(face = "bold", size = 15 * scale),
                   plot.subtitle = ggplot2::element_text(size = 11 * scale,
                                                         margin = ggplot2::margin(0, 0, 6 * scale, 0)),
                   plot.margin = ggplot2::margin(12 * scale, 12 * scale, 12 * scale, 12 * scale),
                   panel.spacing = grid::unit(if (panel_tag) 0.6 else 1.5 * scale, "lines"),
                   strip.text.x = if (panel_tag) ggplot2::element_blank()
                                  else ggplot2::element_text(face = "bold", size = 14 * scale,
                                         margin = ggplot2::margin(8 * scale, 0, 4 * scale, 0)))
}

# plotting window that holds every box, stub and stub label
.flux_limits <- function(hw = 1.45, hh = 0.5, stub_len = 1.1, drop_inputs = FALSE,
                         y_scale = 1) {
  n <- .nodes_at(y_scale)
  top <- if (drop_inputs) max(n$y) + hh + 0.35 else max(n$y) + hh + stub_len + 0.8
  list(x = c(min(n$x) - hw - stub_len - 2.6, max(n$x) + hw + stub_len + 3.0),
       y = c(min(n$y) - hh - 0.5, top))
}
#------------------------------------------------------------------------#
# save_flux_network() ----
#------------------------------------------------------------------------#
# Writes plot_flux_network() to a PNG. `width` is the width of one panel in
# inches; the height follows from the layout's aspect ratio. Also writes the
# flux table next to it (same name, .csv) unless csv = FALSE.
save_flux_network <- function(net, file, panels = c("baseline", "animal", "delta"),
                              ncol = 1, width = 13, dpi = 200, csv = TRUE,
                              hw = 1.45, hh = 0.5, stub_len = 1.1, ...) {
  panels <- match.arg(panels, several.ok = TRUE)
  p <- plot_flux_network(net, panels = panels, ncol = ncol, width = width,
                         hw = hw, hh = hh, stub_len = stub_len, ...)
  lim  <- .flux_limits(hw, hh, stub_len)
  ph   <- width * 0.94 * diff(lim$y) / diff(lim$x)        # one panel's height
  nrow <- ceiling(length(panels) / ncol)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file, p, width = width * ncol, height = 1.1 + (ph + 0.55) * nrow,
                  dpi = dpi, bg = "white", limitsize = FALSE)
  if (csv) {
    out <- data.frame(scenario = net$scenario, net$edges)
    utils::write.csv(out, sub("\\.png$", ".csv", file), row.names = FALSE)
  }
  invisible(p)
}

#------------------------------------------------------------------------#
# plot_flux_network_deltas() / save_flux_network_deltas() ----
#------------------------------------------------------------------------#
# One figure with the CHANGE panel (with - without animals) of several
# scenarios, one panel per scenario labelled with its name, STACKED
# VERTICALLY and sized for a portrait page (default 6.5 in wide, i.e. a
# Letter/A4 page with 1-in margins), with text of ~6 pt at that size.
#   nets  named or unnamed list of flux_network_scenario() outputs
# To keep four panels legible at page width the figure leaves out what does
# not change: arrows and stubs whose change rounds to 0, and the NPP input
# stubs (NPP is a fixed input, identical with and without animals). All pools
# are still drawn. Arrow widths are scaled within each panel (the changes
# differ by orders of magnitude between scenarios), so compare widths within
# a panel and read the values across panels.
#   width       figure width in inches (as placed in the document)
#   label_size, node_size  text sizes in mm (2.1 mm ~ 6 pt)
#   ncol        panels per row (1 = vertical stack)
plot_flux_network_deltas <- function(nets, ncol = 1, width = 6.5,
                                     min_width = 0.25, max_width = 2.6,
                                     label_size = 2.1, node_size = 2.2,
                                     hw = 1.45, hh = 0.5, stub_len = 1.1, y_scale = 0.8,
                                     title = NULL,
                                     subtitle = paste0("Change with minus without animals. Arrows: g C m\u207b\u00b2 yr\u207b\u00b9; ",
                                                       "boxes: g C m\u207b\u00b2 (% of no-animal stock).\n",
                                                       "Blue = larger with animals, red = smaller; ",
                                                       "arrow widths scaled within each panel.")) {
  # unname(): script 2 passes a NAMED list; named labels would turn
  # c(delta = scen_lab[i]) into "delta.<scenario>" and every panel would
  # fall into one unlabelled facet
  scen_lab <- unname(vapply(nets, function(n)
    if (exists("pretty_scenario")) unname(pretty_scenario(n$scenario)) else n$scenario, ""))
  panel_w <- width / ncol
  Ls <- lapply(seq_along(nets), function(i)
    .flux_network_layers(nets[[i]], "delta", panel_w, min_width, max_width,
                         label_size, node_size, hw, hh, stub_len,
                         panel_names = setNames(scen_lab[i], "delta"),
                         drop_zero = TRUE, drop_inputs = TRUE, y_scale = y_scale))
  L <- list(lim = Ls[[1]]$lim)
  for (nm in c("paths", "labs", "boxes", "stub_txt")) {
    d <- do.call(rbind, lapply(Ls, function(x) { v <- x[[nm]]; v$panel <- as.character(v$panel); v }))
    d$panel <- factor(d$panel, levels = scen_lab)
    if (nm == "paths") d$id <- paste(d$panel, d$id)   # ids unique across scenarios
    L[[nm]] <- d
  }
  .flux_network_ggplot(L, ncol, hw, hh, label_size, title, subtitle,
                       scale = min(1, panel_w / 13) * 1.35, panel_tag = TRUE)
}

# save_flux_network_deltas(): write plot_flux_network_deltas() to a PNG.
# `width` is the WHOLE figure width in inches (default 6.5 = page width); the
# height follows from the layout (about 2 in per panel when stacked, so four
# panels fit a portrait page).
save_flux_network_deltas <- function(nets, file, ncol = 1, width = 6.5, dpi = 300,
                                     hw = 1.45, hh = 0.5, stub_len = 1.1, y_scale = 0.8, ...) {
  p <- plot_flux_network_deltas(nets, ncol = ncol, width = width,
                                hw = hw, hh = hh, stub_len = stub_len, y_scale = y_scale, ...)
  lim  <- .flux_limits(hw, hh, stub_len, drop_inputs = TRUE, y_scale = y_scale)
  pw   <- width / ncol
  ph   <- pw * 0.94 * diff(lim$y) / diff(lim$x)
  nrow <- ceiling(length(nets) / ncol)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file, p, width = width, height = 0.45 + (ph + 0.1) * nrow,
                  dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(p)
}
