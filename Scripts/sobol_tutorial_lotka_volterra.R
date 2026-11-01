#------------------------------------------------------------------------#
# Sobol' sensitivity analysis, step by step: a predator-prey tutorial ----
#------------------------------------------------------------------------#
# A small, stand-alone version of Scripts/6_sobol_variance_partition.R.
# Same design and the same object/function names, but for a toy
# Lotka-Volterra predator-prey model that runs in about a minute, so each
# piece can be inspected on its own.
#
# The analogy with the soil food-web project:
#   soil model with no animals   <->  prey alone             (the "baseline")
#   soil model with animals      <->  prey + predator        (the "treatment")
#   animal effect on soil C (%)  <->  predator effect on prey density (%)
#   animal parameters (adj_*)    <->  predator parameters (a, e, m)
#
# The steps are the same as the project:
#   1. Model + equilibrium with and without the "animal" (predator)
#   2. Output: % effect of the predator on the prey equilibrium
#   3. Parameter distributions (what we are uncertain about)
#   4. Sobol' design: matrices A, B and A_B^j (every parameter analysed)
#   5. Run the model on every row
#   6. Sobol' indices S (first order) and ST (total effect) + bootstrap CIs
#   7. Figure
#
# Run from anywhere: source("Scripts/sobol_tutorial_lotka_volterra.R")
# Runs in about a minute.
# Needs: deSolve, rootSolve, ggplot2 (all already used by the project).

library(deSolve); library(rootSolve); library(ggplot2)

#------------------------------------------------------------------------#
# Settings ----
#------------------------------------------------------------------------#
# Same names and meaning as in Scripts/6.
sobol_N      <- 1000         # base rows of A and B; total runs = N (k + 2)
sobol_seed   <- 20261029
sobol_boot_B <- 500          # bootstrap replicates for the 95% CIs

#------------------------------------------------------------------------#
# 1. The model ----
#------------------------------------------------------------------------#
# Logistic prey, linear (type I) predation, predator with its own mortality
# and a weak self-limitation. Two state variables, like the pools in
# millennial_model_wplant(); parameters arrive in a list `parms`.
#   dPrey/dt = r Prey (1 - Prey/K) + imm - a Prey Pred
#   dPred/dt = e a Prey Pred - m Pred - d_pred Pred^2
# Parameters:
#   r       prey growth rate                     (prey, "soil" parameter)
#   K       prey carrying capacity               (prey)
#   imm     prey immigration (tiny on purpose)   (prey)
#   a       predator attack rate                 (predator = "animal")
#   e       predator conversion efficiency       (predator)
#   m       predator mortality                   (predator)
#   d_pred  predator self-limitation (tiny)      (predator)
lv_model <- function(time, state, parms) {
  Prey <- state[["Prey"]]
  Pred <- if ("Pred" %in% names(state)) state[["Pred"]] else 0
  with(parms, {
    dPrey <- r * Prey * (1 - Prey / K) + imm - a * Prey * Pred
    dPred <- e * a * Prey * Pred - m * Pred - d_pred * Pred^2
    list(c(dPrey = dPrey, dPred = dPred)[paste0("d", names(state))])
  })
}

# Nominal ("standard") parameter values: the centre of every distribution.
nominal_parms <- list(r = 1.0, K = 100, imm = 0.1,
                      a = 0.02, e = 0.3, m = 0.3, d_pred = 0.001)

# Which parameters belong to the "animal" (predator). In the project these
# are the fitted feeding rates and the animal-effect slopes; here they are
# only flagged in the results table (is_animal).
animal_params <- c("a", "e", "m", "d_pred")

#------------------------------------------------------------------------#
# 2. Equilibrium with and without the predator ----
#------------------------------------------------------------------------#
# setup_scenario_pair() in the project builds two "arms" that share every
# parameter: a TREATMENT (animals on) and a BASELINE (animals off: their
# pools removed from the state). Same idea here.
setup_pair <- function(parms) {
  list(treatment = list(parms = parms, working_state = c(Prey = 50, Pred = 5)),
       baseline  = list(parms = parms, working_state = c(Prey = 50)))
}

# spinup_equilibrium(): integrate until the derivatives stop changing
# (rootSolve::runsteady), exactly as the project does under constant forcing.
spinup_equilibrium <- function(obj) {
  res <- runsteady(y = obj$working_state, func = lv_model, parms = obj$parms,
                   times = c(0, 1e5), stol = 1e-8)
  obj$init_state_spin <- res$y
  obj
}

# The OUTPUT whose variance we partition. Like sobol_model_output() with
# effect_output = "pct": the % change in prey caused by the predator,
#   effect = 100 (Prey_with - Prey_without) / Prey_without.
# It is usually negative (predators reduce prey); it is 0 when the predator
# cannot persist (when m > e a K the predator dies out at equilibrium).
sobol_model_output <- function(parms) {
  pair <- setup_pair(parms)
  b <- spinup_equilibrium(pair$baseline)$init_state_spin
  t <- spinup_equilibrium(pair$treatment)$init_state_spin
  100 * (t[["Prey"]] - b[["Prey"]]) / b[["Prey"]]
}

cat("Nominal predator effect on prey:",
    round(sobol_model_output(nominal_parms), 1), "%\n")

#------------------------------------------------------------------------#
# 3. Parameter distributions ----
#------------------------------------------------------------------------#
# Every parameter gets a distribution describing our uncertainty about it.
# The project uses log-normal-like ranges set by a CV (cv_default) or by
# expert ranges ("main" analysis). Here: each parameter varies uniformly
# between 0.5x and 1.5x its nominal value (i.e. +-50%).
# A distribution is stored as a QUANTILE FUNCTION q(u): it turns a uniform
# number u in (0,1) into a parameter value. That is all the design needs.
dists <- lapply(nominal_parms, function(v) function(u) qunif(u, 0.5 * v, 1.5 * v))
all_params <- names(nominal_parms)

# turn a vector of named values into a full parameter list
as_parms <- function(values) { p <- nominal_parms; p[names(values)] <- as.list(values); p }

#------------------------------------------------------------------------#
# 4. The Sobol' design (sobol_design) ----
#------------------------------------------------------------------------#
# Sobol' asks: what fraction of the output VARIANCE is due to each parameter?
# The Saltelli scheme estimates this from three kinds of sample matrices,
# each N rows x (all parameters) columns, in unit space (uniform 0-1):
#   A      N independent random parameter sets
#   B      N MORE independent random parameter sets
#   A_B^j  copy of A, but with column j taken from B
# Comparing f(A_B^j) with f(A) shows what happens when ONLY parameter j
# changes; comparing it with f(B) shows what happens when everything EXCEPT
# j changes. One A_B^j per parameter, so the cost is N (k + 2) runs.
set.seed(sobol_seed)
p <- length(all_params)
A <- matrix(runif(sobol_N * p), sobol_N, p, dimnames = list(NULL, all_params))
B <- matrix(runif(sobol_N * p), sobol_N, p, dimnames = list(NULL, all_params))

groups <- setNames(as.list(all_params), all_params)   # factor name -> its column(s)
k <- length(groups)
AB <- lapply(groups, function(cols) { M <- A; M[, cols] <- B[, cols]; M })
cat(sprintf("\nSobol' design: k = %d factors, N = %d -> N (k + 2) = %d model runs\n",
            k, sobol_N, sobol_N * (k + 2)))

#------------------------------------------------------------------------#
# 5. Run the model on every row ----
#------------------------------------------------------------------------#
# Each row of a unit-space matrix is mapped through the quantile functions
# into real parameter values, then run with and without the predator.
# (The project does this in parallel with make_sobol_fun()/parLapply.)
to_par <- function(U) {
  X <- U
  for (nm in colnames(U)) X[, nm] <- dists[[nm]](U[, nm])
  X
}
run_rows <- function(U) {
  X <- to_par(U)
  vapply(seq_len(nrow(X)), function(i) sobol_model_output(as_parms(X[i, ])), numeric(1))
}
t0  <- Sys.time()
yA  <- run_rows(A)
yB  <- run_rows(B)
yAB <- sapply(AB, run_rows)          # N x k matrix: column j = f(A_B^j)
cat(sprintf("Runs finished in %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

#------------------------------------------------------------------------#
# 6. Sobol' indices (sobol_indices) ----
#------------------------------------------------------------------------#
# V = total output variance (over A and B together).
# First-order index (Saltelli 2010):
#   S_j  = mean( f(B) * (f(A_B^j) - f(A)) ) / V
#   share of V explained by factor j ALONE (its "main effect").
# Total-effect index (Jansen 1999):
#   ST_j = mean( (f(A) - f(A_B^j))^2 ) / (2 V)
#   share of V that involves j at all: alone + in interactions with others.
#   = the variance that would disappear if j were fixed.
# Reading them:
#   ST_j ~ 0           j does not matter (could be fixed)
#   S_j ~ ST_j         j acts mostly on its own (additively)
#   ST_j >> S_j        j matters mainly through interactions
#   sum(S) < 1         the rest of the variance is interactions
#   small negative S   just sampling noise around 0 (check that its CI spans 0)
# These are ESTIMATES: nothing forces them into [0, 1]. With well-behaved
# output they land close; with heavy-tailed output they can wander far out.
# 95% CIs: resample the N rows with replacement (bootstrap) and recompute.
sobol_indices <- function(yA, yB, yAB, B_boot = 1000, conf = 0.95, seed = NULL) {
  est <- function(i) {                       # i = which rows to use
    V  <- var(c(yA[i], yB[i]))
    S  <- colMeans(yB[i] * (yAB[i, , drop = FALSE] - yA[i])) / V
    ST <- colMeans((yA[i] - yAB[i, , drop = FALSE])^2) / (2 * V)
    rbind(S = S, ST = ST)
  }
  pt <- est(seq_along(yA))
  if (!is.null(seed)) set.seed(seed)
  bs <- replicate(B_boot, est(sample(length(yA), replace = TRUE)))   # 2 x k x B
  lo <- apply(bs, c(1, 2), quantile, (1 - conf) / 2)
  hi <- apply(bs, c(1, 2), quantile, 1 - (1 - conf) / 2)
  data.frame(factor = colnames(pt),
             S = pt["S", ],  S_lo = lo["S", ],  S_hi = hi["S", ],
             ST = pt["ST", ], ST_lo = lo["ST", ], ST_hi = hi["ST", ],
             row.names = NULL)
}
ind <- sobol_indices(yA, yB, yAB, B_boot = sobol_boot_B, seed = sobol_seed)
ind <- ind[order(-ind$ST), ]
ind$is_animal <- ind$factor %in% animal_params   # predator ("animal") parameters
cat("\nSobol' indices of the predator effect on prey (%):\n")
print(data.frame(lapply(ind, function(x) if (is.numeric(x)) round(x, 3) else x)),
      row.names = FALSE)
cat(sprintf("Sum of first-order indices: %.2f (1 - this = share due to interactions)\n",
            sum(ind$S)))

#------------------------------------------------------------------------#
# 7. Figure (plot_sobol_indices) ----
#------------------------------------------------------------------------#
# Bars = indices, whiskers = 95% bootstrap CIs. The gap between ST (total)
# and S (first order) for a factor is its share through interactions.
long <- rbind(
  data.frame(factor = ind$factor, index = "S (first order)",
             value = ind$S, lo = ind$S_lo, hi = ind$S_hi),
  data.frame(factor = ind$factor, index = "ST (total effect)",
             value = ind$ST, lo = ind$ST_lo, hi = ind$ST_hi))
long$factor <- factor(long$factor, levels = rev(ind$factor))
fig <- ggplot(long, aes(x = value, y = factor, fill = index)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.8),
                 height = 0.25) +
  scale_fill_manual(values = c("grey70", "black"), name = NULL) +
  labs(x = "Share of output variance", y = NULL,
       title = "Sobol' indices: predator effect on prey equilibrium (%)",
       subtitle = "Grey = first order (alone); black = total (incl. interactions)") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
print(fig)
ggsave("sobol_tutorial_lotka_volterra.png", fig, width = 7, height = 4.5, bg = "white")
cat("\nFigure saved to sobol_tutorial_lotka_volterra.png (in the working directory).\n")
