# Rebuild the tree growth forcing for each model to make the simulations simpler:

library(pacman)
p_load(deSolve, FME, tidyverse, yaml)
verbose = F

# ---- Load models ----
source("R/tree_monomolecular.R")
source("R/millennial_model_earthworm.R")
source("R/millennial_model.R")

# ---- Load config utilities ----
source("R/make_tree_forcing.R")
source("R/derive_millennial_parms.R")

parms  <- yaml::read_yaml("config/common.yml")
parms  <- modifyList(parms, yaml::read_yaml("config/tree_monomolecular.yml"))
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))

parms$tree_forcing <- make_tree_forcing(parms)

parms$GPPmax_herb = 4.0
parms$Q10 = 2
parms$Tref = 10
parms$theta_opt = 0.30          # Optimal soil moisture for plant productivity

parms$maint_resp =  0.002        # Maintenance respiration (g C g C-1 d-1)
parms$growth_resp =  0.25        # Fraction of GPP lost to growth respiration

tibble(data.frame(t(sapply(1:365*100, parms$tree_forcing)))) %>%
  mutate(doy = 1:365*100) %>%
  ggplot(aes(x = doy, y = B_tree)) + geom_line()

state <- c(C_shoot = 0, C_wood =0, C_root = 0)
parms

tgm <- function(time, state, parms){
  
  with(c(state, parms), {
  
    # --------------------------------------------------
    # Time + climate
    # --------------------------------------------------
    doy <- (time %% 365) + 1
    
    Temp <- MAT + T_amp * sin(2 * pi * (doy - 110) / 365)
    theta <- MAtheta + theta_amp * cos(4 * pi * (doy - 110) / 365)
    
    # --------------------------------------------------
    # Climate scalars for plant productivity
    # --------------------------------------------------
    f_T     <- Q10 ^ ((Temp - Tref) / 10)
    f_theta <- pmin(1, theta / theta_opt)
    
    # --------------------------------------------------
    # Herbaceous plant carbon fluxes
    # --------------------------------------------------
    GPP_herb <- GPPmax_herb * f_T * f_theta
    
    Ra_herb <- maint_resp * (C_shoot + C_root) +
      growth_resp * GPP_herb
    
    NPP_herb <- pmax(0, GPP_herb - Ra_herb)
    
    shoot_growth <- a_leaf * NPP_herb
    wood_growth  <- a_wood * NPP_herb
    root_growth  <- a_root * NPP_herb
    
    # --------------------------------------------------
    # Seasonal litterfall (shoots)
    # --------------------------------------------------
    wrapped_pdf <- function(t, mu, sigma, L = 365, K = 1) {
      ks <- (-K):K
      rowSums(sapply(ks, function(k)
        dnorm(t - mu - k * L, sd = sigma)))
    }
    
    pdf_vals  <- wrapped_pdf(1:365, litter_peak_doy, litter_width_d)
    prob_vals <- pdf_vals / sum(pdf_vals)
    
    litterfall <- k_litterfall_ann * prob_vals[doy] * C_shoot
    
    # --------------------------------------------------
    # Continuous plant losses
    # --------------------------------------------------
    leaf_mortality <- k_mort_leaf * C_shoot
    wood_mortality <- k_mort_wood * C_wood
    
    root_mortality <- k_mort_root * C_root
    exudates       <- k_exudate* C_root
    
  dC_shoot <- shoot_growth - litterfall - leaf_mortality
  dC_wood <- wood_growth - wood_mortality
  dC_root  <- root_growth  - root_mortality - exudates
  
  list(
    c(
      dC_shoot, dC_wood, dC_root
    )
  )
  }
  )
}

testout = ode(1:365*100, y = state, func = tgm, parms = parms)

tibble(data.frame(testout)) %>%
  pivot_longer(!time) %>%
  ggplot(aes(x = time, y= value)) + geom_line() + facet_wrap(.~name, scales = "free_y")

tibble(data.frame(testout)) %>%
  mutate(C_tree = C_shoot + C_wood + C_root) %>%
  pivot_longer(!time) %>%
  ggplot(aes(x = time, y= value)) + geom_line() + facet_wrap(.~name, scales = "free_y")
