## Using Wang et al 2009; Soil Biology and Biochemistry

source("R/millennial_model.R")

## - Units: All fluxes in g C m^-2 d^-1; states in g C m^-2; moisture in m^3 m^-3; temperature in °C.

y0       <- c(Litter  = 200, 
              CWD     = 1000, 
              Organic = 2000, 
              DOM = 10, 
              MIC = 0.1, 
              P   = 400,  
              L   = 10,    
              A   = 1000,   
              M   = 3000,  
              B   = 40)


# Run model to equilibrium:
parms["force_eqm"] = 1

# Calculate the equilibria across parameter ranges:

sens <- data.frame(
  Parameter = c("k_pa"),
  Baseline = c(0.018),
  High = c(0.018*1.75),
  Low = c(0.018*1.25)
)

# Millennial model sensitivity: 

results <- vector(mode = "list", 1)
names(results) = sens$Parameter

BASELINE = rootSolve::stode(y = y0,
                            func = millennial_model,
                            parms = parms)$y

for(curpar in 1){
  pargrad <- seq(sens[curpar,4], sens[curpar,3], length = 10)
  
  results[[curpar]] <- t(sapply(pargrad, FUN = function(X, nnn = curpar){
    pcur = parms
    pcur[sens[nnn,1]] = X
    rootSolve::stode(y = y0,
                     func = millennial_model,
                     parms = pcur)$y
  }))
  
  results[[curpar]] = as.data.frame(results[[curpar]])
  
  results[[curpar]] <- cbind(results[[curpar]], par = pargrad, parname = names(results)[[curpar]])
}

results$k_pa %>%
  pivot_longer(!par & !parname) %>%
  bind_rows(
    tibble(par = 0.018, parname = "k_pa", name = names(BASELINE), value = BASELINE)
  ) %>%
  ggplot(aes(x = par, y = value, fill = name)) + 
  geom_area() +
  geom_vline(xintercept = 0.018, linetype = 2) +
  geom_vline(xintercept = c(0.0225, 0.0315), linetype = 3) +
  ylab("Pool size (g C m^-2)") + xlab("k_pa") + scale_fill_discrete(name = "Soil pool")
