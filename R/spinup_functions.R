# ---- One-period integrator ----
one_period <- function(x0, times, model, parms, ...) {
  out <- ode(y = x0, times = times, func = model, parms = parms, ...)
  return(out[nrow(out), -1])  # return final state only
}

# ---- Aitken acceleration (vector-safe) ----
aitken_step <- function(x0, x1, x2, eps = 1e-10) {
  denom <- (x2 - 2*x1 + x0)
  
  # avoid division by ~0
  small <- abs(denom) < eps
  denom[small] <- eps
  
  x_acc <- x0 - (x1 - x0)^2 / denom
  
  return(x_acc)
}

# ---- Main solver ----
poincare_aitken <- function(x_init,
                            times,        # one-period time grid
                            model,
                            parms,
                            tol = 1e-6,
                            max_iter = 100,
                            aitken_start = 3,
                            verbose = TRUE,
                            ...) {
  
  # Initial iterations
  x0 <- x_init
  x1 <- one_period(x0, times, model, parms, ...)
  x2 <- one_period(x1, times, model, parms, ...)
  
  for (iter in 1:max_iter) {
    
    # ---- Apply Aitken after a few iterations ----
    if (iter >= aitken_start) {
      x_guess <- aitken_step(x0, x1, x2)
    } else {
      x_guess <- x2
    }
    
    # propagate one period
    x_new <- one_period(x_guess, times, model, parms, ...)
    
    # convergence check (fixed-point error)
    err <- max(abs(x_new - x_guess))
    
    if (verbose) {
      cat(sprintf("Iter %d | error = %.3e\n", iter, err))
    }
    
    if (err < tol) {
      if (verbose) cat("✅ Converged to periodic solution\n")
      return(list(
        x = x_new,
        iterations = iter,
        error = err
      ))
    }
    
    # shift iterates
    x0 <- x1
    x1 <- x2
    x2 <- x_new
  }
  
  warning("⚠️ Did not fully converge")
  
  return(list(
    x = x2,
    iterations = max_iter,
    error = err
  ))
}