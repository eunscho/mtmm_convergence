generate_dat <- function(conditions, condition_number, rep_set, rep) {
  library(simsem)
  set.seed(10000 * condition_number + 100 * rep_set + rep)
  n <- as.integer(conditions[condition_number, 1])
  m <- as.integer(conditions[condition_number, 2]) 
  t <- as.integer(conditions[condition_number, 3]) 
  mload <- as.double(conditions[condition_number, 4]) 
  tload <- as.double(conditions[condition_number, 5]) 
  mcor <- as.double(conditions[condition_number, 6]) 
  tcor <- as.double(conditions[condition_number, 7]) 

  # mloading <- rep(mload, m * t)
  # tloading <- rep(tload, m * t)
  k <- 1
  LY <- bind(LYfree(t, m, k), LYpop(t, m, k, rep(tload, m * t), rep(mload, m * t)))
  tcor <- rep(tcor, t * (t - 1) / 2) # method correlations
  mcor <- rep(mcor, m * (m - 1) / 2) # method correlations

  PS <- binds(PSfree(t, m), PSpop(t, m, tcor, mcor))
  error.cor <- matrix(0, t * m * k, t * m * k)
  diag(error.cor) <- 1
  RTE <- binds(error.cor)
  mod <- model.cfa(LY = LY, PS = PS, RTE = RTE, indLab = item_label(t,m,k), facLab = fac_label(t, m))
  out <- simsem::generate(mod, n)
  
  return(out)
}
