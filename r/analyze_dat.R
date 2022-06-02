analyze_dat <- function(conditions, condition_number, rep_set, rep, data) {
  CUTOFF <- .001
  library(lavaan)
  library(semTools)
  n <- as.integer(conditions[condition_number, 1])
  m <- as.integer(conditions[condition_number, 2]) 
  t <- as.integer(conditions[condition_number, 3]) 
  mload <- as.double(conditions[condition_number, 4]) 
  tload <- as.double(conditions[condition_number, 5]) 
  mcor <- as.double(conditions[condition_number, 6]) 
  tcor <- as.double(conditions[condition_number, 7]) 

  ###############################################################################
  # Make Models
  ###############################################################################
  # ct: correlated trait, common module
  m_ct <- vector("character")
  for (i in 1:t) {
    m_ct <- paste0(m_ct, 'T', i, ' =~ NA * t', i, 'm1')
    for (j in 2:m) {
      m_ct <- paste0(m_ct, '+ t', i, 'm', j)
    }
    m_ct <- paste0(m_ct, ';')
    m_ct <- paste0(m_ct, 'T', i, ' ~~ 1 * T', i, ';')
  }
  for (i in 1:(t - 1)) {
    for (j in (i + 1):t) {
      m_ct <- paste0(m_ct, 'T', i, ' ~~ a', i, j, ' * T', j, ';')
    }
  }
  #ctcm: ct-correlated methods
  m_ctcm <- m_ct
  for (j in 1:m) {
    m_ctcm <- paste0(m_ctcm, 'M', j, ' =~ NA * t1m', j)
    for (i in 2:t) {
      m_ctcm <- paste0(m_ctcm, '+ t', i, 'm', j)
    }
    m_ctcm <- paste0(m_ctcm, ';')
    m_ctcm <- paste0(m_ctcm, 'M', j, ' ~~ 1 * M', j, ';')
  }
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      m_ctcm <- paste0(m_ctcm, 'M', i, ' ~~ b', i, j, ' * M', j, ';')
    }
  }
  for (i in 1:t) {
    for (j in 1:m) {
      m_ctcm <- paste0(m_ctcm, 'T', i, ' ~~ 0 * M', j, ';')
      #m_ctcm <- paste0(m_ctcm, 't', i, 'm', j, ' ~~ c', i, j, ' * t', i, 'm', j, ';')
    }
  }
  m_ctcmuc <- m_ctcm
  # Rindskopf's reparametrization
  m_ctcmr <- m_ctcm
  for (i in 1:t) {
    for (j in 1:m) {
      m_ctcmr <- paste0(m_ctcmr, 'P', i, j, ' =~ NA * t', i, 'm', j, ';')
      m_ctcmr <- paste0(m_ctcmr, 'P', i, j, ' ~~ 1 * P', i, j, ';')
      #m_ctcmr <- paste0(m_ctcmr, 'c', i, j, ' == 0', ';')
      m_ctcmr <- paste0(m_ctcmr, 't', i, 'm', j, ' ~~ 0 * t', i, 'm', j, ';')
      for (k in 1:t) {
        m_ctcmr <- paste0(m_ctcmr, 'P', i, j, ' ~~ 0 * T', k, ';')
      }
      for (l in 1:m) {
        m_ctcmr <- paste0(m_ctcmr, 'P', i, j, ' ~~ 0 * M', l, ';')  
      }
    }
  }
  for (i in 1:t) {
    for (j in 1:m) {
      for (k in i:t) {
        for (l in 1:m) {
          if ((i == k) & (l > j) | k > i)  {
            m_ctcmr <- paste0(m_ctcmr, 'P', i, j, ' ~~ 0 * P', k, l, ';')
          }
        }
      }
    }
  }
  # partial constraining
  m_ctcmpc <- m_ctcm
  for (i in 1:t) {
    for (j in 1:m) {
      m_ctcmpc <- paste0(m_ctcmpc, 't', i, 'm', j, ' ~~ c', i, j, ' * t', i, 'm', j, ';')
      m_ctcmpc <- paste0(m_ctcmpc, 'c', i, j, ' > 0;')
    }
  }
  # full constraining
  m_ctcmfc <- m_ctcmpc
  for (i in 1:(t - 1)) {
    for (j in (i + 1):t) {
      m_ctcmfc <- paste0(m_ctcmfc, 'a', i, j, ' > -1;')
      m_ctcmfc <- paste0(m_ctcmfc, 'a', i, j, ' < 1;')
    }
  }
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      m_ctcmfc <- paste0(m_ctcmfc, 'b', i, j, ' > -1;')
      m_ctcmfc <- paste0(m_ctcmfc, 'b', i, j, ' < 1;')
    }
  }

  # ctcu: ct-correlated uniqueness
  m_ctcu <- m_ct
  for (j in 1:m) {
    for (i in 1:(t - 1)) {
      for (k in (i + 1):t)
        m_ctcu <- paste0(m_ctcu, 't', i, 'm', j, ' ~~ t', k, 'm', j, '; ')
    }
  }
  #ctmin: ct-correlated (m-1)methods
  m_ctmin <- m_ct
  for (j in 2:m) { # not 1, because m-1 
    m_ctmin <- paste0(m_ctmin, 'M', j, ' =~ NA * t1m', j)
    for (i in 2:t) {
      m_ctmin <- paste0(m_ctmin, '+ t', i, 'm', j)
    }
    m_ctmin <- paste0(m_ctmin, ';')
    m_ctmin <- paste0(m_ctmin, 'M', j, ' ~~ 1 * M', j, ';')
  }
  for (i in 1:t) {
    for (j in 2:m) {
      m_ctmin <- paste0(m_ctmin, 'T', i, ' ~~ 0 * M', j, ';')
      #m_ctmin <- paste0(m_ctmin, 't', i, 'm', j, ' ~~ c', i, j, ' * t', i, 'm', j, ';')
    }
  }
  ###############################################################################
  # Estimate Models & Fit indices
  ############################################################################### 
  # estimate output
  o_ctcmuc <- cfa(m_ctcmuc, data, control = list(iter.max = 250))
  o_ctcmr <- cfa(m_ctcmr, data, control = list(iter.max = 250))
  o_ctcmpc <- cfa(m_ctcmpc, data, control = list(iter.max = 250))
  o_ctcmfc <- cfa(m_ctcmfc, data, control = list(iter.max = 250))
  o_ctcu <- cfa(m_ctcu, data, control = list(iter.max = 250))
  o_ctmin <- cfa(m_ctmin, data, control = list(iter.max = 250))
  # convergence
  c_ctcmuc <- lavInspect(o_ctcmuc, "converged")
  c_ctcmr <- lavInspect(o_ctcmr, "converged")
  c_ctcmpc <- lavInspect(o_ctcmpc, "converged")
  c_ctcmfc <- lavInspect(o_ctcmfc, "converged")
  c_ctcu <- lavInspect(o_ctcu, "converged")
  c_ctmin <- lavInspect(o_ctmin, "converged")
  # fit
  if (c_ctcmuc) {
    f_ctcmuc <- lavInspect(o_ctcmuc, "fit")
  } else {
    f_ctcmuc <- NA
  }
  if (c_ctcmr) {
    f_ctcmr <- lavInspect(o_ctcmr, "fit")
  } else {
    f_ctcmr <- NA
  }
  if (c_ctcmpc) {
    f_ctcmpc <- lavInspect(o_ctcmpc, "fit")
  } else {
    f_ctcmpc <- NA
  }
  if (c_ctcmfc) {
    f_ctcmfc <- lavInspect(o_ctcmfc, "fit")
  } else {
    f_ctcmfc <- NA
  }
  if (c_ctcu) {
    f_ctcu <- lavInspect(o_ctcu, "fit")
  } else {
    f_ctcu <- NA
  }
  if (c_ctmin) {
    f_ctmin <- lavInspect(o_ctmin, "fit")
  } else {
    f_ctmin <- NA
  }
  # chisquare
  chi_ctcmuc <- f_ctcmuc['chisq']
  chi_ctcmr <- f_ctcmr['chisq']
  chi_ctcmpc <- f_ctcmpc['chisq']
  chi_ctcmfc <- f_ctcmfc['chisq']
  chi_ctcu <- f_ctcu['chisq']
  chi_ctmin <- f_ctmin['chisq']
  # degree of freedom
  df_ctcmuc <- f_ctcmuc['df']
  df_ctcmr <- f_ctcmr['df']
  df_ctcmpc <- f_ctcmpc['df']
  df_ctcmfc <- f_ctcmfc['df']
  df_ctcu <- f_ctcu['df']
  df_ctmin <- f_ctmin['df']
  # cfi
  cfi_ctcmuc <- f_ctcmuc['cfi']
  cfi_ctcmr <- f_ctcmr['cfi']
  cfi_ctcmpc <- f_ctcmpc['cfi']
  cfi_ctcmfc <- f_ctcmfc['cfi']
  cfi_ctcu <- f_ctcu['cfi']
  cfi_ctmin <- f_ctmin['cfi']
  # tli
  tli_ctcmuc <- f_ctcmuc['tli']
  tli_ctcmr <- f_ctcmr['tli']
  tli_ctcmpc <- f_ctcmpc['tli']
  tli_ctcmfc <- f_ctcmfc['tli']
  tli_ctcu <- f_ctcu['tli']
  tli_ctmin <- f_ctmin['tli']
  # rmsea
  rmsea_ctcmuc <- f_ctcmuc['rmsea']
  rmsea_ctcmr <- f_ctcmr['rmsea']
  rmsea_ctcmpc <- f_ctcmpc['rmsea']
  rmsea_ctcmfc <- f_ctcmfc['rmsea']
  rmsea_ctcu <- f_ctcu['rmsea']
  rmsea_ctmin <- f_ctmin['rmsea']
  ###############################################################################
  # Determining inadmissible solutions & other performance
  ############################################################################### 
  # finding parameter estimates
  e_ctcmuc <- lavInspect(o_ctcmuc, "est")
  e_ctcmr <- lavInspect(o_ctcmr, "est")
  e_ctcmpc <- lavInspect(o_ctcmpc, "est")
  e_ctcmfc <- lavInspect(o_ctcmfc, "est")
  e_ctcu <- lavInspect(o_ctcu, "est")
  e_ctmin <- lavInspect(o_ctmin, "est")
  # inadmissible solutions
  if (min(diag(e_ctcmuc$theta)) < 0 | outrange(e_ctcmuc$psi)) {
    in_ctcmuc <- TRUE
  } else {
    in_ctcmuc <- FALSE
  }
  if (min(diag(e_ctcmr$theta)) < 0 | outrange(e_ctcmr$psi)) {
    in_ctcmr <- TRUE
  } else {
    in_ctcmr <- FALSE
  }
  if (min(diag(e_ctcmpc$theta)) < 0 | outrange(e_ctcmpc$psi)) {
    in_ctcmpc <- TRUE
  } else {
    in_ctcmpc <- FALSE
  }
  if (min(diag(e_ctcmfc$theta)) < 0 | outrange(e_ctcmfc$psi)) {
    in_ctcmfc <- TRUE
  } else {
    in_ctcmfc <- FALSE
  }
  if (min(diag(e_ctcu$theta)) < 0 | outrange(e_ctcu$psi) | outrange(e_ctcu$theta)) {
    in_ctcu <- TRUE
  } else {
    in_ctcu <- FALSE
  }
  if (min(diag(e_ctmin$theta)) < 0 | outrange(e_ctmin$psi)) {
    in_ctmin <- TRUE
  } else {
    in_ctmin <- FALSE
  }
  # boundary solutions
  bnd_ctcmuc <- FALSE
  bnd_ctcmr <- FALSE
  if (min(diag(e_ctcmpc$theta)) < CUTOFF) {
    bnd_ctcmpc <- TRUE
  } else {
    bnd_ctcmpc <- FALSE
  }
  if (min(diag(e_ctcmfc$theta)) < CUTOFF | corone(e_ctcmfc$psi, CUTOFF)) {
    bnd_ctcmfc <- TRUE
  } else {
    bnd_ctcmfc <- FALSE
  }
  bnd_ctcu <- FALSE
  bnd_ctmin <- FALSE
  # rindskopf's reparameterization vs. partial constraining equality
  if (is.na(chi_ctcmr) | is.na(chi_ctcmpc)) {
    equal <- NA
  } else if (chi_ctcmr == chi_ctcmpc) {
    equal <- TRUE
  } else {
    equal <- FALSE
  }
  # bias in trait loadings
  tlbias_ctcmuc <- sum(e_ctcmuc$lambda[, 1:t]) - tload * t * m
  tlbias_ctcmr <- sum(e_ctcmr$lambda[, 1:t]) - tload * t * m
  tlbias_ctcmpc <- sum(e_ctcmpc$lambda[, 1:t]) - tload * t * m
  tlbias_ctcmfc <- sum(e_ctcmfc$lambda[, 1:t]) - tload * t * m
  tlbias_ctcu <- sum(e_ctcu$lambda[, 1:t]) - tload * t * m
  tlbias_ctmin <- sum(e_ctmin$lambda[, 1:t]) - tload * t * m
  #mae in trait loadings
  tlmae_ctcmuc <- mean(abs(rowSums(e_ctcmuc$lambda[, 1:t]) - tload))
  tlmae_ctcmr <- mean(abs(rowSums(e_ctcmr$lambda[, 1:t]) - tload))
  tlmae_ctcmpc <- mean(abs(rowSums(e_ctcmpc$lambda[, 1:t]) - tload))
  tlmae_ctcmfc <- mean(abs(rowSums(e_ctcmfc$lambda[, 1:t]) - tload))
  tlmae_ctcu <- mean(abs(rowSums(e_ctcu$lambda[, 1:t]) - tload))
  tlmae_ctmin <- mean(abs(rowSums(e_ctmin$lambda[, 1:t]) - tload))
  # estimated trait correlations
  etvar_ctcmuc <- cov2cor(e_ctcmuc$psi[1:m, 1:m])
  etvar_ctcmr <- cov2cor(e_ctcmr$psi[1:m, 1:m])
  etvar_ctcmpc <- cov2cor(e_ctcmpc$psi[1:m, 1:m])
  etvar_ctcmfc <- cov2cor(e_ctcmfc$psi[1:m, 1:m])
  etvar_ctcu <- cov2cor(e_ctcu$psi[1:m, 1:m])
  etvar_ctmin <- cov2cor(e_ctmin$psi[1:m, 1:m])
  # bias of trait correlations
  tvbias_ctcmuc <- sum(etvar_ctcmuc[lower.tri(etvar_ctcmuc)] - tcor)
  tvbias_ctcmr <- sum(etvar_ctcmr[lower.tri(etvar_ctcmr)] - tcor)
  tvbias_ctcmpc <- sum(etvar_ctcmpc[lower.tri(etvar_ctcmpc)] - tcor)
  tvbias_ctcmfc <- sum(etvar_ctcmfc[lower.tri(etvar_ctcmfc)] - tcor)
  tvbias_ctcu <- sum(etvar_ctcu[lower.tri(etvar_ctcu)] - tcor)
  tvbias_ctmin <- sum(etvar_ctmin[lower.tri(etvar_ctmin)] - tcor)
  # mae of trait correlations
  tvmae_ctcmuc <- mean(abs(etvar_ctcmuc[lower.tri(etvar_ctcmuc)] - tcor))
  tvmae_ctcmr <- mean(abs(etvar_ctcmr[lower.tri(etvar_ctcmr)] - tcor))
  tvmae_ctcmpc <- mean(abs(etvar_ctcmpc[lower.tri(etvar_ctcmpc)] - tcor))
  tvmae_ctcmfc <- mean(abs(etvar_ctcmfc[lower.tri(etvar_ctcmfc)] - tcor))
  tvmae_ctcu <- mean(abs(etvar_ctcu[lower.tri(etvar_ctcu)] - tcor))
  tvmae_ctmin <- mean(abs(etvar_ctmin[lower.tri(etvar_ctmin)] - tcor))
    
  # the ratio of trait variance
  # tvar_ctcmuc <- simtvar(e_ctcmuc$lambda[, 1:t], 
  #                        e_ctcmuc$lambda[, (t+1):(t+m)],
  #                        e_ctcmuc$psi[1:t, 1:t],
  #                        e_ctcmuc$psi[(t + 1):(t + m), (t + 1):(t + m)],
  #                        e_ctcmuc$theta)
  # tvar_ctcmr <- simtvar(e_ctcmr$lambda[, 1:t], 
  #                        e_ctcmr$lambda[, (t+1):(t+m)],
  #                        e_ctcmr$psi[1:t, 1:t],
  #                        e_ctcmr$psi[(t + 1):(t + m), (t + 1):(t + m)],
  #                        e_ctcmr$theta)
  # tvar_ctcmpc <- simtvar(e_ctcmpc$lambda[, 1:t], 
  #                       e_ctcmpc$lambda[, (t+1):(t+m)],
  #                       e_ctcmpc$psi[1:t, 1:t],
  #                       e_ctcmpc$psi[(t + 1):(t + m), (t + 1):(t + m)],
  #                       e_ctcmpc$theta)
  # tvar_ctcmfc <- simtvar(e_ctcmfc$lambda[, 1:t], 
  #                       e_ctcmfc$lambda[, (t+1):(t+m)],
  #                       e_ctcmfc$psi[1:t, 1:t],
  #                       e_ctcmfc$psi[(t + 1):(t + m), (t + 1):(t + m)],
  #                       e_ctcmfc$theta)
  # tvar_ctcu <- simtvar(tloadmat = e_ctcu$lambda[, 1:t], 
  #                        tcov = e_ctcu$psi[1:t, 1:t],
  #                        errorcov = e_ctcu$theta,
  #                        ctcu = TRUE)
  # tvar_ctmin <- simtvar(e_ctmin$lambda[, 1:t], 
  #                        e_ctmin$lambda[, (t+1):(t+(m - 1))],
  #                        e_ctmin$psi[1:t, 1:t],
  #                        e_ctmin$psi[(t + 1):(t + (m - 1)), (t + 1):(t + (m - 1))],
  #                        e_ctmin$theta)
  
  ###############################################################################
  # Report output
  ###############################################################################
  out <- data.frame(condition_number = condition_number,
                    rep_set = rep_set,
                    rep = rep,
                    n = n,
                    m = m,
                    t = t,
                    mload = mload,
                    tload = tload, 
                    mcor = mcor,
                    tcor = tcor,
                    equal = equal,
                    c_ctcmuc = c_ctcmuc,
                    c_ctcmr = c_ctcmr,
                    c_ctcmpc = c_ctcmpc,
                    c_ctcmfc = c_ctcmfc,
                    c_ctcu = c_ctcu,
                    c_ctmin = c_ctmin,
                    chi_ctcmuc = chi_ctcmuc,
                    chi_ctcmr = chi_ctcmr,
                    chi_ctcmpc = chi_ctcmpc,
                    chi_ctcmfc = chi_ctcmfc,
                    chi_ctcu = chi_ctcu,
                    chi_ctmin = chi_ctmin,
                    df_ctcmuc = df_ctcmuc,
                    df_ctcmr = df_ctcmr,
                    df_ctcmpc = df_ctcmpc,
                    df_ctcmfc = df_ctcmfc,
                    df_ctcu = df_ctcu,
                    df_ctmin = df_ctmin,
                    cfi_ctcmuc = cfi_ctcmuc,
                    cfi_ctcmr = cfi_ctcmr,
                    cfi_ctcmpc = cfi_ctcmpc,
                    cfi_ctcmfc = cfi_ctcmfc,
                    cfi_ctcu = cfi_ctcu,
                    cfi_ctmin = cfi_ctmin,
                    tli_ctcmuc = tli_ctcmuc,
                    tli_ctcmr = tli_ctcmr,
                    tli_ctcmpc = tli_ctcmpc,
                    tli_ctcmfc = tli_ctcmfc,
                    tli_ctcu = tli_ctcu,
                    tli_ctmin = tli_ctmin,
                    rmsea_ctcmuc = rmsea_ctcmuc,
                    rmsea_ctcmr = rmsea_ctcmr,
                    rmsea_ctcmpc = rmsea_ctcmpc,
                    rmsea_ctcmfc = rmsea_ctcmfc,
                    rmsea_ctcu = rmsea_ctcu,
                    rmsea_ctmin = rmsea_ctmin,
                    in_ctcmuc = in_ctcmuc,
                    in_ctcmr = in_ctcmr,
                    in_ctcmpc = in_ctcmpc,
                    in_ctcmfc = in_ctcmfc,
                    in_ctcu = in_ctcu,
                    in_ctmin = in_ctmin,
                    bnd_ctcmuc = bnd_ctcmuc,
                    bnd_ctcmr = bnd_ctcmr,
                    bnd_ctcmpc = bnd_ctcmpc,
                    bnd_ctcmfc = bnd_ctcmfc,
                    bnd_ctcu = bnd_ctcu,
                    bnd_ctmin = bnd_ctmin,
                    tlbias_ctcmuc = tlbias_ctcmuc,
                    tlbias_ctcmr = tlbias_ctcmr,
                    tlbias_ctcmpc = tlbias_ctcmpc,
                    tlbias_ctcmfc = tlbias_ctcmfc,
                    tlbias_ctcu = tlbias_ctcu,
                    tlbias_ctmin = tlbias_ctmin,
                    tlmae_ctcmuc = tlmae_ctcmuc,
                    tlmae_ctcmr = tlmae_ctcmr,
                    tlmae_ctcmpc = tlmae_ctcmpc,
                    tlmae_ctcmfc = tlmae_ctcmfc,
                    tlmae_ctcu = tlmae_ctcu,
                    tlmae_ctmin = tlmae_ctmin,
                    tvbias_ctcmuc = tvbias_ctcmuc,
                    tvbias_ctcmr = tvbias_ctcmr,
                    tvbias_ctcmpc = tvbias_ctcmpc,
                    tvbias_ctcmfc = tvbias_ctcmfc,
                    tvbias_ctcu = tvbias_ctcu,
                    tvbias_ctmin = tvbias_ctmin,
                    tvmae_ctcmuc = tvmae_ctcmuc,
                    tvmae_ctcmr = tvmae_ctcmr,
                    tvmae_ctcmpc = tvmae_ctcmpc,
                    tvmae_ctcmfc = tvmae_ctcmfc,
                    tvmae_ctcu = tvmae_ctcu,
                    tvmae_ctmin = tvbias_ctmin
  )
  
   return(out)
}