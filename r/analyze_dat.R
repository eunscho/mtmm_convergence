analyze_dat <- function(conditions, condition_number, rep_set, rep, data) {
  TOLERANCE <- .0001
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
  m_ctm <- m_ct
  for (j in 1:m) {
    m_ctm <- paste0(m_ctm, 'M', j, ' =~ NA * t1m', j)
    for (i in 2:t) {
      m_ctm <- paste0(m_ctm, '+ t', i, 'm', j)
    }
    m_ctm <- paste0(m_ctm, ';')
  }
  m_ctcm <- modelcm(m_ctm, t, m)
  # ctum - ct-uncorrelated model
  m_ctum <- m_ctcm
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      if (!(i == j)) {
        m_ctum <- paste0(m_ctum, 'b', i, j, ' == 0;')
      }
    }
  }
  # ct - constrained model (fixing all method loadings equal)
  m_ctcom <- m_ct
  for (j in 1:m) {
    m_ctcom <- paste0(m_ctcom, 'M', j, ' =~ NA * t1m', j)
    for (i in 2:t) {
      m_ctcom <- paste0(m_ctcom, ' + equal("M', j, ' =~ t1m', j, '") * t', i, 'm', j)
    }
    m_ctcom <- paste0(m_ctcom, '; ')
  }
  m_ctcom <- modelcm(m_ctcom, t, m)
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
  o_ctum <- cfa(m_ctum, data, control = list(iter.max = 250))
  o_ctcom <- cfa(m_ctcom, data, control = list(iter.max = 250))
  o_ctcu <- cfa(m_ctcu, data, control = list(iter.max = 250))
  o_ctmin <- cfa(m_ctmin, data, control = list(iter.max = 250))
  # convergence
  c_ctcmuc <- lavInspect(o_ctcmuc, "converged")
  c_ctcmr <- lavInspect(o_ctcmr, "converged")
  c_ctcmpc <- lavInspect(o_ctcmpc, "converged")
  c_ctcmfc <- lavInspect(o_ctcmfc, "converged")
  c_ctum <- lavInspect(o_ctum, "converged")
  c_ctcom <- lavInspect(o_ctcom, "converged")
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
  if (c_ctum) {
    f_ctum <- lavInspect(o_ctum, "fit")
  } else {
    f_ctum <- NA
  }
  if (c_ctcom) {
    f_ctcom <- lavInspect(o_ctcom, "fit")
  } else {
    f_ctcom <- NA
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
  chi_ctum <- f_ctum['chisq']
  chi_ctcom <- f_ctcom['chisq']
  chi_ctcu <- f_ctcu['chisq']
  chi_ctmin <- f_ctmin['chisq']
  # degree of freedom
  df_ctcmuc <- f_ctcmuc['df']
  df_ctcmr <- f_ctcmr['df']
  df_ctcmpc <- f_ctcmpc['df']
  df_ctcmfc <- f_ctcmfc['df']
  df_ctum <- f_ctum['df']
  df_ctcom <- f_ctcom['df']
  df_ctcu <- f_ctcu['df']
  df_ctmin <- f_ctmin['df']
  # cfi
  cfi_ctcmuc <- f_ctcmuc['cfi']
  cfi_ctcmr <- f_ctcmr['cfi']
  cfi_ctcmpc <- f_ctcmpc['cfi']
  cfi_ctcmfc <- f_ctcmfc['cfi']
  cfi_ctum <- f_ctum['cfi']
  cfi_ctcom <- f_ctcom['cfi']
  cfi_ctcu <- f_ctcu['cfi']
  cfi_ctmin <- f_ctmin['cfi']
  # tli
  tli_ctcmuc <- f_ctcmuc['tli']
  tli_ctcmr <- f_ctcmr['tli']
  tli_ctcmpc <- f_ctcmpc['tli']
  tli_ctcmfc <- f_ctcmfc['tli']
  tli_ctum <- f_ctum['tli']
  tli_ctcom <- f_ctcom['tli']
  tli_ctcu <- f_ctcu['tli']
  tli_ctmin <- f_ctmin['tli']
  # rmsea
  rmsea_ctcmuc <- f_ctcmuc['rmsea']
  rmsea_ctcmr <- f_ctcmr['rmsea']
  rmsea_ctcmpc <- f_ctcmpc['rmsea']
  rmsea_ctcmfc <- f_ctcmfc['rmsea']
  rmsea_ctum <- f_ctum['rmsea']
  rmsea_ctcom <- f_ctcom['rmsea']
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
  e_ctum <- lavInspect(o_ctum, "est")
  e_ctcom <- lavInspect(o_ctcom, "est")
  e_ctcu <- lavInspect(o_ctcu, "est")
  e_ctmin <- lavInspect(o_ctmin, "est")
  # inadmissible solutions
  if (min(diag(e_ctcmuc$theta)) + TOLERANCE < 0 | outrange(e_ctcmuc$psi, TOLERANCE)) {
    in_ctcmuc <- TRUE
  } else {
    in_ctcmuc <- FALSE
  }
  if (min(diag(e_ctcmr$theta)) + TOLERANCE < 0 | outrange(e_ctcmr$psi, TOLERANCE)) {
    in_ctcmr <- TRUE
  } else {
    in_ctcmr <- FALSE
  }
  if (min(diag(e_ctcmpc$theta)) + TOLERANCE < 0 | outrange(e_ctcmpc$psi, TOLERANCE)) {
    in_ctcmpc <- TRUE
  } else {
    in_ctcmpc <- FALSE
  }
  if (min(diag(e_ctcmfc$theta)) + TOLERANCE < 0 | outrange(e_ctcmfc$psi, TOLERANCE)) {
    in_ctcmfc <- TRUE
  } else {
    in_ctcmfc <- FALSE
  }
  if (min(diag(e_ctum$theta)) + TOLERANCE < 0 | outrange(e_ctum$psi, TOLERANCE)) {
    in_ctum <- TRUE
  } else {
    in_ctum <- FALSE
  }
  if (min(diag(e_ctcom$theta)) + TOLERANCE < 0 | outrange(e_ctcom$psi, TOLERANCE)) {
    in_ctcom <- TRUE
  } else {
    in_ctcom <- FALSE
  }
  if (min(diag(e_ctcu$theta)) + TOLERANCE < 0 | outrange(e_ctcu$psi, TOLERANCE) | outrange(e_ctcu$theta, TOLERANCE)) {
    in_ctcu <- TRUE
  } else {
    in_ctcu <- FALSE
  }
  if (min(diag(e_ctmin$theta)) + TOLERANCE < 0 | outrange(e_ctmin$psi, TOLERANCE)) {
    in_ctmin <- TRUE
  } else {
    in_ctmin <- FALSE
  }
  # boundary solutions
  bnd_ctcmuc <- FALSE
  bnd_ctcmr <- FALSE
  if (min(diag(e_ctcmpc$theta)) < TOLERANCE) {
    bnd_ctcmpc <- TRUE
  } else {
    bnd_ctcmpc <- FALSE
  }
  if (min(diag(e_ctcmfc$theta)) < TOLERANCE | corone(e_ctcmfc$psi, TOLERANCE)) {
    bnd_ctcmfc <- TRUE
  } else {
    bnd_ctcmfc <- FALSE
  }
  bnd_ctum <- FALSE
  bnd_ctcom <- FALSE
  bnd_ctcu <- FALSE
  bnd_ctmin <- FALSE
  # rindskopf's reparameterization vs. partial constraining equality
  if (is.na(chi_ctcmr) | is.na(chi_ctcmpc)) {
   r_pc_equal <- NA
  } else {
   if(abs(chi_ctcmr - chi_ctcmpc) < TOLERANCE) {
     r_pc_equal <- TRUE
   } else
     r_pc_equal <- FALSE
  }

  # bias in trait loadings
  tlbias_ctcmuc <- mean(e_ctcmuc$lambda[, 1:t]) - tload
  tlbias_ctcmr <- mean(e_ctcmr$lambda[, 1:t]) - tload
  tlbias_ctcmpc <- mean(e_ctcmpc$lambda[, 1:t]) - tload
  tlbias_ctcmfc <- mean(e_ctcmfc$lambda[, 1:t]) - tload
  tlbias_ctum <- mean(e_ctum$lambda[, 1:t]) - tload
  tlbias_ctcom <- mean(e_ctcom$lambda[, 1:t]) - tload
  tlbias_ctcu <- mean(e_ctcu$lambda[, 1:t]) - tload
  tlbias_ctmin <- mean(e_ctmin$lambda[, 1:t]) - tload
  #mae in trait loadings
  tlmae_ctcmuc <- mean(abs(rowSums(e_ctcmuc$lambda[, 1:t]) - tload))
  tlmae_ctcmr <- mean(abs(rowSums(e_ctcmr$lambda[, 1:t]) - tload))
  tlmae_ctcmpc <- mean(abs(rowSums(e_ctcmpc$lambda[, 1:t]) - tload))
  tlmae_ctcmfc <- mean(abs(rowSums(e_ctcmfc$lambda[, 1:t]) - tload))
  tlmae_ctum <- mean(abs(rowSums(e_ctum$lambda[, 1:t]) - tload))
  tlmae_ctcom <- mean(abs(rowSums(e_ctcom$lambda[, 1:t]) - tload))
  tlmae_ctcu <- mean(abs(rowSums(e_ctcu$lambda[, 1:t]) - tload))
  tlmae_ctmin <- mean(abs(rowSums(e_ctmin$lambda[, 1:t]) - tload))
  # estimated trait correlations
  etvar_ctcmuc <- cov2cor(e_ctcmuc$psi[1:t, 1:t])
  etvar_ctcmr <- cov2cor(e_ctcmr$psi[1:t, 1:t])
  etvar_ctcmpc <- cov2cor(e_ctcmpc$psi[1:t, 1:t])
  etvar_ctcmfc <- cov2cor(e_ctcmfc$psi[1:t, 1:t])
  etvar_ctum <- cov2cor(e_ctum$psi[1:t, 1:t])
  etvar_ctcom <- cov2cor(e_ctcom$psi[1:t, 1:t])
  etvar_ctcu <- cov2cor(e_ctcu$psi[1:t, 1:t])
  etvar_ctmin <- cov2cor(e_ctmin$psi[1:t, 1:t])
  # bias of trait correlations
  tvbias_ctcmuc <- mean(etvar_ctcmuc[lower.tri(etvar_ctcmuc)] - tcor)
  tvbias_ctcmr <- mean(etvar_ctcmr[lower.tri(etvar_ctcmr)] - tcor)
  tvbias_ctcmpc <- mean(etvar_ctcmpc[lower.tri(etvar_ctcmpc)] - tcor)
  tvbias_ctcmfc <- mean(etvar_ctcmfc[lower.tri(etvar_ctcmfc)] - tcor)
  tvbias_ctum <- mean(etvar_ctum[lower.tri(etvar_ctum)] - tcor)
  tvbias_ctcom <- mean(etvar_ctcom[lower.tri(etvar_ctcom)] - tcor)
  tvbias_ctcu <- mean(etvar_ctcu[lower.tri(etvar_ctcu)] - tcor)
  tvbias_ctmin <- mean(etvar_ctmin[lower.tri(etvar_ctmin)] - tcor)
  # mae of trait correlations
  tvmae_ctcmuc <- mean(abs(etvar_ctcmuc[lower.tri(etvar_ctcmuc)] - tcor))
  tvmae_ctcmr <- mean(abs(etvar_ctcmr[lower.tri(etvar_ctcmr)] - tcor))
  tvmae_ctcmpc <- mean(abs(etvar_ctcmpc[lower.tri(etvar_ctcmpc)] - tcor))
  tvmae_ctcmfc <- mean(abs(etvar_ctcmfc[lower.tri(etvar_ctcmfc)] - tcor))
  tvmae_ctum <- mean(abs(etvar_ctum[lower.tri(etvar_ctum)] - tcor))
  tvmae_ctcom <- mean(abs(etvar_ctcom[lower.tri(etvar_ctcom)] - tcor))
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
                   r_pc_equal =r_pc_equal,
                    c_ctcmuc = c_ctcmuc,
                    c_ctcmr = c_ctcmr,
                    c_ctcmpc = c_ctcmpc,
                    c_ctcmfc = c_ctcmfc,
                    c_ctum = c_ctum,
                    c_ctcom = c_ctcom,
                    c_ctcu = c_ctcu,
                    c_ctmin = c_ctmin,
                    chi_ctcmuc = chi_ctcmuc,
                    chi_ctcmr = chi_ctcmr,
                    chi_ctcmpc = chi_ctcmpc,
                    chi_ctcmfc = chi_ctcmfc,
                    chi_ctum = chi_ctum,
                    chi_ctcom = chi_ctcom,
                    chi_ctcu = chi_ctcu,
                    chi_ctmin = chi_ctmin,
                    df_ctcmuc = df_ctcmuc,
                    df_ctcmr = df_ctcmr,
                    df_ctcmpc = df_ctcmpc,
                    df_ctcmfc = df_ctcmfc,
                    df_ctum = df_ctum,
                    df_ctcom = df_ctcom,
                    df_ctcu = df_ctcu,
                    df_ctmin = df_ctmin,
                    cfi_ctcmuc = cfi_ctcmuc,
                    cfi_ctcmr = cfi_ctcmr,
                    cfi_ctcmpc = cfi_ctcmpc,
                    cfi_ctcmfc = cfi_ctcmfc,
                    cfi_ctum = cfi_ctum,
                    cfi_ctcom = cfi_ctcom,
                    cfi_ctcu = cfi_ctcu,
                    cfi_ctmin = cfi_ctmin,
                    tli_ctcmuc = tli_ctcmuc,
                    tli_ctcmr = tli_ctcmr,
                    tli_ctcmpc = tli_ctcmpc,
                    tli_ctcmfc = tli_ctcmfc,
                    tli_ctum = tli_ctum,
                    tli_ctcom = tli_ctcom,
                    tli_ctcu = tli_ctcu,
                    tli_ctmin = tli_ctmin,
                    rmsea_ctcmuc = rmsea_ctcmuc,
                    rmsea_ctcmr = rmsea_ctcmr,
                    rmsea_ctcmpc = rmsea_ctcmpc,
                    rmsea_ctcmfc = rmsea_ctcmfc,
                    rmsea_ctum = rmsea_ctum,
                    rmsea_ctcom = rmsea_ctcom,
                    rmsea_ctcu = rmsea_ctcu,
                    rmsea_ctmin = rmsea_ctmin,
                    in_ctcmuc = in_ctcmuc,
                    in_ctcmr = in_ctcmr,
                    in_ctcmpc = in_ctcmpc,
                    in_ctcmfc = in_ctcmfc,
                    in_ctum = in_ctum,
                    in_ctcom = in_ctcom,
                    in_ctcu = in_ctcu,
                    in_ctmin = in_ctmin,
                    bnd_ctcmuc = bnd_ctcmuc,
                    bnd_ctcmr = bnd_ctcmr,
                    bnd_ctcmpc = bnd_ctcmpc,
                    bnd_ctcmfc = bnd_ctcmfc,
                    bnd_ctum = bnd_ctum,
                    bnd_ctcom = bnd_ctcom,
                    bnd_ctcu = bnd_ctcu,
                    bnd_ctmin = bnd_ctmin,
                    tlbias_ctcmuc = tlbias_ctcmuc,
                    tlbias_ctcmr = tlbias_ctcmr,
                    tlbias_ctcmpc = tlbias_ctcmpc,
                    tlbias_ctcmfc = tlbias_ctcmfc,
                    tlbias_ctum = tlbias_ctum,
                    tlbias_ctcom = tlbias_ctcom,
                    tlbias_ctcu = tlbias_ctcu,
                    tlbias_ctmin = tlbias_ctmin,
                    tlmae_ctcmuc = tlmae_ctcmuc,
                    tlmae_ctcmr = tlmae_ctcmr,
                    tlmae_ctcmpc = tlmae_ctcmpc,
                    tlmae_ctcmfc = tlmae_ctcmfc,
                    tlmae_ctum = tlmae_ctum,
                    tlmae_ctcom = tlmae_ctcom,
                    tlmae_ctcu = tlmae_ctcu,
                    tlmae_ctmin = tlmae_ctmin,
                    tvbias_ctcmuc = tvbias_ctcmuc,
                    tvbias_ctcmr = tvbias_ctcmr,
                    tvbias_ctcmpc = tvbias_ctcmpc,
                    tvbias_ctcmfc = tvbias_ctcmfc,
                    tvbias_ctum = tvbias_ctum,
                    tvbias_ctcom = tvbias_ctcom,
                    tvbias_ctcu = tvbias_ctcu,
                    tvbias_ctmin = tvbias_ctmin,
                    tvmae_ctcmuc = tvmae_ctcmuc,
                    tvmae_ctcmr = tvmae_ctcmr,
                    tvmae_ctcmpc = tvmae_ctcmpc,
                    tvmae_ctcmfc = tvmae_ctcmfc,
                    tvmae_ctum = tvmae_ctum,
                    tvmae_ctcom = tvmae_ctcom,
                    tvmae_ctcu = tvmae_ctcu,
                    tvmae_ctmin = tvbias_ctmin
  )
   return(out)
}