analyze_dat <- function(conditions, condition_number, rep_set, rep, data) {
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
  # Estimate Models
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
  f_ctcmuc <- ifelse(c_ctcmuc, lavInspect(o_ctcmuc, "fit"), NA)
  f_ctcmr <- ifelse(c_ctcmr, lavInspect(o_ctcmr, "fit"), NA)
  f_ctcmpc <- ifelse(c_ctcmpc, lavInspect(o_ctcmpc, "fit"), NA)
  f_ctcmfc <- ifelse(c_ctcmfc, lavInspect(o_ctcmfc, "fit"), NA)
  f_ctcu <- ifelse(c_ctcu, lavInspect(o_ctcu, "fit"), NA)
  f_ctmin <- ifelse(c_ctmin, lavInspect(o_ctmin, "fit"), NA)
  # chisquare
  chi_ctcmuc <- f_ctcmuc['chisq']
  chi_ctcmr <- f_ctcmr['chisq']
  chi_ctcmpc <- f_ctcmpc['chisq']
  chi_ctcmfc <- f_ctcmfc['chisq']
  chi_ctcu <- f_ctcu['chisq']
  chi_ctmin <- f_ctmin['chisq']
  
  
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
                    chi_ctmin = chi_ctmin
  )
  
   return(out)
}