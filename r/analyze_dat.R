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

  cdata <- composite(data, t, m, k) # transform items scores to composite scores
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
  # Report output
  ###############################################################################
  out <- data.frame(condition_number = condition_number,
                    rep_set = rep_set,
                    rep = rep,
                    n = n,
                    fcor = fcor,
                    nfcor = nfcor,
                    mcor = mcor,
                    tpat = tpat,
                    mpat = mpat,
                    CF1 = CF1,
                    # raw_hthm_5tol = raw_hthm_5tol,
                    # raw_hthm_10tol = raw_hthm_10tol,
                    CF2 = CF2,
                    # raw_htmm_20tol = raw_htmm_20tol,
                    # raw_htmm_30tol = raw_htmm_30tol,
                    CF3 = CF3,
                    CF4 = CF4,
                    GT1 = GT1,
                    GT2 = GT2,
                    DC1 = DC1,
                    # cor_hthm_5tol = cor_hthm_5tol,
                    # cor_hthm_10tol = cor_hthm_10tol,
                    DC2 = DC2,
                    # cor_htmm_20tol = cor_htmm_20tol,
                    # cor_htmm_30tol = cor_htmm_30tol,
                    DC3 = DC3,
                    # cor_gt1 = cor_gt1,
                    # cor_gt2 = cor_gt2,
                    # ctcm_pe <- ctcm_pe,
                    # ctcm_tli <- ctcm_tli,
                    BG = BG,
                    # ctum_p = ctum_p,
                    W_chisq = W_chisq,
                    # ctum_stum_p = ctum_stum_p,
                    W_nfi = W_nfi,
                    # ctum_stum_nfi = ctum_stum_nfi,
                    W_tli = W_tli,
                    # ctum_stum_tli = ctum_stum_tli,
                    ctcm_cicfa = ctcm_cicfa,
                    # ctum_cicfa = ctum_cicfa,
                    ctcm_chisqdif = ctcm_chisqdif,
                    # ctum_chisqdif = ctum_chisqdif,
                    # i_ctcm_pe = i_ctcm_pe,
                    # i_ctcm_tli = i_ctcm_tli,
                    i_ctcm_cicfa = i_ctcm_cicfa,
                    # i_ctum_cicfa = i_ctum_cicfa,
                    i_ctcm_chisqdif = i_ctcm_chisqdif
                    # i_ctum_chisqdif = i_ctum_chisqdif
                    # i_cicfa = i_cicfa,
                    # i_chisqdif = i_chisqdif
  )
  return(out)
}