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
  # ctcm: correlated traits correlated methods 
  mdl <- vector("character")
  if (m == 3) {
    if (t == 3) {
      m_ctcm <- 'T1 =~ NA * t1m1 + t1m2 + t1m3
                 T2 =~ NA * t2m1 + t2m2 + t2m3
                 T3 =~ NA * t3m1 + t3m2 + t3m3
                 M1 =~ NA * t1m1 + t2m1 + t3m1
                 M2 =~ NA * t1m2 + t2m2 + t3m2
                 M3 =~ NA * t1m3 + t2m3 + t3m3
                 T1 ~~ 1 * T1; T2 ~~ 1 * T2; T3 ~~ 1 * T3
                 M1 ~~ 1 * M1; M2 ~~ 1 * M2; M3 ~~ 1 * M3 
                 T1 ~~ 0 * M1; T1 ~~ 0 * M2; T1 ~~ 0 * M3 
                 T2 ~~ 0 * M1; T2 ~~ 0 * M2; T2 ~~ 0 * M3 
                 T3 ~~ 0 * M1; T3 ~~ 0 * M2; T3 ~~ 0 * M3'
      m_ctcm_r <- paste(m_ctcm,
                        't1m1 ~~ c11 * t1m1; t1m2 ~~ c12 * t1m2; t1m3 ~~ c13 * t1m3
                        t2m1 ~~ c21 * t2m1; t2m2 ~~ c22 * t2m2; t2m3 ~~ c23 * t2m3
                        t3m1 ~~ c31 * t3m1; t3m2 ~~ c32 * t3m2; t3m3 ~~ c33 * t3m3
                        c11 > 0; c12 > 0; c13 > 0; c21 > 0; c22 > 0; c23 >0
                        c31 > 0; c32 > 0; c33 >0;')
      m_ctcm_c <- paste(m_ctcm,
                        'T1 ~~ a12 * T2; T1 ~~ a13 * T3; T2 ~~ a23 * T3; 
                        a12 > -1; a12 < 1; a13 > -1; a13 < 1; a23 > -1; a23 < 1;'
                        
                        )
    }
  }

  
  c_score <-    'T1 =~ NA * t1m1 + t1m2 + t1m3
                 T2 =~ NA * t2m1 + t2m2 + t2m3
                 T3 =~ NA * t3m1 + t3m2 + t3m3
                 M1 =~ NA * t1m1 + t2m1 + t3m1
                 M2 =~ NA * t1m2 + t2m2 + t3m2
                 t1m1 ~~ c11 * t1m1; t1m2 ~~ c12 * t1m2; t1m3 ~~ c13 * t1m3
                 t2m1 ~~ c21 * t2m1; t2m2 ~~ c22 * t2m2; t2m3 ~~ c23 * t2m3
                 t3m1 ~~ c31 * t3m1; t3m2 ~~ c32 * t3m2; t3m3 ~~ c33 * t3m3
                 c11 > 0; c12 > 0; c13 > 0; c21 > 0; c22 > 0; c23 >0
                 c31 > 0; c32 > 0; c33 >0;
                 T1 ~~ 1 * T1; T2 ~~ 1 * T2; T3 ~~ 1 * T3
                 M1 ~~ 1 * M1; M2 ~~ 1 * M2; 
                 T1 ~~ 0 * M1; T1 ~~ 0 * M2; 
                 T2 ~~ 0 * M1; T2 ~~ 0 * M2; 
                 T3 ~~ 0 * M1; T3 ~~ 0 * M2;
                 M3 =~ NA * t1m3 + t2m3 + t3m3; M3 ~~ 1 * M3; T1 ~~ 0 * M3; T2 ~~ 0 * M3; T3 ~~ 0 * M3;
                 M1 ~~ b12 * M2; M1 ~~ b13 * M3; M2 ~~ b23 * M3;
                 b12 > -1; b12 < 1; b13 > -1; b13 < 1; b23 > -1; b23 < 1;'
  common <- '    T1 ~~ 1 * T1; T2 ~~ 1 * T2; T3 ~~ 1 * T3
                 M1 ~~ 1 * M1; M2 ~~ 1 * M2; 
                 T1 ~~ 0 * M1; T1 ~~ 0 * M2; 
                 T2 ~~ 0 * M1; T2 ~~ 0 * M2; 
                 T3 ~~ 0 * M1; T3 ~~ 0 * M2; '
  # c_method_add <- 'M3 =~ NA * t1m3 + t2m3 + t3m3; M3 ~~ 1 * M3; T1 ~~ 0 * M3; T2 ~~ 0 * M3; T3 ~~ 0 * M3;'
  # um_add <- 'M1 ~~ 0 * M2; M1 ~~ 0 * M3; M2 ~~ 0 * M3;'
  st_add <- 'T1 ~~ 1 * T2; T1 ~~ 1 * T3; T2 ~~ 1 * T3;'
  const_add <- 'T1 ~~ 1 * T2;'
  ct_add <- 'T1 ~~ a13 * T3; T2 ~~ a23 * T3;
            a13 > -1; a13 < 1; a23 > -1; a23 < 1;'
  unconst_add <- 'T1 ~~ a12 * T2; a12 > -1; a12 < 1;'
  # cm_add <- 'M1 ~~ b12 * M2; M1 ~~ b13 * M3; M2 ~~ b23 * M3;
  #            b12 > -1; b12 < 1; b13 > -1; b13 < 1; b23 > -1; b23 < 1;'
  
  i_score  <- 'T1 =~ NA * t1m1_1 + t1m1_2 + t1m1_3 + t1m2_1 + t1m2_2 + t1m2_3 + t1m3_1 + t1m3_2 + t1m3_3
               T2 =~ NA * t2m1_1 + t2m1_2 + t2m1_3 + t2m2_1 + t2m2_2 + t2m2_3 + t2m3_1 + t2m3_2 + t2m3_3
               T3 =~ NA * t3m1_1 + t3m1_2 + t3m1_3 + t3m2_1 + t3m2_2 + t3m2_3 + t3m3_1 + t3m3_2 + t3m3_3
               M1 =~ NA * t1m1_1 + t1m1_2 + t1m1_3 + t2m1_1 + t2m1_2 + t2m1_3 + t3m1_1 + t3m1_2 + t3m1_3
               M2 =~ NA * t1m2_1 + t1m2_2 + t1m2_3 + t2m2_1 + t2m2_2 + t2m2_3 + t3m2_1 + t3m2_2 + t3m2_3
               t1m1_1 ~~ c111 * t1m1_1; t1m1_2 ~~ c112 * t1m1_2; t1m1_3 ~~ c113 * t1m1_3;
               t1m2_1 ~~ c121 * t1m2_1; t1m2_2 ~~ c122 * t1m2_2; t1m2_3 ~~ c123 * t1m2_3;  
               t1m3_1 ~~ c131 * t1m3_1; t1m3_2 ~~ c132 * t1m3_2; t1m3_3 ~~ c133 * t1m3_3;
               t2m1_1 ~~ c211 * t2m1_1; t2m1_2 ~~ c212 * t2m1_2; t2m1_3 ~~ c213 * t2m1_3;    
               t2m2_1 ~~ c221 * t2m2_1; t2m2_2 ~~ c222 * t2m2_2; t2m2_3 ~~ c223 * t2m2_3;  
               t2m3_1 ~~ c231 * t2m3_1; t2m3_2 ~~ c232 * t2m3_2; t2m3_3 ~~ c233 * t2m3_3;
               t3m1_1 ~~ c311 * t3m1_1; t3m1_2 ~~ c312 * t3m1_2; t3m1_3 ~~ c313 * t3m1_3;  
               t3m2_1 ~~ c321 * t3m2_1; t3m2_2 ~~ c322 * t3m2_2; t3m2_3 ~~ c323 * t3m2_3;  
               t3m3_1 ~~ c331 * t3m3_1; t3m3_2 ~~ c332 * t3m3_2; t3m3_3 ~~ c333 * t3m3_3;
               c111 > 0; c112 > 0; c113 > 0;  c121 > 0; c122 > 0; c123 > 0; 
               c131 > 0; c132 > 0; c133 > 0;  c211 > 0; c212 > 0; c213 > 0; 
               c221 > 0; c222 > 0; c223 > 0;  c231 > 0; c232 > 0; c233 > 0;
               c311 > 0; c312 > 0; c313 > 0;  c321 > 0; c322 > 0; c323 > 0;  
               c331 > 0; c332 > 0; c333 > 0;
               T1 ~~ 1 * T1; T2 ~~ 1 * T2; T3 ~~ 1 * T3
               M1 ~~ 1 * M1; M2 ~~ 1 * M2; 
               T1 ~~ 0 * M1; T1 ~~ 0 * M2; 
               T2 ~~ 0 * M1; T2 ~~ 0 * M2; 
               T3 ~~ 0 * M1; T3 ~~ 0 * M2;
               M3 =~ NA * t1m3_1 + t1m3_2 + t1m3_3 + t2m3_1 + t2m3_2 + t2m3_3 + t3m3_1 + t3m3_2 + t3m3_3;
               M3 ~~ 1 * M3; T1 ~~ 0 * M3; T2 ~~ 0 * M3; T3 ~~ 0 * M3;
               M1 ~~ b12 * M2; M1 ~~ b13 * M3; M2 ~~ b23 * M3;
               b12 > -1; b12 < 1; b13 > -1; b13 < 1; b23 > -1; b23 < 1;'
  # i_method_add <- 'M3 =~ NA * t1m3_1 + t1m3_2 + t1m3_3 + t2m3_1 + t2m3_2 + t2m3_3 + t3m3_1 + t3m3_2 + t3m3_3;
  #                  M3 ~~ 1 * M3; T1 ~~ 0 * M3; T2 ~~ 0 * M3; T3 ~~ 0 * M3;'
  model_ctcm <- paste0(c_score, ct_add, unconst_add)
  model_stcm <- paste0(c_score, st_add)
  model_ctcm_constrained <- paste0(c_score, ct_add, const_add)
  # model_ctum <- paste0(c_score, ct_add, unconst_add, c_method_add, um_add)
  # model_stum <- paste0(c_score, st_add, c_method_add, um_add)
  # model_ctum_constrained <- paste0(c_score, ct_add, const_add, c_method_add, um_add)
  
  
  model_i_ctcm <- paste0(i_score, ct_add, unconst_add )
  model_i_stcm <- paste0(i_score, st_add )
  model_i_ctcm_constrained <- paste0(i_score, ct_add, const_add )
  # model_i_ctum <- paste0(i_score, ct_add, unconst_add, i_method_add, um_add)
  # model_i_stum <- paste0(i_score, st_add, i_method_add, um_add)
  # model_i_ctum_constrained <- paste0(i_score, ct_add, const_add, i_method_add, um_add)
  
  fit_ctcm <- cfa(model = model_ctcm, data = cdata) 
  BG <- ifelse(lavInspect(fit_ctcm, "fit")["pvalue"] < .05, 1, 0) #ctcm_p
  # if (is_noproblem(fit_ctcm)) {
  #   ctcm_p <- ifelse(lavInspect(fit_ctcm, "fit")["pvalue"] < .05, 1, 0)
  # }  else {
  #   ctcm_p <- NA
  # }
  fit_stcm <- cfa(model = model_stcm, data = cdata) 
  ctcm_stcm <- compareFit(fit_ctcm, fit_stcm, nested = TRUE)
  W_chisq <- ifelse(ctcm_stcm@nested$`Pr(>Chisq)`[2] > .05, 1, 0)
  W_nfi <- ifelse(ctcm_stcm@fit.diff$nfi < .01, 1, 0)
  W_tli <- ifelse(ctcm_stcm@fit.diff$tli < .01, 1, 0)
  # if (is_noproblem(fit_ctcm) & is_noproblem(fit_stcm)) {
  #   ctcm_stcm <- compareFit(fit_ctcm, fit_stcm, nested = TRUE)
  #   W_chisq <- ifelse(ctcm_stcm@nested$`Pr(>Chisq)`[2] > .05, 1, 0)
  #   W_nfi <- ifelse(ctcm_stcm@fit.diff$nfi < .01, 1, 0)
  #   W_tli <- ifelse(ctcm_stcm@fit.diff$tli < .01, 1, 0)
  # } else {
  #   W_chisq <- NA
  #   W_nfi <- NA
  #   W_tli <- NA
  # }
  ###############################################################################
  # CICFA, Chi-square technique using composite scores
  ###############################################################################
  cutoff <- 1
  ctcm_param <- parameterEstimates(fit_ctcm)
  maxcor <- find_maxcor(fit_ctcm)
  ctcm_upper <- ctcm_param$ci.upper[ctcm_param$label == maxcor]
  ctcm_lower <- ctcm_param$ci.lower[ctcm_param$label == maxcor]
  ctcm_pe <- ctcm_param$est[ctcm_param$label == maxcor]
  ctcm_cicfa <- ifelse(ctcm_upper > cutoff | ctcm_lower < -cutoff |ctcm_pe > .999 | ctcm_pe < -.999, 1, 0)
  ctcm_tli <- lavInspect(fit_ctcm, "fit")["tli"]
  
  
  fit_ctcm_const <- cfa(model = model_ctcm_constrained, data = cdata)
  ctcm_const_dif <- compareFit(fit_ctcm, fit_ctcm_const, nested = TRUE)
  ctcm_chisqdif <- ifelse(ctcm_const_dif@nested$`Pr(>Chisq)`[2] > .05, 1, 0)
  # if(is_noproblem(fit_ctcm) & is_noproblem(fit_ctcm_const)) {
  #   ctcm_const_dif <- compareFit(fit_ctcm, fit_ctcm_const, nested = TRUE)
  #   ctcm_chisqdif <- ifelse(ctcm_const_dif@nested$`Pr(>Chisq)`[2] > .05, 1, 0)
  # } else {
  #   ctcm_chisqdif <- NA
  # }
  ###############################################################################
  # If there is NA in the result of the ctcm model, the ctum model is estimated.
  ###############################################################################
  # if (is.na(ctcm_cicfa) | is.na(ctcm_chisqdif)) {
  #   fit_ctum <- cfa(model_ctum, cdata)
  #   fit_ctum_const <- cfa(model_ctum_constrained, cdata)
  #   fit_stum <- cfa(model_stum, cdata)
  #   ctum_p <- ifelse(lavInspect(fit_ctum, "fit")["pvalue"] < .05, 1, 0)
  #   # if (is.na(ctcm_p) & is_noproblem(fit_ctum)) {
  #   #   ctum_p <- ifelse(lavInspect(fit_ctum, "fit")["pvalue"] < .05, 1, 0)
  #   # } else {
  #   #   ctum_p <- NA
  #   # }
  #   
  #   if (is.na(W_chisq) & is_noproblem(fit_ctum) & is_noproblem(fit_stum)) {
  #     ctum_stum <- compareFit(fit_ctum, fit_stum, nested = TRUE)
  #     ctum_stum_p <- ifelse(ctum_stum@nested$`Pr(>Chisq)`[2] > .05, 1, 0)
  #     ctum_stum_nfi <- ifelse(ctum_stum@fit.diff$nfi < .01, 1, 0)
  #     ctum_stum_tli <- ifelse(ctum_stum@fit.diff$tli < .01, 1, 0)
  #   } else {
  #     ctum_stum_p <- NA
  #     ctum_stum_nfi <- NA
  #     ctum_stum_tli <- NA
  #   }
  #   
  #   if (is.na(ctcm_cicfa) & is_noproblem(fit_ctum)) {
  #     ctum_param <- parameterEstimates(fit_ctum)
  #     ctum_upper <- ctum_param$ci.upper[ctum_param$label == 'a12']
  #     ctum_cicfa <- ifelse(ctum_upper > cutoff, 1, 0)
  #   } else {
  #     ctum_cicfa <- NA
  #   }
  #   
  #   if (is.na(ctcm_chisqdif) & is_noproblem(fit_ctum)) {
  #     ctum_const_dif <- compareFit(fit_ctum, fit_ctum_const, nested = TRUE)
  #     ctum_chisqdif <- ifelse(ctum_const_dif@nested$`Pr(>Chisq)`[2] > .05, 1, 0)
  #   } else {
  #     ctum_chisqdif <- NA
  #   }
  # } else {
  #   ctum_p <- NA
  #   ctum_stum_p <- NA
  #   ctum_stum_nfi <- NA
  #   ctum_stum_tli <- NA
  #   ctum_cicfa <- NA
  #   ctum_chisqdif <- NA
  # }
  
  ###############################################################################
  # CICFA, Chi-square technique using item scores
  ###############################################################################
  fit_i_ctcm <- cfa(model = model_i_ctcm, data = data) 
  fit_i_stcm <- cfa(model = model_i_stcm, data = data) 
  i_ctcm_param <- parameterEstimates(fit_i_ctcm)
  maxcor <- find_maxcor(fit_i_ctcm)
  i_ctcm_upper <- i_ctcm_param$ci.upper[i_ctcm_param$label == maxcor]
  i_ctcm_lower <- i_ctcm_param$ci.lower[i_ctcm_param$label == maxcor]
  i_ctcm_pe <- i_ctcm_param$est[i_ctcm_param$label == maxcor]
  i_ctcm_cicfa <- ifelse(i_ctcm_upper > cutoff |i_ctcm_lower < -cutoff | i_ctcm_pe > .999 | i_ctcm_pe < -.999 , 1, 0)
  # i_ctcm_tli <- lavInspect(fit_ctcm, "fit")["tli"]
  # if (is_noproblem(fit_i_ctcm)) {
  #   i_ctcm_param <- parameterEstimates(fit_i_ctcm)
  #   i_ctcm_upper <- i_ctcm_param$ci.upper[i_ctcm_param$label == 'a12']
  #   i_ctcm_lower <- i_ctcm_param$ci.lower[i_ctcm_param$label == 'a12']
  #   i_ctcm_pe <- i_ctcm_param$est[i_ctcm_param$label == 'a12']
  #   i_ctcm_cicfa <- ifelse(i_ctcm_upper > cutoff |i_ctcm_lower < -cutoff | i_ctcm_pe > .999 | i_ctcm_pe < -.999 , 1, 0)
  #   i_ctcm_tli <- lavInspect(fit_ctcm, "fit")["tli"]
  # } else {
  #   i_ctcm_cicfa <- NA
  #   i_ctcm_pe <- NA
  #   i_ctcm_tli <- NA
  # }
  
  fit_i_ctcm_const <- cfa(model = model_i_ctcm_constrained, data = data)
  i_ctcm_const_dif <- compareFit(fit_i_ctcm, fit_i_ctcm_const, nested = TRUE)
  i_ctcm_chisqdif <- ifelse(i_ctcm_const_dif@nested$`Pr(>Chisq)`[2] > .05, 1, 0)
  # if(is_noproblem(fit_i_ctcm) & is_noproblem(fit_i_ctcm_const)) {
  #   i_ctcm_const_dif <- compareFit(fit_i_ctcm, fit_i_ctcm_const, nested = TRUE)
  #   i_ctcm_chisqdif <- ifelse(i_ctcm_const_dif@nested$`Pr(>Chisq)`[2] > .05, 1, 0)
  # } else {
  #   i_ctcm_chisqdif <- NA
  # }
  
  # if (is.na(i_ctcm_cicfa) | is.na(i_ctcm_chisqdif)) {
  #   fit_i_ctum <- cfa(model_i_ctum, data)
  #   fit_i_ctum_const <- cfa(model_i_ctum_constrained, data)
  #   fit_i_stum <- cfa(model_i_stum, data)
  # 
  #       if (is.na(i_ctcm_cicfa) & is_noproblem(fit_i_ctum)) {
  #         i_ctum_param <- parameterEstimates(fit_i_ctum)
  #         i_ctum_upper <- i_ctum_param$ci.upper[i_ctum_param$label == 'a12']
  #         i_ctum_cicfa <- ifelse(i_ctum_upper > cutoff, 1, 0)
  #   } else {
  #     i_ctum_cicfa <- NA
  #     
  #   }
  #   
  #   if (is.na(i_ctcm_chisqdif) & is_noproblem(fit_i_ctum)) {
  #     i_ctum_const_dif <- compareFit(fit_i_ctum, fit_i_ctum_const, nested = TRUE)
  #     i_ctum_chisqdif <- ifelse(i_ctum_const_dif@nested$`Pr(>Chisq)`[2] > .05, 1, 0)
  #   } else {
  #     i_ctum_chisqdif <- NA
  #   }
  # } else {
  #   i_ctum_cicfa <- NA
  #   i_ctum_chisqdif <- NA
  # }
  # 
  
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