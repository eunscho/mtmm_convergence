main <- function(start = 1, end = 1296) {
  n <- c(100, 250, 500, 1000)
  fcor <- c(.7, 1 - 1e-10) # focal traits correlation (t1  t2)
  nfcor <- c(.3, .5) # nonfocal traits correlation
  mcor <- c(.2, .4) # method correlation
  tpat <- list(c(.5, .5, .5), c(.7, .3, .5)) # trait loading pattern
  mpat <- list(c(.3, .3, .3), c(.4, .2, .3)) # method loading pattern
  conditions <- tidyr::crossing(n, fcor, nfcor, mcor, tpat, mpat)
  colnames(conditions) <- c("n", "fcor", "nfcor", "mcor", "tpat", "mpat")
  REP_PER_CONDITION <- 1000
  SET_PER_CONDITION <- 10
  rep_sets <- 1:SET_PER_CONDITION
  reps_per_set <- 1:(REP_PER_CONDITION/SET_PER_CONDITION)
  
  for (condition_number in start:end) {
    for (rep_set in rep_sets) {
      tictoc::tic()
      print(paste("Starting condition number", condition_number, "rep", rep_set))
      print(conditions[condition_number, ])
      filename <- paste0("allcor", condition_number, "-", rep_set, ".csv")
      if (!file.exists(filename)) {
        for (rep in reps_per_set) {
          cat("allcor: ", condition_number, "rep set: ", rep_set, "rep: ", rep)
          data <- generate_dat (conditions, condition_number, rep_set, rep)
          temp <- analyze_dat (conditions, condition_number, rep_set, rep, data)
          if (rep == 1) {
            out <- temp
          } else {
            out <- rbind(out, temp)
          }
        } # end of for (rep in reps_per_set)
        
        # eliminate duplicates
        out_edit <- out[1:100,]
        for (i in 1:300) {
          if (i %% 3 == 0) {
            out_edit[round(i / 3), ] <- out[i, ]
          }
        }
        
        readr::write_csv(out_edit, file = filename)
        print(out_edit)
        tictoc::toc()
      } # end of if (!file.exists(filename))
    } # end of  for (rep_set in rep_sets)
  } # end of for (condition_number in condition_numbers)
} # end of function