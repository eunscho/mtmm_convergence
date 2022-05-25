main <- function(start = 1, end = 1296) {
  n <- c(100, 250, 500, 1000)
  m <- c(3, 5)
  t <- c(3, 5)
  mload <- c(.1, .3, .5) # fan & lance c(.17, .28, .48)
  tload <- c(.3, .5, .7) # fan & lance c(.31, .50, .69)
  mcor <- c(0, .3, .6) # fan & lance c(.02, .29, .58)
  tcor <- c(.1, .4, .7) # fan & lance c(.07, .36, .61)
  conditions <- tidyr::crossing(n, m, t, mload, tload, mcor, tcor)
  colnames(conditions) <- c("n", "m", "t", "mload", "tload", "mcor", "tcor")
  REP_PER_CONDITION <- 100
  SET_PER_CONDITION <- 10
  rep_sets <- 1:SET_PER_CONDITION
  reps_per_set <- 1:(REP_PER_CONDITION/SET_PER_CONDITION)
  
  for (condition_number in start:end) {
    for (rep_set in rep_sets) {
      tictoc::tic()
      print(paste("Starting condition number", condition_number, "rep", rep_set))
      print(conditions[condition_number, ])
      filename <- paste0("conv", condition_number, "-", rep_set, ".csv")
      if (!file.exists(filename)) {
        for (rep in reps_per_set) {
          cat("conv: ", condition_number, "rep set: ", rep_set, "rep: ", rep)
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