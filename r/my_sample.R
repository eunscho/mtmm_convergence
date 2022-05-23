my_sample <- function(base, rep) {
  patterns <- c(1,2,3)
  if (rep %% 3 == 0) {
    mod <- sample(patterns, rep / 3)
  }
  return(mod)
}

