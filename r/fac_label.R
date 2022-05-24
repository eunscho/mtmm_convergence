fac_label <- function(trait, method) {
  out <- vector("character", trait + method)
  for (i in 1:length(out)) {
    if (i <= trait) {
      out[i] <- paste0("t", i)
    } else {
      out[i] <- paste0("m", i - trait)
    }
  }
  return(out)
}