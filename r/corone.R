corone <- function(cor, cutoff) {
  diag(cor) <- 0
  cor <- abs(cor)
  if (max(cor) > 1 - cutoff) {
    out <- TRUE
  } else {
    out <- FALSE
  }
  return(out)
}