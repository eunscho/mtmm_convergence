LYfree <- function(trait = 3, method = 3, item = 3) {
  out <- matrix(0, trait * method * item, trait + method)
  for (i in 1:(trait * method * item)) {
    for (j in 1:(trait + method)) {
      tnum <- mydiv(ceiling(i / item), trait)
      mnum <- ceiling(i / (item * trait))
      inum <- mydiv(i, item)
      if (j <= trait) { # trait loadings
        if (tnum == j) {
          out[i, j] <- paste0("t", tnum, "_", "t", tnum, "m", mnum, "-", inum)
        } 
      } else { # method loadings
        if (mnum == j - trait) {
          out[i, j] <- paste0("m", mnum, "_", "t", tnum, "m", mnum, "-", inum)
        }
      } # end of } else { # method loadings
    } # end of for (j in 1:(trait + method)) {
  } # end of for (i in 1:(trait * method * item)) {
  return(out)
}

