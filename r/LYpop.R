LYpop <- function(trait = 3, method = 3, item = 3, trait_loading, method_loading) {
  out <- matrix(0, trait * method * item, trait + method)
  for (i in 1:(trait * method * item)) {
    for (j in 1:(trait + method)) {
      tnum <- mydiv(ceiling(i / item), trait)
      if (j <= trait) { # trait loadings
        if (tnum == j) {
          out[i, j] <- trait_loading[i]
        } 
      } else { # method loadings
        mnum <- ceiling(i / (item * trait))
        if (mnum == j - trait) {
          out[i, j] <- method_loading[i]
        }
      } # end of } else { # method loadings
    } # end of for (j in 1:(trait + method)) {
  } # end of for (i in 1:(trait * method * item)) {
  return(out)
}

