popvar <- function(t, m, tload, mload, tcor, mcor) {
  tloadmat <- matrix(0, t * m, t)
  for (i in 1:t) {
    for (j in 1:(t * m)) {
      if ((j - 1) %/% m == (i - 1)) {
        tloadmat[j, i] <- tload
      }
    }
  }
  for (i in 1:t) {
    if (i == 1) {
      mloadmat <- mloadunit <- diag(m) * mload
    } else {
      mloadmat <- rbind(mloadmat, mloadunit)
    }
  }
  tcormat <- matrix(tcor, t, t)
  diag(tcormat) <- 1
  mcormat <- matrix(mcor, m, m)
  diag(mcormat) <- 1
  
  tvar <- tloadmat %*% tcormat %*% t(tloadmat)
  mvar <- mloadmat %*% mcormat %*% t(mloadmat)
  evar <- diag(diag(1, t * m) - diag(tvar + mvar))
  out <- sum(tvar)/(sum(tvar) + sum(mvar) + sum(evar))
 return(out)
}