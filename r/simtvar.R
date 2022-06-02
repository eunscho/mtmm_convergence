simtvar <- function(tloadmat, mloadmat = 0, tcov, mcov = 0, errorcov, ctcu = FALSE) {
  tvar <- tloadmat %*% tcov %*% t(tloadmat)
  if (ctcu) {
    out <- sum(tvar)/(sum(tvar) + sum(errorcov))
  } else {
    mvar <- mloadmat %*% mcov %*% t(mloadmat)
    out <- sum(tvar)/(sum(tvar) + sum(mvar) + sum(errorcov))
  }
 return(out)
}