outrange <- function(cov, tolerance) { # return if one of correlation is greater than 1 or less than -1
  if (max(cov2cor(abs(cov))) > 1 + tolerance) {
    out <- TRUE
  } else {
    out <- FALSE
  }
return(out)
}