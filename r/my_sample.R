my_sample <- function(base, rep) {
  out <- double(rep)
  qt <- floor (rep / 3) #quotient
  rmd <- rep %% 3 # remainder
  pattern <- sample(c(1, 2, 3, 4, 5, 6), qt, replace = TRUE)
  for (i in 1:length(pattern)) {
    if (pattern[i] == 1) {
      out[(3 * (i - 1) + 1): (3 * (i - 1) + 3)] <- c(0, -.1, .1)
    } else if (pattern[i] == 2) {
      out[(3 * (i - 1) + 1): (3 * (i - 1) + 3)] <- c(0, .1, -.1)
    } else if (pattern[i] == 3) {
      out[(3 * (i - 1) + 1): (3 * (i - 1) + 3)] <- c(.1, -.1, 0)
    } else if (pattern[i] == 4) {
      out[(3 * (i - 1) + 1): (3 * (i - 1) + 3)] <- c(.1, 0, -.1)
    } else if (pattern[i] == 5) {
      out[(3 * (i - 1) + 1): (3 * (i - 1) + 3)] <- c(-.1, 0, .1)
    } else {
      out[(3 * (i - 1) + 1): (3 * (i - 1) + 3)] <- c(-.1, .1, 0)
    }
  }
  if (rmd == 1) {
    out[rep] <- 0
  } else if (rmd == 2) {
    out[(rep - 1): rep ] <- c(.1, -.1)
  }
  out <- base + out
  return(out)
}

