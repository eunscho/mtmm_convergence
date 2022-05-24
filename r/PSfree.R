PSfree <- function(trait = 3, method = 3) {
  len <- (trait + method)
  out <- matrix(0, len, len)
  for (i in 1:trait) { # Trait factor correlation
    for (j in 1:i) {
      if (i == j) {
        out[i, j] <- 1 # variance to 1
      } else if (i > j){
        out[i, j] <- paste0("tr", i, j)
      } 
    }
  }
  for (i in 1:(trait - 1)) {
    for (j in (i + 1): trait) {
      out[i, j] <- out[j, i]
    }
  }
  
  
  for (i in (trait + 1):len) {
    for (j in (trait + 1):i) {
      if (i == j) {
        out[i, j] <- 1 # variance to 1
      } else if (i > j){
        out[i, j] <- paste0("mr", i, j)
      } 
    }
  }
  for (i in (trait + 1):(len - 1)) {
    for (j in (i + 1): len) {
      out[i, j] <- out[j, i]
    }
  }
  
  return(out)
}