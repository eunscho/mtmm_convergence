modelcm <- function(model, t, m) {
  for (j in 1:m) {
    model <- paste0(model, 'M', j, ' ~~ 1 * M', j, ';')
  }
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      model <- paste0(model, 'M', i, ' ~~ b', i, j, ' * M', j, ';')
    }
  }
  for (i in 1:t) {
    for (j in 1:m) {
      model <- paste0(model, 'T', i, ' ~~ 0 * M', j, ';')
    }
  }
  return(model)
}