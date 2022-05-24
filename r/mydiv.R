mydiv <- function(num1, num2) {
  if (num1 %% num2 == 0) {
    out <- num2
  } else {
    out <- num1 %% num2
  }
  return (out)
}