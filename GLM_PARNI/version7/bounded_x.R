bounded_x <- function(x, lowbound = 0.1, upbound = 0.9) {
  return(min(upbound, max(lowbound, x)))
}
