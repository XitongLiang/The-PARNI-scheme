

generate_temperature <- function(taus) {
  cumprod(c(1, exp(-exp(taus))))
}