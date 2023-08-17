# sampling half-Cauchy random variables
# with number of samplers n
# location parameter mu
# scale parameter sigma2
rhcauchy <- function(n, mu, sigma2) {
  return(tan(runif(n)*pi/2))
}