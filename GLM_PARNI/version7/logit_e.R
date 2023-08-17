# (inverse) logistic function with epsilon

logit_e <- function(x, eps){
  
  x[x > 2*(1-eps)] <- 1-2*eps
  x[x < 2*eps] <- 2*eps
  
  return(log(x - eps) - log(1 - x - eps))
}



inv_logit_e <- function(y, eps){
  ey <- exp(-y)
  return((eps*ey - eps + 1)/(ey + 1))
}