mvnorm_trans <- function(v, mu, sigma) {
  
  LS <- chol(sigma)
  p <- nrow(v)
  density <- - p*log(2*pi)/2 - sum(log(diag(LS))) - colSums(v^2)/2
  
  v <- as.numeric(mu) + t(LS) %*% v
  
  return(list(v = v,
              density = density))
  
}
