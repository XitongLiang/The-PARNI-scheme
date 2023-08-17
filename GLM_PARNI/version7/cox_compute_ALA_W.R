cox_riskset <- function(time_of_interest, t) {
  return(t >= time_of_interest)
}


cox_compute_ALA_W <- function(t, d) {
  
  M <- sapply(t, cox_riskset, t)
  n <- length(t)
  
  exp_eta <- rep(1, n)
  
  order_t <- order(t, d, decreasing = c(FALSE, TRUE))
  
  t <- t[order_t]
  d <- d[order_t]
  
  t2 <- rev(cumsum(rev(exp_eta)))
  
  cumsum_d <- cumsum(d)
  t2_inv <- cumsum((1/t2)[d == 1])
  
  t2_inv <- c(0, t2_inv)
  
  Q <- t2_inv[cumsum_d + 1]
  
  W <- Q
  W[order_t] <- W
  W <- diag(W)
  
  d[order_t] <- d
  t2[order_t] <- t2
  
  X <- t(t(M)*d/(t2)^2)
  Y <- t(M) 
  
  # print(X)
  
  # W <- W - (exp_eta %o% exp_eta) * Q2
  W <- W - X %*% Y
  
  
  return(W)
  
}




