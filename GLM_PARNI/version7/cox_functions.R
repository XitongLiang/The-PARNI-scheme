cox_compute_t2 <- function(exp_eta, t) {
  
  order_t <- order(t, decreasing = TRUE)
  
  exp_eta <- exp_eta[order_t]
  t <- t[order_t]
  
  t2 <- cumsum(exp_eta)
  
  t2[order_t] <- t2
  
  return(t2)
  
}




cox_pseudores <- function (exp_eta, t, d) {
  
  # eta <- X %*% beta
  # exp_eta <- exp(eta)
  # p <- ncol(X)
  # n <- nrow(X)
  
  order_t <- order(t, d, decreasing = c(FALSE, TRUE))
  
  exp_eta <- exp_eta[order_t]
  t <- t[order_t]
  d <- d[order_t]
  
  t2 <- rev(cumsum(rev(exp_eta)))
  
  cumsum_d <- cumsum(d)
  t2_inv <- cumsum((1/t2)[d == 1])
  
  t2_inv <- c(0, t2_inv)
  
  Q <- t2_inv[cumsum_d + 1]
  
  
  # Q <- sapply(t, cox_sum_D_j, t, exp_eta, d, t2)
  
  
  tilde_y <- d -  (exp_eta * Q)
  # tilde_y <- t(X_v) %*% d  -  t(X_v) %*% (exp_eta * Q)
  
  tilde_y[order_t] <- tilde_y
  
  return(tilde_y)
  
}



