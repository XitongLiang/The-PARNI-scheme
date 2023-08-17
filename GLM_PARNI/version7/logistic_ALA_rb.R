logistic_ALA_rb <- function(ALL, hyper_par) {
  
  p <- hyper_par$p
  p_z <- hyper_par$p_z
  
  
  tilde_y <- 4*(hyper_par$y - 1/2)
  X <- hyper_par$X
  inv_variances <- hyper_par$inv_variances
  h <- hyper_par$h_exp
  
  
  curr <- ALL$curr
  p_gam <- ALL$p_gam
  g <- ALL$g
  
  includes <- which(curr == 1)
  J_gamma <- cbind(hyper_par$Z, X[, includes])
  
  B <- t(J_gamma) %*% J_gamma
  diag(B) <- diag(B) + 4 * c(inv_variances, rep(1/g, p_gam))
  inv_B <- solve(B)
  L_invB <- chol(inv_B)
  
  tyT_Jg_invB <- t(tilde_y) %*% J_gamma %*% inv_B
  XT_Jg_LinvB <- t(X) %*% J_gamma %*% t(L_invB)
  tyT_Jg_LinvB <- as.numeric(t(tilde_y) %*% J_gamma %*% t(L_invB))
  
  tyT_X <- as.numeric(t(tilde_y) %*% X)
  
  diag_XtX <- colSums(X^2)
  
  # Bayes_factor <- rep(NA, p)
  
  d_vec <- diag_XtX + 4/(g) - rowSums(XT_Jg_LinvB^2)
  
  d_vec[includes] <- 1
  Bayes_factor <- d_vec^{-1/2}*(4/g)^{1/2}*exp((tyT_X - colSums(tyT_Jg_LinvB * t(XT_Jg_LinvB)))^2/(8*d_vec))
  
  q_j <- 0
  
  for (j in includes) {
    q_j <- q_j + 1
    d_j <- 1/inv_B[(p_z+q_j),(p_z+q_j)]
    Bayes_factor[j] <- d_j^{-1/2}*(4/g)^{1/2}*exp(d_j*tyT_Jg_invB[(p_z+q_j)]^2/8)
  }
  
  # print(Bayes_factor)
  
  PIPs <- h*Bayes_factor/(1-h+h*Bayes_factor)
  
  
  return(PIPs)
  
}








