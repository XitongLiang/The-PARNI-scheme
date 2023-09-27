log_llh_DAG <- function(LA, hyper_par) {
  
  gamma <- LA$curr
  
  X <- hyper_par$X
  g <- hyper_par$g
  
  XtX <- hyper_par$XtX
  
  p <- nrow(X)
  n <- ncol(X)
  p_gamma <- sum(gamma)
  
  log_sqrt_det_Sigma <- rep(0, p)
  A <- - n*log(diag(XtX)/2)/2
  
  num_parents <- colSums(gamma)
  
  nodes_connected <- which(num_parents > 0)
  # nodes_unconnected <- which(num_parents == 0)
  
  for (j in nodes_connected) {
    
    # print(j)
    
    parents_j <- which(gamma[,j] == 1)
    # x_j <- X[j,]
    p_paj <- length(parents_j)
    
    x_paj_xj_t <- XtX[parents_j, j]
    
    diag_p_paj <- diag(p_paj)
    
    Vg <- XtX[parents_j, parents_j] + diag_p_paj/g
    # diag(Vg) <- diag(Vg) + 1/g
    L_Vg <- chol(Vg)
    
    A[j] <- - n * log((XtX[j,j] - sum((forwardsolve(t(L_Vg), diag_p_paj) %*%
                                         x_paj_xj_t)^2))/2)/2
    # A[j] <- XtX[j,j] - sum((solve(t(L_Vg)) %*% x_paj_xj_t)^2)
    
    if (p_paj == 1){
      log_sqrt_det_Sigma[j] <- log(L_Vg)
    }
    else {
      log_sqrt_det_Sigma[j] <- sum(log(diag(L_Vg)))
    }
    
  }
  
  # log_llh <- - p_gamma*log(g)/2 - sum(log_sqrt_det_Sigma) - n*p*log(sum(A))/2
  log_llh <- - p_gamma*log(g)/2 - sum(log_sqrt_det_Sigma) + sum(A)
  LA$llh <- log_llh
  # LA$A <- A
  # LA$log_det_sigma <- log_sqrt_det_Sigma
  LA$p_gam <- p_gamma
  
  return(LA)
  
}



log_llh_DAG_table <- function(LA, hyper_par) {
  
  gamma <- LA$curr
  
  X <- hyper_par$X
  g <- hyper_par$g
  permi_pars <- hyper_par$permi_pars
  tables <- hyper_par$tables
  
  XtX <- hyper_par$XtX
  
  p <- nrow(X)
  n <- ncol(X)
  p_gamma <- sum(gamma)
  
  log_sqrt_det_Sigma <- rep(0, p)
  A <- - n*log(diag(XtX)/2)/2
  
  num_parents <- colSums(gamma)
  
  nodes_connected <- which(num_parents > 0)
  # nodes_unconnected <- which(num_parents == 0)
  
  
  
  for (j in nodes_connected) {
    
    # print(j)
    
    parents_j <- which(gamma[,j] == 1)
    # x_j <- X[j,]
    p_paj <- length(parents_j)
    
    check <- !parents_j %in% permi_pars[[j]]
    
    
    if (sum(check) == 0) {
      
      m <- model_encoding(gamma[permi_pars[[j]],j])
      
      A[j] <- tables[[j]][[2]][m]
      
    }
    else if (sum(check) == 1) {
      
      
      table_j <- which(tables[[j]][[1]] == parents_j[check])
      
      m <- model_encoding(gamma[permi_pars[[j]],j])
      
      
      A[j] <- tables[[j]][[table_j]][m]
      # print(j)
      
    }
    else {
      
      x_paj_xj_t <- XtX[parents_j, j]
      
      diag_p_paj <- diag(p_paj)
      
      Vg <- XtX[parents_j, parents_j] + diag_p_paj/g
      # diag(Vg) <- diag(Vg) + 1/g
      L_Vg <- chol(Vg)
      
      A[j] <- - p_paj*log(g)/2 - n * log((XtX[j,j] - sum((forwardsolve(t(L_Vg), diag_p_paj) %*%
                                                            x_paj_xj_t)^2))/2)/2
      # A[j] <- XtX[j,j] - sum((solve(t(L_Vg)) %*% x_paj_xj_t)^2)
      
      if (p_paj == 1){
        log_sqrt_det_Sigma[j] <- log(L_Vg)
        
      }
      else {
        log_sqrt_det_Sigma[j] <- sum(log(diag(L_Vg)))
      }
      
    }
    
    # log_post1 <- usrDAGcorescore(j, parents_j, p, hyper_par)
    # log_llh2 <- A[j] - log_sqrt_det_Sigma[j]
    # lmp2 <- hyper_par$log_m_prior(p_paj, hyper_par$h, hyper_par$max_p)
    
    # cat(abs(log_post1 - log_llh2-lmp2), log_post1, log_llh2+lmp2, "\n")
    
    
  }
  # cat(p_temp, p_gamma, "\n")
  
  # log_llh <- - p_gamma*log(g)/2 - sum(log_sqrt_det_Sigma) - n*p*log(sum(A))/2
  log_llh <- - sum(log_sqrt_det_Sigma) + sum(A)
  LA$llh <- log_llh
  LA$A <- A
  LA$log_det_sigma <- log_sqrt_det_Sigma
  LA$p_gam <- p_gamma
  
  
  return(LA)
  
}

