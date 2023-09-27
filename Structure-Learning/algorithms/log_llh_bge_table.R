log_llh_bge_table <- function(LA, hyper_par) {
  
  scorepar <- hyper_par$scorepar
  
  gamma <- LA$curr
  p_gamma <- sum(gamma)
  num_parents <- colSums(gamma)
  
  permi_pars <- hyper_par$permi_pars
  tables <- hyper_par$tables
  
  p <- hyper_par$p
  n <- hyper_par$n
  
  A <- rep(0, p)
  
  for (j in 1:p) {
    
    parents_j <- which(gamma[,j] == 1)
    # x_j <- X[j,]
    p_paj <- length(parents_j)
    
    if (p_paj == 0) {
      
      A[[j]] <- tables[[j]][[2]][1]
      
    }
    else {
      
      
      
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
        
        A[j] <- BiDAG:::DAGcorescore(j, parentnodes = parents_j, p, scorepar)
        
      }
      
    }
    
  }
  
  
  log_llh <- sum(A)
  LA$llh <- log_llh
  LA$A <- A
  LA$p_gam <- p_gamma
  
  return(LA)
  
}









log_llh_bge_update_table <- function(changes, LA_old, LA, hyper_par) {
  
  A <- LA_old$A
  
  gamma <- LA$curr
  
  permi_pars <- hyper_par$permi_pars
  tables <- hyper_par$tables
  scorepar <- hyper_par$scorepar
  p <- hyper_par$p
  n <- hyper_par$n
  
  p_gamma <- sum(gamma)
  
  
  # log_sqrt_det_Sigma <- rep(0, p)
  # A <- - n*log(diag(XtX)/2)/2
  
  num_parents <- colSums(gamma)
  
  # nodes_connected <- which(num_parents > 0)
  # nodes_unconnected <- which(num_parents == 0)
  
  
  for (j in changes) {
    
    if (num_parents[j] == 0) {
      
      A[j] <- tables[[j]][[2]][1]
      
    }
    else {
      
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
        
        A[j] <- BiDAG:::DAGcorescore(j, parentnodes = parents_j, p, scorepar)
        
      }
      
      
    }
    
  }
  
  
  log_llh <- sum(A)
  LA$llh <- log_llh
  LA$A <- A
  LA$p_gam <- p_gamma
  
  return(LA)
  
}



