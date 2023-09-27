# independent prior
marPIPs_DAG_H <- function(hyper_par, kappa = 0){
  
  n <- hyper_par$n
  p <- hyper_par$p
  permi_pars <- hyper_par$permi_pars
  
  PIPs <- matrix(0, p, p)
  tables <- list()
  
  
  for (j in 1:p) {
    
    Pa_j <- permi_pars[[j]]
    
    D <- enumerate_full_models_C(j, Pa_j, hyper_par)
    
    # D <- enumerate_full_models(j, Pa_j, hyper_par)
    
    # if (any(is.nan(D$PIPs))) {
    #   D <- enumerate_full_models(j, Pa_j, hyper_par)
    # }
    
    PIPs[,j] <- D$PIPs
    tables[[j]] <- D$table
    
    
    # cat(j ,sum(abs(D$PIPs - enumerate_full_models(j, Pa_j, hyper_par)$PIPs)), "\n")
    
  }
  
  
  
  PIPs <- PIPs * (1-2*kappa) + kappa
  ratios <- PIPs/(1-PIPs)
  PIPs <- ratios/(ratios + t(ratios) + 1)
  # PIPs[which(H) != 1] = 0
  
  
  return(list(tables = tables, 
              PIPs = PIPs))
}




# table[[j]][[1]]
# index


enumerate_full_models_C <- function(j, Pa_j, hyper_par) {
  
  
  
  X <- hyper_par$X
  g <- hyper_par$g
  n <- hyper_par$n
  p <- hyper_par$p
  h <- hyper_par$h
  log_m_prior <- hyper_par$log_m_prior
  max_p <- hyper_par$max_p
  
  inv_g <- 1/g
  inv_sqrt_g <- sqrt(inv_g)
  n_power <- -n/2
  
  PIPs <- rep(0, p)
  index <- c(0,0, setdiff(1:p, c(j, Pa_j)))
  table <- list()
  table[[1]] <- index
  
  x_j <- X[j,]
  xj_t_xj <- hyper_par$XtX[j, j] 
  
  log_llh_null <- - n*log(xj_t_xj/2)/2
  lmp_null <- log_m_prior(0, h, max_p)
  
  prob_curr_null <- exp(log_llh_null + lmp_null)
  
  
  C <- 1
  max_C <- log_llh_null + lmp_null
  
  p_Pa <- length(Pa_j)
  k_power <- 2^p_Pa
  
  for (i in 2:length(index)) {
    
    k <- index[i]
    
    if (k == 0) {
      Pa_j_new <- Pa_j
      p_Pa_new <- length(Pa_j_new)
      if (p_Pa_new == 0) {
        table[[i]] <- log_llh_null
        next
      }
    }
    else {
      Pa_j_new <- c(Pa_j, k)
      p_Pa_new <- length(Pa_j_new)
    }
    
    
    
    x_pa_j <- X[Pa_j_new,]
    
    V <- hyper_par$XtX[Pa_j_new, Pa_j_new]
    x_paj_xj_t <- hyper_par$XtX[Pa_j_new, j]
    
    if (p_Pa_new == 1) {
      V <- as.matrix(V + inv_g)
    }
    else {
      diag(V) <- diag(V) + inv_g
    }
    
    
    log_llhs <- rep(log_llh_null, k_power)
    
    curr <- rep(0, p_Pa_new)
    
    if (p_Pa_new != p_Pa) {
      curr[p_Pa_new] <- 1
    }
    
    k_matrix_old <- rep(0, p_Pa)
    
    
    
    for (idx in 0:(k_power-1)) {
      
      k_matrix <- Gray_code(idx, k_matrix_old)
      temp_j <- which(k_matrix != k_matrix_old)
      curr[temp_j] <- 1 - curr[temp_j]
      k_matrix_old <- k_matrix
      
      curr_incl <- which(curr == 1)
      Pa_j_incl <- Pa_j_new[curr_incl]
      p_curr <- length(curr_incl)
      
      if (p_curr == 0) {
        next
      }
      
      
      x_paj_xj_t_curr <- x_paj_xj_t[curr_incl]
      L_Vg <- chol(V[curr_incl, curr_incl])
      A <- - n * log((xj_t_xj - sum((forwardsolve(t(L_Vg), diag(p_curr)) %*%
                                       x_paj_xj_t_curr)^2))/2)/2
      
      if (p_curr == 1){
        log_sqrt_det_Sigma <- log(L_Vg)
      }
      else {
        log_sqrt_det_Sigma <- sum(log(diag(L_Vg)))
      }
      
      log_llh <- - p_curr*log(g)/2 - sum(log_sqrt_det_Sigma) + A 
      
      lmp <- log_m_prior(p_curr, h, max_p)
      
      prob_curr <- exp(log_llh + lmp - max_C)
      
      
      if (prob_curr > 1) {
        
        PIPs <- PIPs/prob_curr
        C <- C/prob_curr
        
        max_C <- log_llh + lmp # log(prob_curr)
        prob_curr <- exp(log_llh + lmp - max_C)
        
        # print(max_C)
      }
      
      
      
      PIPs[Pa_j_incl] <- PIPs[Pa_j_incl] + prob_curr
      
      if (p_Pa_new != p_Pa) {
        m <- model_encoding(curr[-p_Pa_new])
      }
      else {
        m <- model_encoding(curr)
      }
      
      log_llhs[m] <- log_llh
      
      
      C <- C + prob_curr
      
      
    }
    
    
    table[[i]] <- log_llhs
    
    
  }
  
  return(list(table = table,
              PIPs = PIPs/C))
  
}




enumerate_full_models <- function(j, Pa_j, hyper_par) {
  
  
  
  X <- hyper_par$X
  g <- hyper_par$g
  n <- hyper_par$n
  p <- hyper_par$p
  h <- hyper_par$h
  log_m_prior <- hyper_par$log_m_prior
  max_p <- hyper_par$max_p
  
  inv_g <- 1/g
  inv_sqrt_g <- sqrt(inv_g)
  n_power <- -n/2
  
  PIPs <- rep(0, p)
  index <- c(0,0, setdiff(1:p, c(j, Pa_j)))
  table <- list()
  table[[1]] <- index
  
  x_j <- X[j,]
  xj_t_xj <- hyper_par$XtX[j, j] 
  
  log_llh_null <- - n*log(xj_t_xj/2)/2
  lmp_null <- log_m_prior(0, h, max_p)
  
  prob_curr_null <- exp(log_llh_null + lmp_null)
  
  
  C <- prob_curr_null
  # max_C <- log_llh_null + lmp_null
  
  
  p_Pa <- length(Pa_j)
  k_power <- 2^p_Pa
  
  for (i in 2:length(index)) {
    
    k <- index[i]
    
    if (k == 0) {
      Pa_j_new <- Pa_j
      p_Pa_new <- length(Pa_j_new)
      if (p_Pa_new == 0) {
        table[[i]] <- log_llh_null
        next
      }
    }
    else {
      Pa_j_new <- c(Pa_j, k)
      p_Pa_new <- length(Pa_j_new)
    }
    
    
    
    x_pa_j <- X[Pa_j_new,]
    
    V <- hyper_par$XtX[Pa_j_new, Pa_j_new]
    x_paj_xj_t <- hyper_par$XtX[Pa_j_new, j]
    
    if (p_Pa_new == 1) {
      V <- as.matrix(V + inv_g)
    }
    else {
      diag(V) <- diag(V) + inv_g
    }
    
    
    log_llhs <- rep(log_llh_null, k_power)
    
    curr <- rep(0, p_Pa_new)
    
    if (p_Pa_new != p_Pa) {
      curr[p_Pa_new] <- 1
    }
    
    k_matrix_old <- rep(0, p_Pa)
    
    
    
    for (idx in 0:(k_power-1)) {
      
      k_matrix <- Gray_code(idx, k_matrix_old)
      temp_j <- which(k_matrix != k_matrix_old)
      curr[temp_j] <- 1 - curr[temp_j]
      k_matrix_old <- k_matrix
      
      curr_incl <- which(curr == 1)
      Pa_j_incl <- Pa_j_new[curr_incl]
      p_curr <- length(curr_incl)
      
      if (p_curr == 0) {
        next
      }
      
      
      x_paj_xj_t_curr <- x_paj_xj_t[curr_incl]
      L_Vg <- chol(V[curr_incl, curr_incl])
      A <- - n * log((xj_t_xj - sum((forwardsolve(t(L_Vg), diag(p_curr)) %*%
                                       x_paj_xj_t_curr)^2))/2)/2
      
      if (p_curr == 1){
        log_sqrt_det_Sigma <- log(L_Vg)
      }
      else {
        log_sqrt_det_Sigma <- sum(log(diag(L_Vg)))
      }
      
      log_llh <- - p_curr*log(g)/2 - sum(log_sqrt_det_Sigma) + A 
      
      lmp <- log_m_prior(p_curr, h, max_p)
      
      prob_curr <- exp(log_llh + lmp) # - max_C
      
      PIPs[Pa_j_incl] <- PIPs[Pa_j_incl] + prob_curr
      
      if (p_Pa_new != p_Pa) {
        m <- model_encoding(curr[-p_Pa_new])
      }
      else {
        m <- model_encoding(curr)
      }
      
      log_llhs[m] <- log_llh
      
      
      C <- C + prob_curr
      
      # if (prob_curr > 1) {
      #   
      #   PIPs <- PIPs/prob_curr
      #   C <- C/prob_curr
      #   
      #   max_C <- log_llh + lmp # log(prob_curr)
      #   print(max_C)
      #   # print(max_C)
      # }
      
      
    }
    
    
    table[[i]] <- log_llhs
    
    
  }
  
  
  return(list(table = table,
              PIPs = PIPs/C))
  
}






convert_to_binary <- function(n, k = c()) {
  if(n > 1) {
    k <- convert_to_binary(as.integer(n/2), k)
  }
  return(c(k, n %% 2))
  
}

Gray_code <- function(idx, vector) {
  b <- convert_to_binary(idx)
  g <- b
  if (length(g) > 1) {
    for (i in 2:length(b)) {
      g[i] <- (b[i] + b[i-1]) %% 2
    }
  }
  vector[(length(vector)-length(g)+1) : length(vector)] <- g
  
  return(vector)
}






