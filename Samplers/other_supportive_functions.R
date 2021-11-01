# bayes factor

# g-prior
bf_g_L <- function(gamma, hyper_par){
  
  includes <- which(gamma == 1)
  
  n <- hyper_par$n
  p <- hyper_par$p
  yty <- hyper_par$yty
  # XtX <- hyper_par$XtX
  ytX <- hyper_par$ytX
  g <- hyper_par$g
  X <- hyper_par$X
  diag_V <- hyper_par$diag_V
  # y <- hyper_par$y
  # h <- hyper_par$h
  
  g_ratio <- g/(g+1)
  inv_sqrtg1 <- 1/sqrt(g+1)
  
  n_power <- n/2
  
  BF <- rep(NA, p)
  p_gam <- sum(gamma)
  
  # empty model
  if(p_gam == 0){
    A <- yty
    for(j in 1:p){
      tilda_A <- yty - g_ratio * ytX[j]^2/diag_V[j]
      
      A_ratio <- A/tilda_A
      
      BF[j] <- (A_ratio)^n_power * inv_sqrtg1
    }
  }
  else{
    
    zj <- 0
    
    Xg <- X[, includes]
    XgtX <- t(Xg) %*% X
    XgtXg <- XgtX[, includes]  # XtX[includes, includes]
    L_Xg <- chol(XgtXg)
    
    ytXg <- matrix(ytX[includes],nrow = 1,ncol = p_gam)
    
    F <- solve(XgtXg)
    L_F <- chol(F)
    
    ytXgFXgty <- sum((solve(t( L_Xg )) %*% t(ytXg))^2)
    
    A <- yty - ytXgFXgty * g_ratio
    
    # XgtX <- matrix(XtX[includes,], nrow = p_gam, ncol = p)
    # XgtX <- XtX[includes,]
    
    d_vec <- 1 / (diag_V - rowSums((t(XgtX) %*% t(L_F))^2))
    
    ytXFXtxj_vec <- ytXg %*% F %*% XgtX
    tilda_A_vec <- A - d_vec * (ytXFXtxj_vec - ytX)^2 * g_ratio
    
    # print(( d_vec * (ytXFXtxj_vec - ytX)^2 * g_ratio)[1])
    
    BF <- (A / tilda_A_vec)^n_power * inv_sqrtg1
    
    for (j in includes) {
      
      if(p_gam == 1){
        
        tilde_A <- yty
        A_ratio <- tilde_A/A
        
      }
      else{
        
        zj <- zj+1
        
        # dj <- 1/F[zj,zj]
        # tilde_A <- A + (ytXg %*% F[,zj] )^2 * g_ratio * dj
        # A_ratio <- tilde_A/A
        
        A_ratio <- 1 + (ytXg %*% F[,zj] )^2 * g_ratio / (A * F[zj,zj])
        
      }
      
      
      BF[j] <- A_ratio^n_power * inv_sqrtg1
      
    }
  }
  
  return(BF)
}


# independent prior
bf_ind_L <- function(gamma, hyper_par, L_gamma = NULL){
  
  includes <- which(gamma == 1)
  
  n <- hyper_par$n
  p <- hyper_par$p
  yty <- hyper_par$yty
  XtX <- hyper_par$XtX
  ytX <- hyper_par$ytX
  g <- hyper_par$g
  # y <- hyper_par$y
  # h <- hyper_par$h
  X <- hyper_par$X
  # diag_XtX <- hyper_par$diag_XtX
  diag_V <- hyper_par$diag_V # + inv_g
  inv_g <- 1/g
  inv_sqrt_g <- sqrt(inv_g)
  n_power <- n/2
  
  BF <- rep(NA, p)
  p_gam <- sum(gamma)
  
  # empty model
  if(p_gam == 0){
    A <- yty
    for(j in 1:p){
      
      dj <- 1 / diag_V[j]
      
      tilda_A <- yty - ytX[j]^2 * dj
      
      A_ratio <- A/tilda_A
      
      BF[j] <- (A_ratio)^n_power * sqrt(dj) * inv_sqrt_g
    }
  }
  else{
    
    zj <- 0
    
    # Vg <- V[includes, includes]
    
    if (is.null(L_gamma)) {
      Xg <- X[, includes]
      XgtX <- t(Xg) %*% X
      Vg <- XgtX[, includes]
      
      if (p_gam == 1){
        Vg <- Vg + inv_g
      }
      else {
        diag(Vg) <- diag(Vg) + inv_g
      }
      L_Vg <- chol(Vg)
    }
    else {
      L_Vg <- L_gamma
    }
    
    ytXg <- matrix(ytX[includes],nrow = 1,ncol = p_gam)
    
    F <- solve(Vg)
    L_F <- chol(F)
    
    ytXgFXgty <- sum((solve(t( L_Vg )) %*% t(ytXg))^2)
    
    A <- yty - ytXgFXgty
    
    # XgtX <- matrix(XtX[includes,], nrow = p_gam, ncol = p)
    # XgtX <- XtX[includes,]
    
    
    inv_d_vec <- 1 / (diag_V - rowSums((t(XgtX) %*% t(L_F))^2))
    # inv_d_vec <- 1 / (diag(XtX) + inv_g - diag(t(XgtX) %*% F %*%  XgtX))
    inv_d_vec[includes] <- 0
    
    #if (any(inv_d_vec < 0)){
    #  print(inv_d_vec)
    #}
    
    ytXFXtxj_vec <- ytXg %*% F %*% XgtX
    tilda_A_vec <- A - inv_d_vec * (ytXFXtxj_vec - ytX)^2 
    
    
    BF <- sqrt(inv_d_vec) * inv_sqrt_g * (A / tilda_A_vec)^n_power 
    
    
    for (j in includes) {
      
      if(p_gam == 1){
        
        tilde_A <- yty
        A_ratio <- tilde_A/A
        
        BF[j] <-  sqrt(F) * inv_sqrt_g * A_ratio^n_power
        
      }
      else{
        
        zj <- zj+1
        
        inv_dj <- F[zj,zj]
        
        tilde_A <- A + (ytXg %*% F[,zj] )^2 / inv_dj
        A_ratio <- tilde_A/A
        
        # A_ratio <- 1 + (ytXg %*% F[,zj] )^2 / (A * F[zj,zj])
        
        BF[j] <- sqrt(inv_dj) * inv_sqrt_g * A_ratio^n_power
        
      }
      
      
    }
  }
  
  return(BF)
}

compute_Delta <- function(PIPs) {
  2 * (sum(PIPs[PIPs >= 0.5]) + sum(1 - PIPs[PIPs < 0.5]))
}


compute_LA <- function(gamma, hyper_par, t = 1){
  
  p_gam <- sum(gamma)
  
  log_llh <- hyper_par$log_llh(gamma, hyper_par)
  
  log_m_prior <- hyper_par$log_m_prior(p_gam, hyper_par$h, hyper_par$p)
  
  return(list(# curr = gamma,
    llh = log_llh,
    lmp = log_m_prior,
    log_post = t*log_llh + log_m_prior,
    p_gam = p_gam))
}

extract_llh <- function(LAs, n_temp, n_chain) {
  llhs <- matrix(NA, nrow = n_chain, ncol = n_temp)
  
  for (i in 1:n_chain) {
    for (t in 1:n_temp) {
      llhs[i, t] <- LAs[[t]][[i]]$llh
    }
  }
  return(llhs)
}



generate_temperature <- function(taus) {
  cumprod(c(1, exp(-exp(taus))))
}




log_llh_g_L <- function(gamma, hyper_par){
  
  
  p_gam <- sum(gamma)
  
  n <- hyper_par$n
  p <- hyper_par$p
  h <- hyper_par$h
  g <- hyper_par$g
  yty <- hyper_par$yty
  
  includes <- which(gamma == 1)
  
  if(p_gam == 0){
    A <- yty
  }
  else{
    
    ytX <- hyper_par$ytX
    
    Xg <- hyper_par$X[, includes]
    L_Xg <- chol(t(Xg) %*% Xg)
    
    ytXgFXgty <- sum((forwardsolve(t(L_Xg), diag(p_gam)) %*% 
                        as.matrix(ytX[includes]))^2)
    
    A <-  yty - g/(1+g) * ytXgFXgty
    
  }
  
  log_llh <-  - p_gam/2 * log(1+g) - n * log(A)/2
  
  # log_m <- hyper_par$log_m_prior(p_gam, h, p)
  
  # log_post <- log_llh + log_m
  
  return(log_llh)
  
}



log_llh_ind_L <- function(gamma, hyper_par){
  
  p_gam <- sum(gamma)
  
  n <- hyper_par$n
  p <- hyper_par$p
  h <- hyper_par$h
  # h_alpha <- hyper_par$h_alpha
  # h_beta <- hyper_par$h_beta
  g <- hyper_par$g
  yty <- hyper_par$yty
  
  
  includes <- which(gamma == 1)
  
  if(p_gam == 0){
    A <- yty
    L_Vg <- NULL
    sqrt_det_Vg <- 0
  }
  else{
    
    X <- hyper_par$X
    # V <- hyper_par$V
    ytX <- hyper_par$ytX
    
    # Vg <- V[includes, includes]
    Xg <- X[, includes]
    Vg <- t(Xg) %*% Xg
    diag(Vg) <- diag(Vg) + 1/g
    L_Vg <- chol(Vg)
    
    ytXgFXgty <- sum((forwardsolve(t(L_Vg), diag(p_gam)) %*% 
                        as.matrix(ytX[includes]))^2)
    
    A <-  yty - ytXgFXgty
    
    if (p_gam == 1){
      sqrt_det_Vg <- log(L_Vg)
    }
    else {
      sqrt_det_Vg <- sum(log(diag(L_Vg)))
    }
    
  }
  
  log_llh <- - p_gam/2 * log(g) - sqrt_det_Vg - n * log(A)/2
  
  # fixed h
  # log_m <- p_gam * (log(h) - log(1-h))
  
  # binomial-beta prior
  # log_m <- sum(log(1:p_gam)) - sum(log(p - p_gam + h_beta - 0:p_gam))
  # log_m <- hyper_par$log_m_prior(p_gam, h, p) # lbeta(p_gam+h_alpha, p-p_gam+h_beta)
  
  # log_post <- log_llh + log_m
  
  return(log_llh)
}


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






# sampling independently
sample_ind <- function(whe_sam, probs, n_sam = NULL, samples = NULL, log = FALSE){
  
  n <- length(probs)
  
  if (whe_sam) {
    samples <- runif(n) < probs
  }
  
  if (log) {
    prob <- sum(log(probs[samples]))
  }
  else {
    prob <- prod(probs[samples])
  }
  
  return(list(prob = prob,
              sample = samples))
  
}



update_temp_in_LA <- function(LA, new_temp) {
  LA$log_post <- new_temp*LA$llh + LA$lmp
  return(LA)
}




# g/independent prior
update_LA <- function(curr, LA, k, hyper_par, alg_par, PIPs, zeta, t = 1) {
  
  temp <- curr
  
  # p_temp <- LA$p_gam
  log_post_temp <- LA$log_post
  llh_temp <- LA$llh
  lmp_temp <- LA$lmp
  
  # X <- hyper_par$X
  # ytX <- hyper_par$ytX
  # yty <- hyper_par$yty
  # n <- hyper_par$n
  p <- hyper_par$p
  h <- hyper_par$h
  # h_odd <- hyper_par$h_odd
  log_llh <- hyper_par$log_llh
  log_m_prior <- hyper_par$log_m_prior
  bal_fun <- alg_par$bal_fun
  
  k_size <- length(k)
  
  # inv_g <- 1/g
  # n_power <- (n-1)/2
  # inv_sqrtg <- sqrt(inv_g)
  # sqrtg <- sqrt(g)
  
  # initialisation
  prob_prop <- 0
  rev_prob_prop <- 0
  # prob_k_odds <- 1
  JD <- 0
  
  prod_bal_con <- 0
  rev_prod_bal_con <- 0
  
  for (k_j in k) {
    
    temp_2 <- temp
    temp_kj <- temp_2[k_j]
    temp_2[k_j] <- 1 - temp_kj
    
    llh_temp2 <- log_llh(temp_2, hyper_par)
    lmp_temp2 <- log_m_prior(sum(temp_2), h, p)
    # log_post_temp_2 <- log_pi(temp_2, hyper_par)
    log_post_temp_2 <- t*llh_temp2 + lmp_temp2
    
    prob_temp2_temp <- exp(log_post_temp_2-log_post_temp)
    
    mar_eff <- PIPs[k_j]
    odd_k <- (mar_eff/(1-mar_eff))^(2*temp_kj-1)
    
    prob_change <- zeta * bal_fun(prob_temp2_temp * odd_k)
    prob_keep <- (1-zeta) * bal_fun(1)
    
    bal_Const <- prob_change + prob_keep
    prob_change <- prob_change/bal_Const
    prob_keep <- prob_keep/bal_Const
    
    if (runif(1) < prob_change) {
      # chang k_j
      rev_prob_change <- zeta * bal_fun(1 / prob_temp2_temp / odd_k)
      rev_prob_keep <- (1-zeta) * bal_fun(1)
      
      rev_bal_Const <- rev_prob_change + rev_prob_keep
      rev_prob_change <- rev_prob_change/rev_bal_Const
      rev_prob_keep <- rev_prob_keep/rev_bal_Const
      
      prob_prop <- prob_prop + log(prob_change)
      rev_prob_prop <- rev_prob_prop + log(rev_prob_change)
      
      
      # update temp
      temp <- temp_2
      
      # update JD
      JD <- JD + 1
      
      # update contents LA
      # p_temp <- p_temp2
      log_post_temp <- log_post_temp_2
      llh_temp <- llh_temp2
      lmp_temp <- lmp_temp2
      
      # update k odds
      # prob_k_odds <- prob_k_odds * odd_k
      
    }
    else {
      # keep k_j
      rev_prob_change <- zeta * bal_fun(prob_temp2_temp * odd_k)
      rev_prob_keep <- (1-zeta) * bal_fun(1)
      
      rev_bal_Const <- rev_prob_change + rev_prob_keep
      rev_prob_change <- rev_prob_change/rev_bal_Const
      rev_prob_keep <- rev_prob_keep/rev_bal_Const
      
      prob_prop <- prob_prop + log(prob_keep)
      rev_prob_prop <- rev_prob_prop + log(rev_prob_keep)
      
    }
    
    prod_bal_con <- prod_bal_con + log(bal_Const)
    rev_prod_bal_con <- rev_prod_bal_con + log(rev_bal_Const)
    
  }
  
  
  # update LA
  LA_prop <- list(# curr = temp,
    llh = llh_temp,
    lmp = lmp_temp,
    log_post = log_post_temp,
    p_gam = sum(temp))
  
  prop <- temp
  
  # acc_prob <- log_post + rev_prob_prop - LA$llh - prob_prop
  # acc_rate <- min(1, exp(acc_prob) * prob_k_odds)
  # cat(acc_rate -  min(1, exp(prod_bal_con-rev_prod_bal_con)), "\n")
  
  acc_rate <- min(1, exp(prod_bal_con-rev_prod_bal_con))
  # cat(exp(LA_prop$llh-LA$llh), exp(prob_prop), exp(rev_prob_prop), exp(prod_bal_con-rev_prod_bal_con), "\n")
  
  return(list(prop = temp, 
              LA_prop = LA_prop,
              JD = JD,
              acc_rate = acc_rate))
}


