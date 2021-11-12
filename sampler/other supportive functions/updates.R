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


