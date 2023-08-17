PARNI_proposal <- function(ALL, k, hyper_par, alg_par, PIPs, omega, zeta, eta, t = 1) {
  
  # barplot(PIPs, ylim = c(0,1))
  temp <- ALL$curr
  ALL_temp <- ALL
  
  # p_temp <- ALL$p_gam
  # llh_temp <- ALL$llh
  # lmp_temp <- ALL$lmp
  # log_post_temp <- t * llh_temp + lmp_temp
  
  p_z <- hyper_par$p_z
  p <- hyper_par$p
  
  k_size <- length(k)
  
  thinned_k <- k[runif(k_size) < zeta]
  thinned_k_size <- length(thinned_k)
  
  approx_log_post <- hyper_par$approx_log_post
  
  
  
  ALL_temp <- try(approx_log_post(ALL_temp, hyper_par, eta), silent = TRUE)
  
  if (is.list(ALL_temp)) {
    
    llh_temp <- ALL_temp$approx_llh
    lmp_temp <- ALL_temp$lmp
    log_post_temp <- t * llh_temp + lmp_temp
    
    
    
    bal_fun <- alg_par$bal_fun
    
    
    
    # initialisation
    prob_prop <- 0
    rev_prob_prop <- 0
    prob_k_odds <- 1
    JD <- 0
    
    prod_bal_con <- 0
    rev_prod_bal_con <- 0
    
    for (k_j in thinned_k) {
      
      temp_2 <- temp
      temp_kj <- temp_2[k_j]
      temp_2[k_j] <- 1 - temp_kj
      
      ALL_temp2 <- ALL_temp
      ALL_temp2$curr <- temp_2
      ALL_temp2$p_gam <- sum(temp_2) # can be facilitate
      
      # cat(k_j, ALL_temp2$p_gam, "\n")
      
      # update log posterior approximated by ALA
      ALL_temp2 <- try(approx_log_post(ALL_temp2, hyper_par, eta), silent = TRUE)
      if (is.list(ALL_temp2)) {
        llh_temp2 <- ALL_temp2$approx_llh
        lmp_temp2 <- ALL_temp2$lmp
      }
      else {
        llh_temp2 <- -Inf
        lmp_temp2 <- 0
      }
      
      
      
      
      log_post_temp_2 <- t*llh_temp2 + lmp_temp2
      
      prob_temp2_temp <- exp(log_post_temp_2-log_post_temp)
      
      # print(prob_temp2_temp)
      
      mar_eff <- PIPs[k_j]
      odd_k <- (mar_eff/(1-mar_eff))^(2*temp_kj-1)
      
      prob_change <- omega * bal_fun(prob_temp2_temp * odd_k)
      prob_keep <- (1-omega) * bal_fun(1)
      
      
      # if (is.infinite(prob_change)) {
      #   print(prob_change)
      #   bal_Const <- 1
      #   prob_change <- 1
      #   prob_keep <- 0
      # }
      # else {
      #   bal_Const <- prob_change + prob_keep
      #   prob_change <- prob_change/bal_Const
      #   prob_keep <- prob_keep/bal_Const
      # }
      
      bal_Const <- prob_change + prob_keep
      prob_change <- prob_change/bal_Const
      prob_keep <- prob_keep/bal_Const
      
      
      # print(prob_change)
      
      
      if (runif(1) < prob_change) {
        # chang k_j
        rev_prob_change <- omega * bal_fun(1 / prob_temp2_temp / odd_k)
        rev_prob_keep <- (1-omega) * bal_fun(1)
        
        rev_bal_Const <- rev_prob_change + rev_prob_keep
        rev_prob_change <- rev_prob_change/rev_bal_Const
        rev_prob_keep <- rev_prob_keep/rev_bal_Const
        
        prob_prop <- prob_prop + log(prob_change)
        rev_prob_prop <- rev_prob_prop + log(rev_prob_change)
        
        
        # update temp
        temp <- temp_2
        
        # update JD
        JD <- JD + 1
        
        # update contents ALL
        # p_temp <- p_temp2
        log_post_temp <- log_post_temp_2
        llh_temp <- llh_temp2
        lmp_temp <- lmp_temp2
        ALL_temp <- ALL_temp2
        # update k odds
        prob_k_odds <- prob_k_odds * odd_k
        
      }
      else {
        # keep k_j
        rev_prob_change <- omega * bal_fun(prob_temp2_temp * odd_k)
        rev_prob_keep <- (1-omega) * bal_fun(1)
        
        rev_bal_Const <- rev_prob_change + rev_prob_keep
        rev_prob_change <- rev_prob_change/rev_bal_Const
        rev_prob_keep <- rev_prob_keep/rev_bal_Const
        
        prob_prop <- prob_prop + log(prob_keep)
        rev_prob_prop <- rev_prob_prop + log(rev_prob_keep)
        
      }
      
      prod_bal_con <- prod_bal_con + log(bal_Const)
      rev_prod_bal_con <- rev_prod_bal_con + log(rev_bal_Const)
      
    }
    
    
    # ALL_prop
    ALL_prop <- ALL_temp
    
    ALL_prop <- hyper_par$log_post(ALL_prop, hyper_par)
    # log_post_prop <- t*ALL_prop$llh + ALL_prop$lmp
    
    # prop <- temp
    # acc_prob <- log_post + rev_prob_prop - ALL$llh - prob_prop
    # acc_rate <- min(1, exp(acc_prob) * prob_k_odds)
    # cat(acc_rate -  min(1, exp(prod_bal_con-rev_prod_bal_con)), "\n")
    
    # print(ALL_prop$ala2_llh)
    
    # cat(ALL_prop$llh, ALL_temp$approx_llh, "\n")
    
    # acc_rate <- min(1, exp(prod_bal_con-rev_prod_bal_con))
    acc_rate <- min(1, exp(t*ALL_prop$llh + ALL_prop$lmp - t*ALL$llh - ALL$lmp + rev_prob_prop - prob_prop) * prob_k_odds)
    # cat(exp(prod_bal_con-rev_prod_bal_con), exp(ALL_prop$llh + ALL_prop$lmp - ALL$llh - ALL$lmp + rev_prob_prop - prob_prop) * prob_k_odds, "\n")
    
    # print("123")
    # print(which(ALL$curr == 1))
    # print(which(ALL_temp$curr == 1))
    
    
    # print(JD)
    # print(acc_rate)
    
  }
  else {
    
    ALL_prop <- ALL
    JD <- 0
    acc_rate <- 1
    rev_prob_prop <- 0
    prob_prop <- 0
    
  }
  
  
  
  return(list(ALL_prop = ALL_prop,
              JD = JD,
              acc_rate = acc_rate,
              prop_prob_ratio = rev_prob_prop - prob_prop,
              thinned_k_size = thinned_k_size))
}
