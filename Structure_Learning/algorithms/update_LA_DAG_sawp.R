# add reversing move
update_LA_DAG <- function(LA, k, hyper_par, bal_fun, PIPs, thinning_rate, omega) {
  
  temp <- LA$curr
  LA_temp <- LA
  
  # p_temp <- LA$p_gam
  log_post_temp <- LA$log_post
  llh_temp <- LA$llh
  lmp_temp <- LA$lmp
  
  max_p <- hyper_par$max_p
  p <- hyper_par$p
  h <- hyper_par$h
  
  log_llh_update <- hyper_par$log_llh_update
  log_m_prior <- hyper_par$log_m_prior
  swap_idx <- hyper_par$swap_idx
  
  
  k_swap <- swap_idx[k]
  swaps <- k_swap %in% k
  num_swaps <- as.integer(sum(swaps)/2)
  
  
  k_size <- length(k)
  M <- get_moves()
  
  # thinned_k <- list()
  # thinned_k_size <- 0
  # 
  # used_idx <- c()
  # 
  # for (j in 1:k_size) {
  # 
  #   if (k[j] %in% used_idx) { next }
  # 
  #   if (swaps[j]) {
  #     k_temp <- c(k[j], k_swap[j])
  #     used_idx <- c(used_idx, k_swap[j])
  #   }
  #   else {
  #     k_temp <- k[j]
  #   }
  # 
  #   if (runif(1) < thinning_rate) {
  # 
  #     thinned_k_size <- thinned_k_size+1
  #     thinned_k[[thinned_k_size]] <- k_temp
  # 
  #   }
  # 
  # }
  
  
  # thinned_k <- k[which(runif(k_size) < thinning_rate)]
  # thinned_k_size <- length(thinned_k)
  
  
  grouped_k <- matrix(k, nrow = k_size, ncol = 2)
  
  grouped_k[,2] <- ifelse(swaps, k_swap, Inf)
  
  grouped_k <- grouped_k[grouped_k[,2] > grouped_k[,1],]
  
  thinned_k <- grouped_k[runif(nrow(grouped_k)) < thinning_rate,]

  thinned_k_size <- nrow(thinned_k)
  
  if (is.null(thinned_k_size)) {
    thinned_k <- matrix(thinned_k, nrow = 1, ncol = 2)
    thinned_k_size <- 1
  }
  
  
  # initialisation
  prob_prop <- 0
  rev_prob_prop <- 0
  prob_k_odds <- 1
  JD <- 0
  
  prod_bal_con <- 0
  rev_prod_bal_con <- 0
  
  chage_times <- 0
  swap_times <- 0
  
  if (thinned_k_size >= 1) {
    new_order <- sample(thinned_k_size)
    
    for (j in 1:thinned_k_size) {
      
      # (0,0) -> change
      # (1,0) -> swap/change
      # (0,1) -> swap
      
      new_j <- new_order[j]
      k_j <- thinned_k[new_j,]
      
      log_post_temp_old <- log_post_temp
      
      
      if (!is.infinite(k_j[2])) {
        
        temp_kj <- temp[k_j[1]]
        temp_kj_swap <- temp[k_j[2]]
        
        # link k_j and kj_swap
        
        LA_temps <- list()
        prob_change_k_ratio <- rep(1,4)
        odd_k_change <- rep(1, 4)
        
        mar_eff_kj <- PIPs[k_j[1]]
        mar_eff_kj_swap <- PIPs[k_j[2]]
        
        prob_change <- rep((1-omega)^2 * bal_fun(1), 4)
        
        # 4 moves
        for (i in 2:4) {
          
          temp_change <- temp
          temp_change[k_j[1]] <- abs(M[i,1] - temp_kj)
          temp_change[k_j[2]] <- abs(M[i,2] - temp_kj_swap)
          
          # check fair graph
          if ((temp_change[k_j[1]] + temp_change[k_j[2]]) < 2) {
            if (is.DAG_adjmat(temp_change)) {
              
              # LA_temp_change <- hyper_par$log_llh(list(curr = temp_change), hyper_par)
              LA_temp_change <- log_llh_update(changes = ceiling(k_j/p), LA_temp, list(curr = temp_change), hyper_par)
              
              
              LA_temps[[i]] <- LA_temp_change
              llh_temp_change <- LA_temp_change$llh
              lmp_temp_change <- log_m_prior(LA_temp_change$p_gam, h, max_p)
              
              
              log_post_temp_change <- llh_temp_change + lmp_temp_change
              
              prob_change_ratio <- exp(log_post_temp_change-log_post_temp)
              
              
              odd_k_change[i] <- (mar_eff_kj/(1-mar_eff_kj))^((2*temp_kj-1)*M[i,1]) * (mar_eff_kj_swap/(1-mar_eff_kj_swap))^((2*temp_kj_swap-1)*M[i,2])
              
              
              
              # eta <- (1-temp) * hyper_par$A + temp * hyper_par$D
              # p1 <- sample_ind_DAG(FALSE, probs = eta, samples = k, log = TRUE)$prob
              # eta_prop <- (1-temp_change) * hyper_par$A + temp_change * hyper_par$D
              # p2 <- sample_ind_DAG(FALSE, probs = eta_prop, samples = k, log = TRUE)$prob
              # cat("k_prob", odd_k_change[i], exp(p2-p1), "\n")
              
              
              
              prob_change_k_ratio[i] <- prob_change_ratio*odd_k_change[i]
              
              prob_change[i] <- omega^sum(M[i,]) * (1-omega)^(2-sum(M[i,])) * bal_fun(prob_change_k_ratio[i])
              
              
              
            }
            else {
              llh_temp_change <- 0
              lmp_temp_change <- -Inf
              
              prob_change[i] <- 0
            }
          }
          else {
            llh_temp_change <- 0
            lmp_temp_change <- -Inf
            
            prob_change[i] <- 0
          }
          
        }
        
        
        bal_Const <- sum(prob_change)
        prob_change <- prob_change/bal_Const
        
        
        change_idx <- which(runif(1) < cumsum(prob_change))[1]   
        
        
        if (change_idx == 1) {
          
          # rev_prob_change <- omega * bal_fun(prob_temp2_temp * odd_k)
          # rev_prob_keep <- (1-omega) * bal_fun(1)
          # 
          # rev_bal_Const <- rev_prob_change + rev_prob_keep
          # rev_prob_change <- rev_prob_change/rev_bal_Const
          # rev_prob_keep <- rev_prob_keep/rev_bal_Const
          
          rev_bal_Const <- bal_Const
          
          prob_prop <- prob_prop + log(prob_change[1])
          rev_prob_prop <- rev_prob_prop + log(prob_change[1])
          
          
        }
        else {
          
          # change gamma
          # update contents LA
          
          LA_temp <- LA_temps[[change_idx]]
          temp <- LA_temp$curr
          
          llh_temp <- LA_temp$llh
          lmp_temp <- log_m_prior(LA_temp$p_gam, h, max_p)
          log_post_temp <- llh_temp + lmp_temp
          
          # update JD
          JD <- JD + sum(M[change_idx,])
          
          # update change count
          if (sum(M[change_idx,]) == 1) {
            chage_times <- chage_times + 1
          }
          else {
            swap_times <- swap_times + 1
          }
          
          
          # update k odds
          prob_k_odds <- prob_k_odds * odd_k_change[change_idx]
          
          rev_prob_change_k_ratio <- prob_change_k_ratio/prob_change_k_ratio[change_idx]
          
          if (sum(is.nan(rev_prob_change_k_ratio)) > 0) {
            rev_prob_change_k_ratio[which(is.nan(rev_prob_change_k_ratio))] <- 1
          }
                    
          omega_vec <- get_omega_vec(omega, change_idx, M)
          rev_prob_change <- bal_fun(rev_prob_change_k_ratio) * omega_vec
          rev_bal_Const <- sum(rev_prob_change)
          rev_prob_change <- rev_prob_change/rev_bal_Const
          
          prob_prop <- prob_prop + log(prob_change[change_idx])
          rev_prob_prop <- rev_prob_prop + log(rev_prob_change[1])
          
          
          odd_k_now <- odd_k_change[change_idx]
          prop_prob_temp <- log(prob_change[change_idx])
          rev_prop_prob_temp <- log(rev_prob_change[1])
          
        }
        
        prod_bal_con <- prod_bal_con + log(bal_Const)
        rev_prod_bal_con <- rev_prod_bal_con + log(rev_bal_Const)
        
        
        
        
      }
      else {
        
        k_j <- k_j[1]
        
        temp_change <- temp
        temp_kj <- temp_change[k_j]
        temp_change[k_j] <- 1 - temp_kj
        
        if (is.DAG_adjmat(temp_change)) {
          
          # LA_temp_change <- list(curr = temp_change)
          LA_temp_change <- log_llh_update(changes = ceiling(k_j/p), LA_temp, list(curr = temp_change), hyper_par)
          
          # LA_temp_change <- hyper_par$log_llh(list(curr = temp_change), hyper_par)
          
          llh_temp_change <- LA_temp_change$llh
          lmp_temp_change <- log_m_prior(sum(temp_change), h, max_p)
          # print(lmp_temp_change)
          
        }
        else {
          llh_temp_change <- 0
          lmp_temp_change <- -Inf
        }
        
        
        log_post_temp_change <- llh_temp_change + lmp_temp_change
        
        
        prob_change_ratio <- exp(log_post_temp_change-log_post_temp)
        
        
        mar_eff <- PIPs[k_j]
        odd_k_change <- (mar_eff/(1-mar_eff))^(2*temp_kj-1)
        
        prob_change <- omega * bal_fun(prob_change_ratio * odd_k_change)
        prob_keep <- (1-omega) * bal_fun(1)
        
        bal_Const <- prob_change + prob_keep
        prob_change <- prob_change/bal_Const
        prob_keep <- prob_keep/bal_Const
        
        # cat(k_j, omega, odd_k_change, prob_change, prob_keep, "\n")
        
        if (runif(1) < prob_change) {
          
          # chang k_j
          rev_prob_change <- omega * bal_fun(1 / prob_change_ratio / odd_k_change)
          rev_prob_keep <- (1-omega) * bal_fun(1)
          
          rev_bal_Const <- rev_prob_change + rev_prob_keep
          rev_prob_change <- rev_prob_change/rev_bal_Const
          rev_prob_keep <- rev_prob_keep/rev_bal_Const
          
          prob_prop <- prob_prop + log(prob_change)
          rev_prob_prop <- rev_prob_prop + log(rev_prob_change)
          
          
          # update temp
          temp <- temp_change
          
          # update JD
          JD <- JD + 1
          
          # update contents LA
          # p_temp <- p_temp2
          log_post_temp <- log_post_temp_change
          llh_temp <- llh_temp_change
          lmp_temp <- lmp_temp_change
          
          LA_temp <- LA_temp_change
          
          # update k odds
          prob_k_odds <- prob_k_odds * odd_k_change
          
          
          # odd_k_now <- odd_k_change
          # prop_prob_temp <- log(prob_change)
          # rev_prop_prob_temp <- log(rev_prob_change)
          
        }
        else {
          
          # keep k_j
          # rev_prob_change <- omega * bal_fun(prob_change_ratio * odd_k_change)
          # rev_prob_keep <- (1-omega) * bal_fun(1)
          # 
          # rev_bal_Const <- rev_prob_change + rev_prob_keep
          # rev_prob_change <- rev_prob_change/rev_bal_Const
          # rev_prob_keep <- rev_prob_keep/rev_bal_Const
          
          rev_bal_Const <- bal_Const
          
          prob_prop <- prob_prop + log(prob_keep)
          rev_prob_prop <- rev_prob_prop + log(prob_keep)
          
          # odd_k_now <- 1
          # prop_prob_temp <- log(prob_keep)
          # rev_prop_prob_temp <- log(prob_keep)
          
        }
        
        prod_bal_con <- prod_bal_con + log(bal_Const)
        rev_prod_bal_con <- rev_prod_bal_con + log(rev_bal_Const)
        
        
        # cat(bal_Const/rev_bal_Const, exp(log_post_temp - log_post_temp_old) * odd_k_now * exp(rev_prop_prob_temp-prop_prob_temp),"\n")
        
        
      }
      
      
      
      
    }
  }
  
  
  
  
  # update LA
  # LA_prop <- list(
  #   curr = temp,
  #   llh = llh_temp,
  #   lmp = lmp_temp,
  #   log_post = log_post_temp,
  #   p_gam = sum(temp))
  LA_prop <- LA_temp
  LA_prop$lmp <- lmp_temp
  LA_prop$log_post <- log_post_temp
  
  
  # prop <- temp
  
  # acc_prob <- log_post_temp + rev_prob_prop - LA$log_post - prob_prop
  # acc_rate <- min(1, exp(acc_prob) * prob_k_odds)
  # cat(acc_rate -  min(1, exp(prod_bal_con-rev_prod_bal_con)), "\n")
  
  
  
  # eta <- (1-LA$curr) * hyper_par$A + LA$curr * hyper_par$D
  # p1 <- sample_ind_DAG(FALSE, probs = eta, samples = k, log = TRUE)$prob
  # eta_prop <- (1-LA_prop$curr) * hyper_par$A + LA_prop$curr * hyper_par$D
  # p2 <- sample_ind_DAG(FALSE, probs = eta_prop, samples = k, log = TRUE)$prob
  # cat(prob_k_odds/ exp(p2-p1), "\n")
  # 
  # print(swap_times)
  # if (is.nan(rev_prod_bal_con)) {
  #   return(LA_prop)
  # }
  acc_rate <- min(1, exp(prod_bal_con-rev_prod_bal_con))
  # cat("accrate", exp(prod_bal_con-rev_prod_bal_con), "\n")
  # cat(exp(LA_prop$llh-LA$llh), exp(prob_prop), exp(rev_prob_prop), exp(prod_bal_con-rev_prod_bal_con), "\n")
  
  # cat(acc_rate, prod_bal_con, rev_prod_bal_con, "\n")
  
  # cat(exp(prod_bal_con-rev_prod_bal_con), exp(LA_prop$log_post - LA$log_post)*exp(rev_prob_prop-prob_prop)*prob_k_odds, "\n")
  
  # print(omega)
  # print(JD)
  # print(acc_rate)
  
  # stop()
  
  
  return(list(# prop = temp, 
    LA_prop = LA_prop,
    JD = JD,
    acc_rate = acc_rate, 
    thinned_k_size = thinned_k_size))
  
}



get_moves <- function(){
  return(matrix(c(rep(c(0,1), each = 2), rep(c(0,1), times = 2)), 4, 2))
  # return(matrix(c(0,1,1,0,1,1), byrow = TRUE, 3, 2))
}


get_omega_vec <- function(omega, change_idx, M) {
  
  M <- colSums( abs(t(M) - M[change_idx,]) )
  return(omega^M *(1-omega)^(2-M))
  
}


