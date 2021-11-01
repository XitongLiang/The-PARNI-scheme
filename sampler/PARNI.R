# Pointwise implementation of Adaptive Random Neighbourhood Informed Sampler (Parallel Tempering)
PARNI <- function(alg_par, hyper_par){
  
  # initialisation of alg_par
  N <- alg_par$N   # number of iterations
  Nb <- alg_par$Nb     # number of burn-in interations
  Nl <- alg_par$Nl
  n_chain <- alg_par$n_chain    # number of multiple chain
  n_temp <- alg_par$n_temp
  store_chains <- alg_par$store_chains
  verbose <- alg_par$verbose
  f <- alg_par$eval_f
  taus <- alg_par$init_taus
  target_swap <- alg_par$target_swap
  phi_swap <- alg_par$phi_swap
  omega_adap <- alg_par$omega_adap
  omega <- alg_par$omega_init
  omega_par <- alg_par$omega_par
  eps <- alg_par$eps
  full_adap <- alg_par$full_adap
  kappa <- alg_par$kappa
  bal_fun <- alg_par$bal_fun
  use_rb <- alg_par$use_rb
  
  
  # full adaptation
  if (full_adap) {
    Nl <- N
  }
  else {
    if (is.null(Nl)) {
      Nl <- Nb
    }
  }
  
  if (is.null(use_rb)){
    use_rb <- TRUE
  }
  
  
  
  # initialisation of hyper_para
  # n <- hyper_par$n       # number of datapoint
  p <- hyper_par$p       # number of regressors
  g <- hyper_par$g       # g
  h <- hyper_par$h       # model prior parameter
  
  
  # h prior type
  if (length(h) == 1) {
    # fixed h
    h_exp <- h
    hyper_par$h_odd <- function(I, h, ...) {ifelse(I, h/(1-h), (1-h)/h)}
    h_til <- function(h, ...) {h}
    log_m_prior <- function(p_gam, h, p) {p_gam * (log(h) - log(1-h))}
  }
  else {
    # binomial-beta
    h_alpha <- h[1]
    h_beta <- h[2]
    h_exp <- h_alpha/(h_alpha+h_beta)
    hyper_par$h_odd <- function(I, h, p_gam, p){
      ifelse(I, (p_gam+h[1])/(p-p_gam-1+h[2]), (p-p_gam+h[2])/(p_gam-1+h[1]))}
    h_til <- function(h, curr, p_gam) {(p_gam - curr + h[1])/(p + h[1] + h[2] - 1)}
    log_m_prior <- function(p_gam, h, p) {lbeta(p_gam+h[1],p-p_gam+h[2])}
  }
  
  hyper_par$log_m_prior <- log_m_prior
  
  # initialisation
  PIPs <- rep(h_exp, p)
  A <- sapply(PIPs/(1-PIPs), min, 1)
  A <- matrix(A, byrow = TRUE, nrow = n_temp, ncol = p)
  D <- sapply((1-PIPs)/PIPs, min, 1)
  D <- matrix(D, byrow = TRUE, nrow = n_temp, ncol = p)
  
  PIPs <- matrix(PIPs, byrow = TRUE, nrow = n_temp, ncol = p)
  
  # g prior type
  if (hyper_par$g_prior_type == "g") {
    
    bf <- bf_g_L
    log_llh <- log_llh_g_L
    
  }
  else if (hyper_par$g_prior_type == "ind") {
    inv_g <- 1/g
    
    bf <- bf_ind_L
    log_llh <- log_llh_ind_L
    
  }
  
  hyper_par$log_llh <- log_llh
  
  # temps   vector of temperatures
  # taus    vector of controlling the temperatures
  
  if (is.null(taus)) {
    taus <- rep(-4, n_temp-1)
  }
  
  if (n_temp > 1) {
    temps <- generate_temperature(taus)
    temps_trace <- matrix(temps, byrow = TRUE, nrow = N+1, ncol = n_temp)
    taus_trace <- matrix(taus, byrow = TRUE, nrow = N+1, ncol = n_temp-1)
    temp_infs <- matrix(0, nrow = 2, ncol = n_temp-1)
  }
  else {
    temps <- 1
    temps_trace <- NULL
    taus_trace <- NULL
    temp_infs <- NULL
  }
  
  
  # omegas
  omegas <- list()
  
  if (omega_adap == "f") {
    
    omega_curr <- omega
    
  }
  else if (omega_adap == "kw") {
    
    phi_a <- omega_par[1]
    phi_c <- omega_par[2]
    
    n_pos <- floor(n_chain/2)
    n_neg <- n_chain - n_pos
    citer <- 1
    
    omega <- matrix(inv_logit_e(logit_e(omega, eps) + c(0, -citer, citer), eps), 
                    nrow = n_temp, ncol = 3, byrow = TRUE)
    omega_temp <- matrix(omega[1,], byrow = TRUE, nrow = N+1, ncol = 3)
  
  }
  else if (omega_adap == "rm") {
    
    phi_a <- omega_par[1]
    target_omega <- omega_par[2]
    
    omega <- rep(omega, n_temp)
    omega_temp <- rep(omega[1], N+1)
    
  }
  
  
  
  
  
  # swap_acc_times
  swap_acc_times <- 0
  
  # initial model
  currs <- list() # matrix(runif(n = p*n_chain) < h_exp, nrow = n_chain, ncol = p)# rep(0, p)
  chains <- list()
  LAs <- list()
  # llhs_curr <- matrix(NA, nrow = n_chain, ncol = n_temp)
  # prior_curr <- matrix(NA, nrow = n_chain, ncol = n_temp)
  log_posts <- list() # matrix(NA, nrow = N+1, ncol = n_chain)
  model_sizes <- list() # matrix(NA, nrow = N+1, ncol = n_chain)
  infs <- matrix(0, byrow = TRUE, nrow = 2, ncol = p+1)
  
  if (use_rb) {
    BF_temp <- list()
    h_tils <- list()
  }
  else {
    sum_PIPs_temp <- list()
  }
  
  
  for (t in 1:n_temp) {
    if (is.null(alg_par$gamma_init)) {
      currs_temp <- matrix(runif(n = p*n_chain) < h_exp, nrow = n_chain, ncol = p)
    }
    else {
      currs_temp <- matrix(alg_par$gamma_init, nrow = n_chain, ncol = p)
    }
    log_posts_temp <- matrix(NA, nrow = N+1, ncol = n_chain)
    model_sizes_temp <- matrix(NA, nrow = N+1, ncol = n_chain)
    
    if (use_rb){
      BF_temp_temp <- matrix(NA, nrow = n_chain, ncol = p)
      h_tils_temp <- matrix(NA, nrow = n_chain, ncol = p)
    }
    else {
      sum_PIPs_temp_temp <- matrix(NA, nrow = n_chain, ncol = p)
    }
    
    LAs_temp <- list()
    
    for (i in 1:n_chain) {
      curr <- currs_temp[i,]
      LAs_temp_temp <- compute_LA(curr, hyper_par, t = temps[t])
      LAs_temp[[i]] <- LAs_temp_temp
      p_gam <- LAs_temp_temp$p_gam # sum(curr)
      model_sizes_temp[1, i] <- p_gam
      # llhs_curr[i, t] <- LA_temps_temp$llh # log_llh(curr, hyper_par) 
      # prior_curr[i, t] <- LA_temps_temp$lmp # log_m_prior(p_gam, h, p)
      log_posts_temp[1, i] <- LAs_temp_temp$log_post # temps[t] * llhs_curr[i, t] + prior_curr[i, t]
      
      if (use_rb) {
        BF_temp_temp[i, ] <- bf(curr, hyper_par)
        h_tils_temp[i, ] <- h_til(h, curr, p_gam)
      }
      else {
        sum_PIPs_temp_temp[i, ] <- curr
      }
      
      
      if ((t == 1) && store_chains) {
        chains[[i]] <- matrix(curr, nrow = N+1, ncol = p, byrow = TRUE)
      }
    }
    
    currs[[t]] <- currs_temp
    log_posts[[t]] <- log_posts_temp
    model_sizes[[t]] <- model_sizes_temp
    
    if (use_rb) {
      BF_temp[[t]] <- BF_temp_temp
      h_tils[[t]] <- h_tils_temp
    }
    else {
      sum_PIPs_temp[[t]] <- sum_PIPs_temp_temp
    }
    
    LAs[[t]] <- LAs_temp
    
    if (omega_adap != "f") {
      omegas[[t]] <- omega_temp
    }
    
  }
  
  
  # initialisation
  acc_times <- rep(0, n_temp)
  mut <- rep(0, n_temp)
  ESJD <- rep(0, n_temp)
  ESJD_swap <- 0
  ESJD_temp <- rep(0, n_temp-1)
  acc_rate <- 0
  estm_PIPs <- matrix(0, nrow = n_temp, ncol = p)
  JD_swap <- matrix(0, nrow = n_chain, ncol = n_temp-1)
  
  # eval_f
  eval_f <- !is.null(f)
  sum_f <- 0
  
  if (use_rb) {
    # store sum of BF
    Bayes_fac <- matrix(0, nrow = n_temp, ncol = p)
  }
  else {
    sum_PIPs <- matrix(0, nrow = n_temp, ncol = p)
  }
  
  
  start.time1 <- Sys.time()
  
  if ((verbose)) {
    pb <- txtProgressBar(min = 0, max = N, style = 3)
  }
  
  for (iter in 1:N) {
    
    if (verbose) {
      setTxtProgressBar(pb, iter)
    }
    
    if (iter == (Nb + 1)) {
      start.time2 <- Sys.time()
    }
    
    
    if (omega_adap == "kw") {
      # group separation
      omega_idx <- 1:n_chain %in% sample(n_chain, n_pos)
    
      ESJD_pos <- rep(0, n_temp)
      ESJD_neg <- rep(0, n_temp)
      
    }
    else if (omega_adap == "rm") {
      
      omega_acc_rate <- rep(0, n_temp)
      
    }
    
    
    
    for (i in 1:n_chain) {
      
      if (n_temp > 1) {
        # swap move
        swap_idx <- sample(1:(n_temp-1), size = 1)
        
        # acc_rate_swap <- min(1, exp((temps[swap_idx] - temps[swap_idx+1]) * (llhs_curr[i,swap_idx+1] - llhs_curr[i,swap_idx])))
        acc_rate_swap <- min(1, exp((temps[swap_idx] - temps[swap_idx+1]) * (LAs[[swap_idx+1]][[i]]$llh - LAs[[swap_idx]][[i]]$llh)))
        
        if (runif(1) < acc_rate_swap) {
          
          curr_temp <- currs[[swap_idx]][i,]
          currs[[swap_idx]][i,] <- currs[[swap_idx+1]][i,]
          currs[[swap_idx+1]][i,] <- curr_temp
          
          # swap LAs
          LAs_temp <- LAs[[swap_idx]][[i]]
          LAs[[swap_idx]][[i]] <- LAs[[swap_idx+1]][[i]]
          LAs[[swap_idx+1]][[i]] <- LAs_temp
          
          # update log post in LAs
          LAs[[swap_idx]][[i]] <- update_temp_in_LA(LAs[[swap_idx]][[i]], temps[swap_idx])
          LAs[[swap_idx+1]][[i]] <- update_temp_in_LA(LAs[[swap_idx+1]][[i]], temps[swap_idx+1])
          
          if (use_rb) {
            # swap bayes factor and h_til
            BF_temp_temp <- BF_temp[[swap_idx]][i,]
            BF_temp[[swap_idx]][i,] <- BF_temp[[swap_idx+1]][i,]
            BF_temp[[swap_idx+1]][i,] <- BF_temp_temp
            
            h_tils_temp <- h_tils[[swap_idx]][i,]
            h_tils[[swap_idx]][i,] <- h_tils[[swap_idx+1]][i,]
            h_tils[[swap_idx+1]][i,] <- h_tils_temp
          }
          else {
            sum_PIPs_temp_temp <- sum_PIPs_temp[[swap_idx]][i,]
            sum_PIPs_temp[[swap_idx]][i,] <- sum_PIPs_temp[[swap_idx+1]][i,]
            sum_PIPs_temp[[swap_idx+1]][i,] <- sum_PIPs_temp_temp
          }
          
          
          if (iter > Nb) {
            swap_acc_times <- swap_acc_times + 1
          }
        }
        
        # update ESJD
        if (iter > Nb) {
          ESJD_swap <- ESJD_swap + acc_rate_swap * sum(abs(currs[[swap_idx]][i,] - currs[[swap_idx+1]][i,]))
          temp_infs[1, swap_idx] <- temp_infs[1, swap_idx] + 1
          temp_infs[2, swap_idx] <- temp_infs[2, swap_idx] + acc_rate_swap
        }
      }
      
      for (t in 1:n_temp) {
        curr <- currs[[t]][i,]
        LA <- LAs[[t]][[i]]
        
        
        if (omega_adap == "kw") {
          omega_idx_curr <- omega_idx[i]
          omega_curr <- ifelse(omega_idx_curr, omega[t, 3], omega[t, 2])
        }
        else if (omega_adap == "rm") {
          omega_curr <- omega[t]
        }
        
        log_post_curr <- LA$log_post 
        p_gam <- LA$p_gam
        
        eta <- abs(curr-1) * A[t,] + curr * D[t,]
        neighs <- sample_ind(TRUE, n_sam = NULL, probs = eta, samples = NULL, log = TRUE)
        k <- which(neighs$sample)
        k_size <- length(k)
        
        if (k_size > 0) {
          
          if (k_size == 1) {
            k_ordered <- k
          }
          else {
            k_ordered <- sample(k)
          }
          
          updates <- update_LA(curr, LA, k_ordered, hyper_par, alg_par, PIPs[t,], omega_curr, temps[t])
          prop <- updates$prop
          JD <- updates$JD
          acc_rate <- updates$acc_rate
          
          if (JD > 0) {
            if (runif(1) < acc_rate){
              
              curr <- prop
              currs[[t]][i,] <- curr
              LA_prop <- updates$LA_prop
              log_post_curr <- LA_prop$log_post
              p_gam <- LA_prop$p_gam
              LAs[[t]][[i]] <- LA_prop
              
              # update bayes factors
              if (iter <= Nl){
                
                if (use_rb) {
                  BF_temp[[t]][i, ] <- bf(curr, hyper_par)
                  h_tils[[t]][i, ] <- h_til(h, curr, p_gam)
                }
                else {
                  sum_PIPs_temp[[t]][i, ] <- curr
                }
                
              }
              
              # update acceptance times after burn-in
              if (iter > Nb) {
                acc_times[t] <- acc_times[t] + 1
                mut[t] <- mut[t] + 1
              }
            }
          }
          else {
            if (iter > Nb) {
              acc_times[t] <- acc_times[t] + 1
            }
          }
        }
        else {
          JD <- 0
          acc_rate <- 1
          if (iter > Nb){
            acc_times[t] <- acc_times[t] + 1
          }
        }
        
        model_sizes[[t]][iter+1,i] <- p_gam
        log_posts[[t]][iter+1,i] <- log_post_curr
        # acc_probs[t] <- acc_probs[t] + acc_rate
        
        if (iter > Nb) {
          estm_PIPs[t,] <- estm_PIPs[t,] + curr
          if (t == 1) {
            infs[1, JD+1] <- infs[1, JD+1] + 1
            infs[2, JD+1] <- infs[2, JD+1] + acc_rate
            if (eval_f) {
              sum_f <- sum_f + f(curr)
            }
          }
          
          ESJD[t] <- ESJD[t] + acc_rate * JD
          if (t > 1) {
            JD_swap[i, t-1] <- sum(abs(curr - currs[[t-1]][i,]))
          }
        }
        
        if ((t==1) && store_chains) {
          chains[[i]][iter+1,] <- curr
        } 
        
        if (omega_adap == "kw") {
          # update ESJD to positive group or negative group
          if (omega_idx_curr) {
            ESJD_pos[t] <- ESJD_pos[t] + acc_rate * JD
          }
          else {
            ESJD_neg[t] <- ESJD_neg[t] + acc_rate * JD
          }
        }
        else if (omega_adap == "rm") {
          omega_acc_rate[t] <- omega_acc_rate[t] + acc_rate
        }
        
      }
    }
    
    # acc_prob <- acc_probs/n_chain
    
    if (n_temp > 1) {
      # update temperature adaptively
      psiter <- iter^(phi_swap)
      llhs_curr <-extract_llh(LAs, n_temp, n_chain)
      H_swap_full <- matrix(apply(exp(t((temps[1:(n_temp-1)] - temps[2:n_temp]) * t(llhs_curr[,2:n_temp] - llhs_curr[,1:(n_temp-1)]))), 
                                  MARGIN = 1:2, min, 1), nrow = n_chain, ncol = n_temp-1)
      
      if (iter > Nb) {
        ESJD_temp <- ESJD_temp + colSums(JD_swap * H_swap_full)
      }
      
      H_swap <- colMeans(H_swap_full)
      
      taus <- taus + psiter * (H_swap  - target_swap)
      temps <- generate_temperature(taus)
      
      temps_trace[iter+1,] <- temps
      taus_trace[iter+1,] <- taus
    }
    
    
    # update PIPs and omegas
    if (iter<= Nl) {
      for (t in 1:n_temp) {
        
        if (use_rb) {
          BF_temp_temp <- BF_temp[[t]]^temps[t]
          Bayes_fac[t,] <- Bayes_fac[t,] + colMeans((h_tils[[t]]*BF_temp_temp)/(1-h_tils[[t]]+h_tils[[t]]*BF_temp_temp))
          
          PIPs[t,] <- kappa + (1 - 2*kappa) * Bayes_fac[t,]/iter
        }
        else {
          sum_PIPs[t,] <- sum_PIPs[t,] + colMeans(sum_PIPs_temp[[t]])
          PIPs[t,] <- kappa + (1 - 2*kappa) * sum_PIPs[t,]/iter
        }
        
        
        A[t,] <- sapply(PIPs[t,]/(1-PIPs[t,]), min, 1)
        D[t,] <- sapply((1-PIPs[t,])/PIPs[t,], min, 1)
      }
    }
    
    if (omega_adap == "kw") {
      # update omega
      aiter <- iter^phi_a
      citer <- iter^phi_c
      omega_pre <- matrix(logit_e(omega[,1], eps) + aiter*(ESJD_pos/n_pos - ESJD_neg/n_neg)/(2*citer), nrow = n_temp, ncol = 3)
      
      citer <- (iter+1)^phi_c
      omega_curr <- inv_logit_e(omega_pre + matrix(c(0, -citer, citer), byrow = TRUE, nrow = n_temp, ncol = 3), eps)
      
      for (t in 1:n_temp) {
        
        max_inc <- 0.2
        
        if (abs(omega_curr[t,1] - omega[t,1]) <= max_inc) {
          omega[t,] <- omega_curr[t,]
        }
        else if ((omega_curr[t,1] - omega[t,1]) > max_inc){
          omega[t,] <- inv_logit_e(logit_e(omega[t,1] + max_inc, eps) + c(0, -citer, citer), eps)
        }
        else {
          omega[t,] <- inv_logit_e(logit_e(omega[t,1] - max_inc, eps) + c(0, -citer, citer), eps)
        }
        
        # omega[t,] <- omega_curr[t,]
        omegas[[t]][iter+1, ] <- omega[t,]
      }
      # omegas <- lapply(omega, margin, ...)
      # function(t, omega, omegas, iter) {omegas[[t]][iter+1, ] <- omega[t,]}
    }
    else if(omega_adap == "rm") {
      
      aiter <- iter^phi_a
      omega <- inv_logit_e(logit_e(omega, eps) + aiter*(omega_acc_rate/n_chain - target_omega), eps)
      
      for (t in 1:n_temp) {
        omegas[[t]][iter+1] <- omega[t]
      }
      
    }
    
    if (iter == Nb) {
      end.time1 <- Sys.time()
    }
    
  }
  
  end.time2 <- Sys.time()
  
  if (verbose) {
    close(pb)
  }
  
  c <- (N-Nb) * n_chain
  
  # infs
  infs[2,] <- infs[2,]/infs[1,]
  infs[1,] <- infs[1,]/c
  
  infs <- data.frame(infs, row.names = c("aver_prop_chances",
                                         "aver_acc_rates"))
  colnames(infs) <- 0:p
  infs <- infs[,-which(infs[1,]==0)]
  
  if (n_temp > 1) {
    # temp_infs
    temp_infs[2,] <- temp_infs[2,]/temp_infs[1,]
    temp_infs[1, ] <- temp_infs[1, ]/c
    
    temp_infs <- data.frame(temp_infs, row.names = c("aver_prop_swap",
                                                     "aver_swap_acc_rates"))
    colnames(temp_infs) <- 1:(n_temp-1)
  }
  
  if (use_rb) {
    sum_PIPs <- NULL
  }
  else {
    Bayes_fac <- NULL
  }
  
  
  return(list(chains = chains,
              infs = infs,
              temp_infs = temp_infs,
              log_post_trace = log_posts, 
              model_size_trace = model_sizes,
              omegas = omegas, omega = omega,
              taus_trace = taus_trace,
              temps_trace = temps_trace,
              acc_rate = acc_times/c, 
              mut_rate = mut/c,
              swap_acc_rate = swap_acc_times/c,
              estm_PIPs = estm_PIPs/c,
              rb_PIPs = Bayes_fac/Nl,
              emp_PIPs = sum_PIPs/Nl,
              ad_PIPs = PIPs,
              eval_f = sum_f/c,
              ESJD = ESJD/c,
              ESJD_swap = ESJD_swap/c,
              ESJD_temp = ESJD_temp/c,
              CPU_time = c(difftime(end.time2, start.time1, units = "mins"),
                           difftime(end.time1, start.time1, units = "mins"),
                           difftime(end.time2, start.time2, units = "mins"))
  ))
  
}





