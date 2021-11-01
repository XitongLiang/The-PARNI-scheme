# Aadd-delete-swap with multiple chain and parallel tempering
ADS <- function(alg_par, hyper_par){
  
  # initialisation of alg_par
  N <- alg_par$N   # number of iterations
  Nb <- alg_par$Nb     # number of burn-in interations
  n_chain <- alg_par$n_chain    # number of multiple chain
  n_temp <- alg_par$n_temp
  store_chains <- alg_par$store_chains
  verbose <- alg_par$verbose
  f <- alg_par$f
  taus <- alg_par$init_taus
  phi_swap <- alg_par$phi
  target_swap <- alg_par$target_swap
  
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
    # h_til <- function(h, ...) {h}
    log_m_prior <- function(p_gam, h, p) {p_gam * (log(h) - log(1-h))}
  }
  else {
    # binomial-beta
    h_alpha <- h[1]
    h_beta <- h[2]
    h_exp <- h_alpha/(h_alpha+h_beta)
    hyper_par$h_odd <- function(I, h, p_gam, p){
      ifelse(I, (p_gam+h[1])/(p-p_gam-1+h[2]), (p-p_gam+h[2])/(p_gam-1+h[1]))}
    # h_til <- function(h, curr, p_gam) {(p_gam - curr + h[1])/(p + h[1] + h[2] - 1)}
    log_m_prior <- function(p_gam, h, p) {lbeta(p_gam+h[1],p-p_gam+h[2])}
  }
  
  
  # g prior type
  if (hyper_par$g_prior_type == "g") {
    log_llh <- log_llh_g_L
  }
  else if (hyper_par$g_prior_type == "ind") {
    log_llh <- log_llh_ind_L
  }
  
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

  # swap_acc_times
  swap_acc_times <- 0
  
  # initial model
  currs <- list() # matrix(runif(n = p*n_chain) < h_exp, nrow = n_chain, ncol = p)# rep(0, p)
  chains <- list()
  llhs_curr <- matrix(NA, nrow = n_chain, ncol = n_temp)
  prior_curr <- matrix(NA, nrow = n_chain, ncol = n_temp)
  log_posts <- list() # matrix(NA, nrow = N+1, ncol = n_chain)
  model_sizes <- list() # matrix(NA, nrow = N+1, ncol = n_chain)
  infs <- matrix(0, byrow = TRUE, nrow = 2, ncol = p+1)
  
  for (t in 1:n_temp) {
    currs_temp <- matrix(runif(n = p*n_chain) < h_exp, nrow = n_chain, ncol = p)
    log_posts_temp <- matrix(NA, nrow = N+1, ncol = n_chain)
    model_sizes_temp <- matrix(NA, nrow = N+1, ncol = n_chain)
    
    for (i in 1:n_chain) {
      curr <- currs_temp[i,]
      p_gam <- sum(curr)
      model_sizes_temp[1, i] <- p_gam
      llhs_curr[i, t] <- log_llh(curr, hyper_par) 
      prior_curr[i, t] <- log_m_prior(p_gam, h, p)
      log_posts_temp[1, i] <- temps[t] * llhs_curr[i, t] + prior_curr[i, t]
      
      if ((t == 1) && store_chains) {
        chains[[i]] <- matrix(curr, nrow = N+1, ncol = p, byrow = TRUE)
      }
    }
    
    currs[[t]] <- currs_temp
    log_posts[[t]] <- log_posts_temp
    model_sizes[[t]] <- model_sizes_temp
    
  }
  
  
  # initialisation
  acc_times <- rep(0, n_temp)
  ESJD <- rep(0, n_temp)
  ESJD_swap <- 0
  ESJD_temp <- rep(0, n_temp-1)
  acc_rate <- 0
  estm_PIPs <- matrix(0, nrow = n_temp, ncol = p)
  JD_swap <- matrix(0, nrow = n_chain, ncol = n_temp-1)
  
  eval_f <- ifelse(is.null(f), FALSE, TRUE)
  sum_f <- 0
  
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
    
    for (i in 1:n_chain) {
      
      if (n_temp > 1) {
        # swap move
        swap_idx <- sample(1:(n_temp-1), size = 1)
        
        acc_rate_swap <- min(1, exp((temps[swap_idx] - temps[swap_idx+1]) * (llhs_curr[i,swap_idx+1] - llhs_curr[i,swap_idx])))
        # print(acc_rate_swap)
        if (runif(1) < acc_rate_swap) {
          
          curr_temp <- currs[[swap_idx]][i,]
          currs[[swap_idx]][i,] <- currs[[swap_idx+1]][i,]
          currs[[swap_idx+1]][i,] <- curr_temp
          
          llhs_curr_temp <- llhs_curr[i,swap_idx]
          llhs_curr[i,swap_idx] <- llhs_curr[i,swap_idx+1]
          llhs_curr[i,swap_idx+1] <- llhs_curr_temp
          
          prior_curr_temp <- prior_curr[i,swap_idx]
          prior_curr[i,swap_idx] <- prior_curr[i,swap_idx+1]
          prior_curr[i,swap_idx+1] <- prior_curr_temp
          
          
          
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
        p_gam <- sum(curr) # model_sizes[[t]][iter, i]
        
        if ((p_gam > 1) && (p_gam < (p-1))) {
          
          incl <- which(curr == 1)
          excl <- which(curr == 0)
          
          if (runif(1) < 1/3) {
            # swap
            change <- c(sample(x = incl, size = 1), sample(x = excl, size = 1))
            log_prop <- 0
            log_reverse <- 0
          }
          else if (runif(1) < 1/2){
            # add
            change <- sample(x = excl, size = 1)
            log_prop <- -log(length(excl))
            log_reverse <- -log(length(incl)+1)
          }
          else {
            # delete
            change <- sample(x = incl, size = 1)
            log_prop <- -log(length(incl))
            log_reverse <- -log(length(excl)+1)
          }
        }
        else if (p_gam == 1) {
          
          incl <- which(curr == 1)
          excl <- which(curr == 0)
          
          if (runif(1) < 1/2) {
            # add
            change <- sample(x = excl, size = 1)
            log_prop <- -log(length(excl))
            log_reverse <- -log(length(incl)+1)
          }
          else {
            # swap
            change <- c(incl, sample(x = excl, size = 1))
            log_prop <- 0
            log_reverse <- 0
          }
        }
        else if ((p_gam == 0) || (p_gam == p)){
          change <- sample(x = p, size = 1)
          log_prop <- -log(p)
          log_reverse <- 0
        }
        else if (p_gam == (p-1)) {
          
          incl <- which(curr == 1)
          excl <- which(curr == 0)
          
          if (runif(1) < 1/2) {
            # delete
            change <- sample(x = incl, size = 1)
            log_prop <- -log(length(incl))
            log_reverse <- -log(length(excl)+1)
          }
          else {
            # swap
            change <- c(excl, sample(x = incl, size = 1))
            log_prop <- 0
            log_reverse <- 0
          }
        }
        
        
        prop <- curr
        prop[change] <- 1 - prop[change]
        log_llh_prop <- log_llh(prop, hyper_par)
        p_gam_prop <- sum(prop)
        log_m_prior_prop <- log_m_prior(p_gam_prop, h, p)
        log_pi_prop <- temps[t]*log_llh_prop + log_m_prior_prop
        log_pi_curr <- temps[t]*llhs_curr[i, t] + prior_curr[i, t]
        
        JD <- length(change)
        log_ratio <- log_pi_prop + log_reverse - log_pi_curr - log_prop
        acc_rate <- min(1, exp(log_ratio))
        
        # cat(t, i, temps[t], "\n")
        
        if (runif(1) < acc_rate){
          
          curr <- prop
          currs[[t]][i,] <- curr
          p_gam <- p_gam_prop
          log_pi_curr <- log_pi_prop
          
          llhs_curr[i, t] <- log_llh_prop
          prior_curr[i, t] <- log_m_prior_prop
          #if(t==1){cat(log_pi_prop == log_post_ind_L(prop, hyper_par))}
          #if(t==1){cat(log_pi_curr == log_post_ind_L(curr, hyper_par),"\n")}
          # update acceptance times after burn-in
          if (iter > Nb) {
            acc_times[t] <- acc_times[t] + 1
          }
        }
        
        model_sizes[[t]][iter+1,i] <- p_gam
        log_posts[[t]][iter+1,i] <- log_pi_curr
        
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
      }
    }
    
   
    
    if (iter == Nb) {
      end.time1 <- Sys.time()
    }
    
    if (n_temp > 1) {
      # update temperature adaptively
      psiter <- iter^(phi_swap)
      
      H_swap_full <- matrix(apply(exp(t((temps[1:(n_temp-1)] - temps[2:n_temp]) * t(llhs_curr[,2:n_temp] - llhs_curr[,1:(n_temp-1)]))), 
                            MARGIN = 1:2, min, 1), nrow = n_chain, ncol = n_temp-1)
      
      
      if (iter > Nb) {
        ESJD_temp <- ESJD_temp + colSums(JD_swap * H_swap_full)
      }
      
      H_swap <- colMeans(H_swap_full)
      
      taus <- taus + psiter * (H_swap - target_swap)
      
      temps <- generate_temperature(taus)
      
      temps_trace[iter+1,] <- temps
      taus_trace[iter+1,] <- taus
    }

  }
  
  if (verbose) {
    close(pb)
  }
  
  end.time2 <- Sys.time()
  
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
  
  
  return(list(chains = chains,
              infs = infs,
              temp_infs = temp_infs,
              log_post_trace = log_posts, 
              model_size_trace = model_sizes,
              taus_trace = taus_trace,
              temps_trace = temps_trace,
              acc_rate = acc_times/c, 
              swap_acc_rate = swap_acc_times/c,
              estm_PIPs = estm_PIPs/c,
              f = sum_f/c,
              ESJD = ESJD/c,
              ESJD_swap = ESJD_swap/c,
              ESJD_temp = ESJD_temp/c,
              CPU_time = c(difftime(end.time2, start.time1, units = "mins"),
                           difftime(end.time1, start.time1, units = "mins"),
                           difftime(end.time2, start.time2, units = "mins"))
  ))
  
}




