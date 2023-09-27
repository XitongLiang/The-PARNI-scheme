# Pointwise implementation of Adaptive Random Neighbourhood Informed Sampler (Parallel Tempering)
PARNI <- function(alg_par, hyper_par){
  
  
  
  
  # initialisation of alg_par
  N <- alg_par$N   # number of iterations
  Nb <- alg_par$Nb     # number of burn-in interations
  n_chain <- alg_par$n_chain    # number of multiple chain
  store_chains <- alg_par$store_chains
  verbose <- alg_par$verbose
  f <- alg_par$eval_f
  omega_adap <- alg_par$omega_adap
  omega <- alg_par$omega_init
  omega_par <- alg_par$omega_par
  eps <- alg_par$eps
  kappa <- alg_par$kappa
  bal_fun <- alg_par$bal_fun
  PIPs_update <- alg_par$PIPs_update
  
  
  # proportion of using marginal pips in tuning parameters
  p_mar_PIPs <- alg_par$p_mar_PIPs
  
  if (is.null(p_mar_PIPs)) {
    p_mar_PIPs <- 0
  }
  
  use_logit_e <- alg_par$use_logit_e
  if (is.null(use_logit_e)) {
    use_logit_e <- TRUE
  }
  if (use_logit_e) {
    omega_map <- logit_e
    inv_omega_map <- inv_logit_e
  }
  else {
    omega_map <- function(x, eps) {pmax(pmin(x, 1-eps), eps)}
    inv_omega_map <- function(x, eps) {pmax(pmin(x, 1-eps), eps)}
  }
  
  
  # initialisation of hyper_para
  # n <- hyper_par$n       # number of datapoint
  p <- hyper_par$p       # number of regressors
  g <- hyper_par$g       # g
  h <- hyper_par$h       # model prior parameter
  max_p <- p*(p-1)/2
  hyper_par$max_p <- max_p
  
  use_bge <- hyper_par$use_bge
  
  if (is.null(use_bge)) {
    use_bge <- FALSE
  }
  
  
  
  
  hyper_par$permi_pars <- H_to_permi_pars(unname(alg_par$H))
  
  if (use_bge) {
    
    hyper_par$log_llh <- log_llh_bge_table
    hyper_par$log_llh_update <- log_llh_bge_update_table
    hyper_par$log_m_prior <- function(...) { 0 }
    
    
    D <- marPIPs_bge_H(hyper_par, alg_par$kappa)
    
  }
  else {
    
    # h prior type
    if (length(h) == 1) {
      # fixed h
      h_exp <- h
      hyper_par$h_odd <- function(I, h, ...) {ifelse(I, h/(1-h), (1-h)/h)}
      h_til <- function(h, ...) {h}
      hyper_par$log_m_prior <- function(p_gam, h, p) {p_gam * (log(h) - log(1-h))}
    }
    else {
      # binomial-beta
      h_alpha <- h[1]
      h_beta <- h[2]
      h_exp <- h_alpha/(h_alpha+h_beta)
      hyper_par$h_odd <- function(I, h, p_gam, p){
        ifelse(I, (p_gam+h[1])/(p-p_gam-1+h[2]), (p-p_gam+h[2])/(p_gam-1+h[1]))}
      h_til <- function(h, curr, p_gam, p) {(p_gam - curr + h[1])/(p + h[1] + h[2] - 1)}
      hyper_par$log_m_prior <- function(p_gam, h, p) {lbeta(p_gam+h[1],p-p_gam+h[2])}
    }
    
    
    # initialisation
    
    # log likelihood type
    hyper_par$log_llh <- log_llh_DAG_table
    # hyper_par$log_llh <- log_llh_DAG
    # hyper_par$XtX <- hyper_par$X %*% t(hyper_par$X)
    D <- marPIPs_DAG_H(hyper_par, alg_par$kappa)
    
    hyper_par$log_llh_update <- log_llh_DAG_update_table
    
  }
  
  if (is.null(PIPs_update)) {
    PIPs_update <- TRUE
  }
  
  
  
  
  
  # D <- alg_par$D
  hyper_par$tables <- D$tables
  approx_PIPs <- D$PIPs
  # approx_PIPs <- alg_par$approx_PIPs * (1-2*alg_par$kappa) + alg_par$kappa
  
  PIPs <- approx_PIPs
  
  # remove(D)
  # print(PIPs)
  # plot(DAG_heatmap(PIPs))
  
  A <- pmin(PIPs/(1-PIPs), 1)
  diag(A) <- 0
  
  D <- pmin((1-PIPs)/PIPs, 1)
  diag(D) <- 0
  
  
  # swap index
  # swap_idx <- matrix(1:(p^2), byrow = TRUE, p, p)
  hyper_par$swap_idx <- matrix(1:(p^2), byrow = TRUE, p, p)
  
  #  how to initalise gamma
  gamma_init <- alg_par$gamma_init
  
  if (is.null(alg_par$random_gamma_init)){
    random_gamma_init <- FALSE
  }
  else {
    random_gamma_init <- alg_par$random_gamma_init
  }
  
  
  
  
  # omegas
  if (is.null(omega)) {
    omega <- 0.9
  }
  
  if (omega_adap == "f") {
    
    omega_curr <- omega
    
  }
  else if (omega_adap == "kw") {
    
    phi_a <- omega_par[1]
    phi_c <- omega_par[2]
    
    n_pos <- floor(n_chain/2)
    n_neg <- n_chain - n_pos
    
    if (use_logit_e) {
      c_coeff <- 1
    }
    else {
      c_coeff <- 0.1
    }
    
    citer <- 1 * c_coeff
    
    omega <- inv_omega_map(omega_map(omega, eps) + c(0, -citer, citer), eps)
    omegas <- matrix(omega, byrow = TRUE, nrow = N+1, ncol = 3)
    
  }
  else if (omega_adap == "rm") {
    
    phi_a <- omega_par[1]
    # neighbourhood target size
    target_omega <- omega_par[2]
    
    omega <- omega
    omegas <- rep(omega, N+1)
    
  }
  
  
  # initial model
  # currs <- list() # matrix(runif(n = p*n_chain) < h_exp, nrow = n_chain, ncol = p)# rep(0, p)
  chains <- list()
  LAs <- list()
  
  log_posts <- matrix(NA, nrow = N+1, ncol = n_chain)
  model_sizes <- matrix(NA, nrow = N+1, ncol = n_chain)
  infs <- matrix(0, byrow = TRUE, nrow = 2, ncol = hyper_par$max_p+1)
  
  
  for (i in 1:n_chain) {
    
    if (is.null(gamma_init)) {
      
      if (random_gamma_init) {
        gamma <- diag(p)
        while(!is.DAG_adjmat(gamma)) {
          gamma <-  matrix(runif(n = p*p) < h_exp, nrow = p, ncol = p)
          diag(gamma) <- 0
        }
        curr <- gamma
      } 
      else{
        curr <- matrix(0, ncol = p, nrow = p)
      }
      
    }
    else {
      curr <- gamma_init
    }
    
    
    
    # currs[[i]] <- curr
    LAs_temp <- compute_LA_DAG(curr, hyper_par)
    LAs[[i]] <- LAs_temp
    p_gam <- LAs_temp$p_gam # sum(curr)
    model_sizes[1, i] <- p_gam
    
    log_posts[1, i] <- LAs_temp$log_post 
    
    if (store_chains) {
      # chains[[i]] <- matrix(as.vector(curr), nrow = N+1, ncol = p*p, byrow = TRUE)
      
      chains[[i]] <- list()
      chains[[i]][[1]] <- as(curr, "sparseMatrix")
      
    }
    
  }
  
  
  
  # initialisation
  acc_times <- 0
  mut <- 0
  ESJD <- 0
  acc_rate <- 0
  estm_PIPs <- matrix(0, nrow = p, ncol = p)
  k_sizes <- 0
  
  # eval_f
  eval_f <- !is.null(f)
  sum_f <- 0
  
  sum_PIPs <- matrix(0, nrow = p, ncol = p)
  
  
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
      
      ESJD_pos <- 0
      ESJD_neg <- 0
      
    }
    else if (omega_adap == "rm") {
      
      thinned_k_sizes <- 0
      
    }
    
    
    
    for (i in 1:n_chain) {
      
      # curr <- currs[[i]]
      LA <- LAs[[i]]
      curr <- LA$curr
      
      if (omega_adap == "kw") {
        omega_idx_curr <- omega_idx[i]
        omega_curr <- ifelse(omega_idx_curr, omega[3], omega[2])
      }
      else {
        omega_curr <- omega
      }
      
      
      log_post_curr <- LA$log_post 
      p_gam <- LA$p_gam
      
      # hyper_par$A <- A
      # hyper_par$D <- D
      eta <- (1-curr) * A + curr * D
      neighs <- sample_ind_DAG(TRUE, probs = eta, samples = NULL, log = TRUE)
      k <- neighs$sample
      k_size <- length(k)
      
      
      if (k_size > 0) {
        
        # if (k_size == 1) {
        #   k_ordered <- k
        # }
        # else {
        #   # k_ordered <- sample_k(k, hyper_par$swap_idx)
        #   k_ordered <- sample(k)
        #   
        # }
        
        
        updates <- update_LA_DAG(LA, k, hyper_par, bal_fun, PIPs, thinning_rate = omega_curr, omega = 1/2)
        
        
        # cat(updates$LA_prop$log_post, 
        #     hyper_par$log_llh(updates$LA_prop, hyper_par)$llh +
        #       hyper_par$log_m_prior(sum(updates$LA_prop$curr), hyper_par$h, hyper_par$max_p), "\n")
        # 
        
        
        
        JD <- updates$JD
        acc_rate <- updates$acc_rate
        thinned_k_size <- updates$thinned_k_size
        
        k_sizes <- k_sizes + thinned_k_size
        
        # return(k_ordered)
        
        if (JD > 0) {
          if (runif(1) < acc_rate){
            
            LA_prop <- updates$LA_prop
            curr <- LA_prop$curr
            log_post_curr <- LA_prop$log_post
            p_gam <- LA_prop$p_gam
            LAs[[i]] <- LA_prop
            
            
            
            
            
            # update acceptance times after burn-in
            if (iter > Nb) {
              acc_times <- acc_times + 1
              mut <- mut + 1
            }
          }
        }
        else {
          if (iter > Nb) {
            acc_times <- acc_times + 1
          }
        }
      }
      else {
        JD <- 0
        acc_rate <- 1
        thinned_k_size <- 0
        if (iter > Nb){
          acc_times <- acc_times + 1
        }
      }
      
      model_sizes[iter+1,i] <- p_gam
      log_posts[iter+1,i] <- log_post_curr
      
      
      if (iter > Nb) {
        estm_PIPs <- estm_PIPs + curr
        
        infs[1, JD+1] <- infs[1, JD+1] + 1
        infs[2, JD+1] <- infs[2, JD+1] + acc_rate
        if (eval_f) {
          sum_f <- sum_f + f(curr)
        }
        
        ESJD <- ESJD + acc_rate * JD
      }
      
      if (store_chains) {
        chains[[i]][[iter+1]] <- as(curr, "sparseMatrix")
      } 
      
      if (omega_adap == "kw") {
        # update ESJD to positive group or negative group
        if (omega_idx_curr) {
          ESJD_pos <- ESJD_pos + acc_rate * JD / (thinned_k_size + 1)
        }
        else {
          ESJD_neg <- ESJD_neg + acc_rate * JD / (thinned_k_size + 1)
        }
      }
      else if (omega_adap == "rm") {
        thinned_k_sizes <- thinned_k_sizes + thinned_k_size
      }
      
    }
    
    # acc_prob <- acc_probs/n_chain
    
    
    # update PIPs and omegas
    sum_PIPs_temp <- LAs[[1]]$curr
    if (n_chain > 1) {
      for (i in 2:n_chain) {
        sum_PIPs_temp <- sum_PIPs_temp +  LAs[[i]]$curr
      }
    }
    sum_PIPs <- sum_PIPs + sum_PIPs_temp/n_chain
    
    
    if (PIPs_update) {
      
      
      if (iter <= Nb) {
        a <- 1 - 1/(2*(Nb-iter+1)^(0.2))
      }
      else {
        a <- (iter-Nb)^(-0.5)/2
      }
      
      
      PIPs <- kappa + (1 - 2*kappa) * (a*approx_PIPs + (1-a)*sum_PIPs/iter)
      diag(PIPs) <- 0
      
      
      # compute A and D
      A <- pmin(PIPs/(1-PIPs), 1)
      D <- pmin((1-PIPs)/PIPs, 1)
      
    }
    
    
    
    
    if (omega_adap == "kw") {
      # update omega
      aiter <- iter^phi_a
      citer <- c_coeff * iter^phi_c
      omega_pre <- omega_map(omega[1], eps) + aiter*(ESJD_pos/n_pos - ESJD_neg/n_neg)/(2*citer)
      
      citer <- c_coeff * (iter+1)^phi_c
      omega_curr <- inv_omega_map(omega_pre + c(0, -citer, citer), eps)
      
      max_inc <- 0.2
      
      if (abs(omega_curr[1] - omega[1]) <= max_inc) {
        omega <- omega_curr
      }
      else if ((omega_curr[1] - omega[1]) > max_inc){
        omega <- inv_omega_map(omega_map(omega[1] + max_inc, eps) + c(0, -citer, citer), eps)
      }
      else {
        omega <- inv_omega_map(omega_map(omega[1] - max_inc, eps) + c(0, -citer, citer), eps)
      }
      
      
      omegas[iter+1, ] <- omega
      
    }
    else if(omega_adap == "rm") {
      
      aiter <- iter^phi_a
      omega <- inv_omega_map(omega_map(omega, eps) + aiter*(target_omega - thinned_k_sizes/n_chain), eps)
      # cat(thinned_k_sizes/n_chain, target_omega, omega, "\n")
      
      omega <- ifelse(omega > 0.99, 0.99,
                      ifelse(omega < 0.01, 0.01, omega))
      
      omegas[iter+1] <- omega
      
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
  colnames(infs) <- 0:hyper_par$max_p
  infs <- infs[,-which(infs[1,]==0)]
  
  
  
  return(list(chains = chains,
              infs = infs,
              log_post_trace = log_posts, 
              model_size_trace = model_sizes,
              omegas = omegas, omega = omega,
              acc_rate = acc_times/c, 
              mut_rate = mut/c,
              estm_PIPs = estm_PIPs/c,
              emp_PIPs = sum_PIPs/N,
              approx_PIPs = approx_PIPs,
              ad_PIPs = PIPs,
              eval_f = sum_f/c,
              ESJD = ESJD/c,
              k_sizes = k_sizes/N/n_chain,
              CPU_time = c(difftime(end.time2, start.time1, units = "mins"),
                           difftime(end.time1, start.time1, units = "mins"),
                           difftime(end.time2, start.time2, units = "mins"))
  ))
  
}





