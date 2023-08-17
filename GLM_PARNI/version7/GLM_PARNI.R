# Pointwise implementation of Adaptive Random Neighbourhood Informed Sampler (Parallel Tempering)
GLM_PARNI <- function(alg_par, hyper_par){
  
  # initialisation of alg_par
  N <- alg_par$N   # number of iterations
  Nb <- alg_par$Nb     # number of burn-in interations
  # Nl <- alg_par$Nl
  # Nrb <- alg_par$Nrb
  n_chain <- alg_par$n_chain    # number of multiple chain
  store_chains <- alg_par$store_chains
  verbose <- alg_par$verbose
  f <- alg_par$eval_f
  
  
  # full_adap <- alg_par$full_adap
  kappa <- alg_par$kappa
  bal_fun <- alg_par$bal_fun
  # pg_cpp <- alg_par$pg_cpp
  
  
  
  
  # initialisation of hyper_para
  n <- hyper_par$n       # number of datapoint
  p <- hyper_par$p       # number of regressors
  g <- hyper_par$g       # g
  h <- hyper_par$h       # model prior parameter
  hyper_par <- make_fixed_variances(hyper_par)
  
  
  # adapting eta
  eta <- rep(0, n)
  # eta <- alg_par$eta
  sum_eta <- rep(0, n)
  whether_eta <- TRUE
  
  
  
  # h prior type
  if (length(h) == 1) {
    # fixed h
    h_exp <- h
    hyper_par$h_exp <- h_exp
    hyper_par$h_odd <- function(I, h, ...) {ifelse(I, h/(1-h), (1-h)/h)}
    # h_til <- function(h, ...) {h}
    log_m_prior <- function(p_gam, h, p) {p_gam * (log(h) - log(1-h))}
  }
  else {
    # binomial-beta
    h_alpha <- h[1]
    h_beta <- h[2]
    h_exp <- h_alpha/(h_alpha+h_beta)
    hyper_par$h_exp <- h_exp
    hyper_par$h_odd <- function(I, h, p_gam, p){
      ifelse(I, (p_gam+h[1])/(p-p_gam-1+h[2]), (p-p_gam+h[2])/(p_gam-1+h[1]))}
    # h_til <- function(h, curr, p_gam) {(p_gam - curr + h[1])/(p + h[1] + h[2] - 1)}
    log_m_prior <- function(p_gam, h, p) {lbeta(p_gam+h[1],p-p_gam+h[2])}
  }
  
  hyper_par$log_m_prior <- log_m_prior
  
  
  model <- hyper_par$model
  method <- alg_par$method
  
  
  if (model == "logistic") {
    
    
    if (method == "DA") {
      
      pg_cpp <- alg_par$pg_cpp
      
      if (is.null(pg_cpp)) {
        pg_cpp <- TRUE
      } 
      
      if (pg_cpp) {
        update_omega <- logistic_update_omega_cpp
      }
      else {
        update_omega <- logistic_update_omega_r
      }
      
      hyper_par$kappa <- hyper_par$y - 1/2
      log_post <- logistic_update_matrices
      
      PIPs_rb <- logistic_DA_rb
      hyper_par$approx_log_post <- logistic_update_matrices
      
      
    }
    else {
      
      PIPs_rb <- logistic_ALA_rb
      
      if (method == "LA") {
        
        log_post <- logistic_laplace
        
      }
      else if (method == "CPM") {
        
        hyper_par$n_particles <- alg_par$n_particles # number of particles used in importance estimation
        hyper_par$rho <- alg_par$paricle_corr
        
        log_post <- logistic_cpm
        # hyper_par$log_fixed_const <- sum(log(choose(1, hyper_par$y)))
        hyper_par$log_fixed_const <- 1
      }
      
      
      use_ALA2 <- ifelse(is.null(alg_par$use_ALA2), TRUE, alg_par$use_ALA2)
      use_ALA <- ifelse(is.null(alg_par$use_ALA), FALSE, alg_par$use_ALA)
      if (use_ALA2) {
        hyper_par$approx_log_post <- logistic_ALA2
      }
      else if (use_ALA) {
        hyper_par$approx_log_post <- logistic_ALA
      }
      else if (TRUE) {
        hyper_par$approx_log_post <- logistic_laplace
      }
      
    }
    
    
    
    
  }
  else if (model == "Cox") {
    
    hyper_par$tilde_y <- cox_pseudores(rep(1, n), hyper_par$t, hyper_par$d)
    hyper_par$W <- cox_compute_ALA_W(hyper_par$t, hyper_par$d)
    PIPs_rb <- cox_ALA_rb
    
    # hyper_par$sqrt_W <- sqrt(diag(cox_compute_ALA_W(hyper_par$t, hyper_par$d)))
    # # # print(hyper_par$sqrt_W)
    # PIPs_rb <- cox_ALA_rb2
    
    if (method == "LA") {
      
      log_post <- cox_laplace
      
      
    }
    else if (method == "CPM") {
      
      # number of particles used in importance estimation
      hyper_par$n_particles <- alg_par$n_particles 
      hyper_par$rho <- alg_par$paricle_corr
      
      log_post <- cox_cpm
      hyper_par$log_fixed_const <- 1
      
    }
    
    use_ALA2 <- ifelse(is.null(alg_par$use_ALA2), TRUE, alg_par$use_ALA2)
    use_ALA <- ifelse(is.null(alg_par$use_ALA), FALSE, alg_par$use_ALA)
    if (use_ALA2) {
      hyper_par$approx_log_post <- cox_ALA2
    }
    else if (use_ALA) {
      hyper_par$approx_log_post <- cox_ALA
    }
    else if (TRUE) {
      hyper_par$approx_log_post <- cox_laplace
    }
    
    
  }
  else if (model == "Weibull") {
    
    k_rw_var <- 1
    ks <- matrix(0, nrow = N, ncol = n_chain)
    k_acc_times <- 0
    
    PIPs_rb <- weibull_ALA_rb
    
    if (method == "LA") {
      
      log_post <- weibull_laplace
      
    }
    else if (method == "CPM") {
      
      # number of particles used in importance estimation
      hyper_par$n_particles <- alg_par$n_particles 
      hyper_par$rho <- alg_par$paricle_corr
      
      log_post <- weibull_cpm
      hyper_par$log_fixed_const <- 1
      
    }
    
    use_ALA2 <- ifelse(is.null(alg_par$use_ALA2), TRUE, alg_par$use_ALA2)
    use_ALA <- ifelse(is.null(alg_par$use_ALA), FALSE, alg_par$use_ALA)
    if (use_ALA2) {
      hyper_par$approx_log_post <- weibull_ALA2
    }
    else if (use_ALA) {
      hyper_par$approx_log_post <- weibull_ALA
    }
    else if (TRUE) {
      hyper_par$approx_log_post <- weibull_laplace
    }
    
    
  }
  
  hyper_par$log_post <- log_post
  
  
  # initialisation
  
  if (model == "logistic" && method == "DA") {
    
    approx_PIPs <- rep(h_exp, p)
    sum_approx_PIPs <- rep(0, p)
    counts_approx <- 0
    
  }
  else {
    
    g_temp <- ifelse(is.numeric(g), g, 
                     ifelse(is.null(alg_par$g_init), 1, alg_par$g_init))
    
    approx_PIPs <- PIPs_rb(list(curr = rep(0, p), p_gam = 0, k = 1, g = g_temp), hyper_par)
    sum_approx_PIPs <- approx_PIPs
    counts_approx <- 1
    
  }
  
  
  # barplot(approx_PIPs[1:30], ylim = c(0,1))
  
  PIPs <- kappa + (1-2*kappa) * approx_PIPs
  # PIPs <- kappa + (1-2*kappa) * alg_par$PIPs
  
  PIPs_update <- TRUE
  update_approx <- TRUE
  iter_approx <- alg_par$iter_approx
  
  
  if (is.null(iter_approx)) {
    iter_approx <- 5
  }
  
  A <- sapply(PIPs/(1-PIPs), min, 1)
  D <- sapply((1-PIPs)/PIPs, min, 1)
  
  
  
  
  
  
  # omegas
  
  omega_adap <- alg_par$omega_adap
  omega <- alg_par$omega_init
  omega_par <- alg_par$omega_par
  eps <- alg_par$eps
  
  if (is.null(omega)) {
    omega <- 0.5
  }
  
  if (omega_adap == "f") {
    
    omega_curr <- omega
    omegas <- omega
    
  }
  else if (omega_adap == "kw") {
    
    phi_a <- omega_par[1]
    phi_c <- omega_par[2]
    
    n_pos <- floor(n_chain/2)
    n_neg <- n_chain - n_pos
    citer <- 1
    
    omega <- inv_logit_e(logit_e(omega, eps) + c(0, -citer, citer), eps)
    
    omegas <- matrix(omega[1,], byrow = TRUE, nrow = N+1, ncol = 3)
    
  }
  else if (omega_adap == "rm") {
    
    phi_a <- omega_par[1]
    target_omega <- omega_par[2]
    
    omegas <- rep(omega[1], N+1)
    
  }
  
  
  
  
  # zeta
  zeta_adap <- alg_par$zeta_adap
  zeta <- alg_par$zeta_init
  zeta_par <- alg_par$zeta_par
  
  if (is.null(zeta_adap)) {
    zeta_adap <- FALSE
  }
  
  if (zeta_adap) {
    
    phi_a <- zeta_par[1]
    target_k_size <- zeta_par[2]
    
    if (is.null(zeta)) {
      
      # zeta <- 0.95 # 1 - 2*eps
      
      init_sizes <- sum(A + h_exp*(1-2*PIPs)/PIPs)
      zeta <- target_k_size/init_sizes 
      zeta <- ifelse(zeta > (1 - 2*eps), 1 - 2*eps, 
                     ifelse(zeta < (2*eps), 2*eps, zeta))
      
    }
    
    # print(zeta)
    zetas <- rep(zeta, N+1)
    
  }
  else {
    
    if (is.null(zeta)) {
      zeta <- 1
    }
    zetas <- zeta
    
  }
  
  
  sum_k_sizes <- 0
  
  
  
  
  # if g is random
  if (is.null(g)) {
    g <- "random"
  }
  
  if (g == "random") {
    
    update_g <- TRUE
    
    if (is.null(alg_par$g_init)) {
      g_curr <- rhcauchy(n_chain, 0, 1)
    }
    else {
      g_curr <- rep(alg_par$g_init, n_chain)
    }
    
    gs <- matrix(0, nrow = N, ncol = n_chain)
    g_rw_var <- 1
    g_acc_times <- 0
    
    
  }
  else {
    update_g <- FALSE
    g_curr <- rep(g, n_chain)
  }
  
  
  
  
  
  
  # initial model
  ALL <- list()
  chains <- list()
  log_posts <- list() # matrix(NA, nrow = N+1, ncol = n_chain)
  model_sizes <- list() # matrix(NA, nrow = N+1, ncol = n_chain)
  infs <- matrix(0, byrow = TRUE, nrow = 2, ncol = p+1)
  
  
  ALL_temp <- list()
  
  if (is.null(alg_par$gamma_init)) {
    currs_temp <- matrix(runif(n = p*n_chain) < h_exp, nrow = n_chain, ncol = p)
  }
  else {
    currs_temp <- matrix(alg_par$gamma_init, byrow = TRUE, nrow = n_chain, ncol = p)
  }
  
  log_posts_temp <- matrix(NA, nrow = N+1, ncol = n_chain)
  model_sizes_temp <- matrix(NA, nrow = N+1, ncol = n_chain)
  
  
  for (i in 1:n_chain) {
    
    ALL_temp_0 <- list()
    
    curr <- currs_temp[i,]
    p_gam <- sum(curr)
    
    ALL_temp_0$curr <- curr
    ALL_temp_0$p_gam <- p_gam
    ALL_temp_0$g <- g_curr[i]
    
    if (model == "logistic" && method == "DA") {
      # generate the first set of omega
      
      ALL_temp_0 <- update_omega(ALL_temp_0, hyper_par, initn = TRUE)
      
    }
    else if (model == "Weibull") {
      # initialise k
      ALL_temp_0$k <- 1
    }
    
    ALL_temp_0 <- log_post(ALL_temp_0, hyper_par)
    
    
    # store in a large list
    ALL_temp[[i]] <- ALL_temp_0
    
    log_posts_temp[1, i] <- ALL_temp_0$llh + ALL_temp_0$lmp
    model_sizes_temp[1, i] <- p_gam
    
    if (store_chains) {
      chains[[i]] <- matrix(curr, nrow = N+1, ncol = p, byrow = TRUE)
    }
  }
  
  # store in a large list
  ALL <- ALL_temp
  log_posts <- log_posts_temp
  model_sizes <- model_sizes_temp
  
  
  
  # initialisation
  acc_times <- 0
  mut <- 0
  ESJD <- 0
  
  acc_rate <- 0
  estm_PIPs <- 0
  
  # eval_f
  eval_f <- !is.null(f)
  sum_f <- 0
  
  # h_tils <- matrix(NA, nrow = n_chain, ncol = p)
  # BF_temp <- matrix(NA, nrow = n_chain, ncol = p)
  # Bayes_fac <- matrix(0, nrow = n_temp, ncol = p)
  # Nrb_iters <- 0
  
  sum_PIPs_temp <- matrix(NA, nrow = n_chain, ncol = p)
  sum_PIPs <- rep(0, p)
  # S <- list()
  
  
  start.time1 <- Sys.time()
  
  if ((verbose)) {
    pb <- txtProgressBar(min = 0, max = N, style = 3)
  }
  
  
  # return(list(ALL,hyper_par))
  
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
      
      omega_acc_rate <- 0
      
    }
    
    thinned_k_sizes <- 0
    
    # barplot(approx_PIPs, ylim = c(0,1))
    
    
    for (i in 1:n_chain) {
      
      
      ALL_curr <- ALL[[i]]
      
      curr <- ALL_curr$curr
      p_gam <- ALL_curr$p_gam
      
      if (omega_adap == "kw") {
        omega_idx_curr <- omega_idx[i]
        omega_curr <- ifelse(omega_idx_curr, omega[3], omega[2])
      }
      else if (omega_adap == "rm") {
        omega_curr <- omega
      }
      
      log_pi_curr <- ALL_curr$llh + ALL_curr$lmp
      
      
      prob_AD <- abs(curr-1) * A + curr * D
      neighs <- sample_ind(TRUE, n_sam = NULL, probs = prob_AD, samples = NULL, log = TRUE)
      k <- which(neighs$sample)
      k_size <- length(k)
      
      if (k_size > 0) {
        
        if (k_size == 1) {
          k_ordered <- k
        }
        else {
          k_ordered <- sample(k)
        }
        
        
        proposals <- PARNI_proposal(ALL_curr, k_ordered, hyper_par, alg_par, PIPs, omega_curr, zeta, eta)
        
        # prop <- proposals$prop
        JD <- proposals$JD
        acc_rate <- proposals$acc_rate
        thinned_k_sizes <- thinned_k_sizes + proposals$thinned_k_size
        sum_k_sizes <- sum_k_sizes + proposals$thinned_k_size
        
        if (JD > 0) {
          if (runif(1) < acc_rate){
            
            ALL_curr <- proposals$ALL_prop
            curr <- ALL_curr$curr
            log_pi_curr <- ALL_curr$llh + ALL_curr$lmp
            p_gam <- ALL_curr$p_gam
            
            
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
        if (iter > Nb){
          acc_times <- acc_times + 1
        }
      }
      
      model_sizes[iter+1,i] <- p_gam
      log_posts[iter+1,i] <- log_pi_curr
      
      
      sum_eta <- sum_eta + ALL_curr$eta
      
      
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
        chains[[i]][iter+1,] <- curr
      } 
      
      if (omega_adap == "kw") {
        # update ESJD to positive group or negative group
        if (omega_idx_curr) {
          ESJD_pos <- ESJD_pos + acc_rate * JD
        }
        else {
          ESJD_neg <- ESJD_neg + acc_rate * JD
        }
      }
      else if (omega_adap == "rm") {
        omega_acc_rate <- omega_acc_rate + acc_rate
      }
      
      
      
      
      if (model == "logistic" && method == "DA") {
        # generate the first set of omega
        
        ALL_curr <- update_omega(ALL_curr, hyper_par, initn = FALSE)
        ALL_curr <- log_post(ALL_curr, hyper_par)
        
      }
      else if (model == "Weibull") {
        # update k
        
        k_lists <- weibull_update_k(ALL_curr, hyper_par, k_rw_var, iter, n_chain)
        ALL_curr <- k_lists$ALL
        ks[iter, i] <- k_lists$ALL$k
        
        k_rw_var <- k_lists$k_rw_var
        k_acc_times <- k_acc_times + k_lists$k_acc
        
      }
      
      if (update_g) {
        
        # update_g
        g_lists <- hc_update_g(ALL_curr, hyper_par, g_rw_var, iter, n_chain)
        
        ALL_curr <- g_lists$ALL
        gs[iter, i] <- g_lists$ALL$g
        
        g_rw_var <- g_lists$g_rw_var
        g_acc_times <- g_acc_times + g_lists$g_acc
        
      }
      
      ALL[[i]] <- ALL_curr
      
      
    }
    
    # acc_prob <- acc_probs/n_chain
    
    
    
    
    # update PIPs and omegas
    if (PIPs_update) {
      
      if (update_approx) {
        
        counts_approx <- counts_approx + 1
        
        for (i in 1:n_chain) {
          # sum_approx_PIPs <- sum_approx_PIPs + logstic_advanced_ALA_rb(ALL[[i]], hyper_par, eta)/n_chain
          sum_approx_PIPs <- sum_approx_PIPs + PIPs_rb(ALL[[i]], hyper_par)/n_chain
        }
        
        approx_PIPs <- sum_approx_PIPs/counts_approx
        
        # barplot(approx_PIPs[1:50], ylim = c(0,1))
        
        if (iter >= iter_approx) {
          
          update_approx <- FALSE
          
        }
        
      }
      
      
      for (i in 1:n_chain) {
        sum_PIPs_temp[i,] <- ALL[[i]]$curr
      }
      sum_PIPs <- sum_PIPs + colMeans(sum_PIPs_temp)
      
      if (iter <= Nb) {
        a <- 1 - 1/(2*(Nb-iter+1)^(0.2))
      }
      else {
        a <- (iter-Nb)^(-0.5)/2
      }
      
      
      PIPs <- kappa + (1 - 2*kappa) * (a*approx_PIPs + (1-a)*sum_PIPs/iter)
      
      
      A <- sapply(PIPs/(1-PIPs), min, 1)
      D <- sapply((1-PIPs)/PIPs, min, 1)
      
    }
    
    
    # if (whether_eta && (iter >= Nb)) {
    if (whether_eta) {
      
      eta <- sum_eta/(iter*n_chain)
      
    }
    
    
    if (zeta_adap) {
      
      aiter <- iter^phi_a
      
      zeta <- inv_logit_e(logit_e(zeta, eps) - aiter*(thinned_k_sizes/n_chain - target_k_size), eps)
      
      zetas[iter + 1] <- zeta
      
    }
    
    
    # barplot(PIPs[1,], ylim = c(0,1))
    
    if (omega_adap == "kw") {
      # update omega
      aiter <- iter^phi_a
      citer <- iter^phi_c
      omega_pre <- logit_e(omega[,1], eps) + aiter*(ESJD_pos/n_pos - ESJD_neg/n_neg)/(2*citer)
      
      citer <- (iter+1)^phi_c
      omega_curr <- inv_logit_e(omega_pre + c(0, -citer, citer), eps)
      
      max_inc <- 0.2
      
      if (abs(omega_curr[1] - omega[1]) <= max_inc) {
        omega <- omega_curr
      }
      else if ((omega_curr[1] - omega[1]) > max_inc){
        omega <- inv_logit_e(logit_e(omega + max_inc, eps) + c(0, -citer, citer), eps)
      }
      else {
        omega <- inv_logit_e(logit_e(omega[1] - max_inc, eps) + c(0, -citer, citer), eps)
      }
      
      omegas[iter+1, ] <- omega
      
      
    }
    else if(omega_adap == "rm") {
      
      aiter <- iter^phi_a
      omega <- inv_logit_e(logit_e(omega, eps) + aiter*(omega_acc_rate/n_chain - target_omega), eps)
      
      omegas[iter+1] <- omega
      
    }
    # plot(PIPs[1,], ylim = c(0,1))
    
    if (iter == Nb) {
      end.time1 <- Sys.time()
      
      # print(sum(S_neigh_counts)/p)
      
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
  
  if (model != "Weibull") {
    ks <- NULL
    k_acc_times <- NULL
    k_rw_var <- NULL
  }
  
  if (!update_g) {
    gs <- NULL
    g_acc_times <- NULL
    g_rw_var <- NULL
  }
  
  
  return(list(chains = chains,
              infs = infs,
              log_post_trace = log_posts, 
              model_size_trace = model_sizes,
              omegas = omegas, omega = omega,
              zetas = zetas, zeta = zeta, thinned_k_sizes = sum_k_sizes/(N*n_chain),
              acc_rate = acc_times/c, 
              mut_rate = mut/c,
              estm_PIPs = estm_PIPs/c,
              # rb_PIPs = Bayes_fac,
              mar_PIPs = sum_PIPs/N,
              approx_PIPs = approx_PIPs,
              ad_PIPs = PIPs,
              eval_f = sum_f/c,
              ESJD = ESJD/c,
              weibull_ks = ks,
              weibull_k_acc_times = k_acc_times/(n_chain*N),
              weibull_k_rw_var = k_rw_var,
              random_gs = gs,
              random_g_acc_times = g_acc_times/(n_chain*N),
              random_g_rw_var = g_rw_var,
              CPU_time = c(difftime(end.time2, start.time1, units = "mins"),
                           difftime(end.time1, start.time1, units = "mins"),
                           difftime(end.time2, start.time2, units = "mins")),
              eta = eta
  ))
  
}





