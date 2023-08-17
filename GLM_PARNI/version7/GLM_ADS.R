# Aadd-delete-swap with multiple chain and parallel tempering
# Laplace approximation
GLM_ADS <- function(alg_par, hyper_par){
  
  # initialisation of alg_par
  N <- alg_par$N   # number of iterations
  Nb <- alg_par$Nb     # number of burn-in interations
  n_chain <- alg_par$n_chain    # number of multiple chain
  store_chains <- alg_par$store_chains
  verbose <- alg_par$verbose
  f <- alg_par$f
  
  
  # initialisation of hyper_para
  n <- hyper_par$n       # number of datapoint
  p_z <- hyper_par$p_z     # number of fixed regressors
  p <- hyper_par$p       # number of regressors
  g <- hyper_par$g       # g
  h <- hyper_par$h       # model prior parameter
  # hyper_par$kappa <- hyper_par$y - hyper_par$n_trials/2
  # X <- hyper_par$X
  # Z <- hyper_par$Z
  hyper_par <- make_fixed_variances(hyper_par)
  
  
  # h prior type
  if (length(h) == 1) {
    # fixed h
    h_exp <- h
    # hyper_par$h_odd <- function(I, h, ...) {ifelse(I, h/(1-h), (1-h)/h)}
    # h_til <- function(h, ...) {h}
    hyper_par$log_m_prior <- function(p_gam, h, p) {p_gam * (log(h) - log(1-h))}
  }
  else {
    # binomial-beta
    h_alpha <- h[1]
    h_beta <- h[2]
    h_exp <- h_alpha/(h_alpha+h_beta)
    # hyper_par$h_odd <- function(I, h, p_gam, p){
    #   ifelse(I, (p_gam+h[1])/(p-p_gam-1+h[2]), (p-p_gam+h[2])/(p_gam-1+h[1]))}
    # h_til <- function(h, curr, p_gam) {(p_gam - curr + h[1])/(p + h[1] + h[2] - 1)}
    hyper_par$log_m_prior <- function(p_gam, h, p) {lbeta(p_gam+h[1],p-p_gam+h[2])}
  }
  
  model <- hyper_par$model
  method <- alg_par$method
  
  if (model == "logistic") {
    
    if (method == "LA") {
      
      log_post <- logistic_laplace
      
    }
    else if (method == "DA") {
      
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
      
      
    }
    else if (method == "CPM") {
      
      hyper_par$n_particles <- alg_par$n_particles # number of particles used in importance estimation
      hyper_par$rho <- alg_par$paricle_corr
      
      log_post <- logistic_cpm
      # hyper_par$log_fixed_const <- sum(log(choose(1, hyper_par$y)))
      hyper_par$log_fixed_const <- 1
      
    }
    
  }
  else if (model == "Cox") {
    
    
    
    if (method == "LA") {
      
      log_post <- cox_laplace
      
    }
    else if (method == "CPM") {
      
      hyper_par$n_particles <- alg_par$n_particles # number of particles used in importance estimation
      hyper_par$rho <- alg_par$paricle_corr
      
      log_post <- cox_cpm
      hyper_par$log_fixed_const <- 1
      
    }
    
  }
  else if (model == "Weibull") {
    
    k_rw_var <- 1
    k_acc_times <- 0
    ks <- matrix(0, nrow = N, ncol = n_chain)
    
    if (method == "LA") {
      
      log_post <- weibull_laplace
      
    }
    else if (method == "CPM") {
      
      hyper_par$n_particles <- alg_par$n_particles 
      hyper_par$rho <- alg_par$paricle_corr
      
      log_post <- weibull_cpm
      # hyper_par$log_fixed_const <- sum(log(choose(1, hyper_par$y)))
      hyper_par$log_fixed_const <- 1
      
    }
    
  }
  
  hyper_par$log_post <- log_post
  
  
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
  
  ALL <- list()
  
  log_posts <- matrix(NA, nrow = N+1, ncol = n_chain)
  model_sizes <- matrix(NA, nrow = N+1, ncol = n_chain)
  
  for (i in 1:n_chain) {
    
    ALL_temp_0 <- list()
    
    curr <- runif(n = p) < h_exp
    p_gam <- sum(curr)
    ALL_temp_0$curr <- curr
    ALL_temp_0$p_gam <- p_gam
    ALL_temp_0$g <- g_curr[i]
    
    if (model == "logistic" && method == "DA") {
      # generate the first set of omega
      ALL_temp_0 <- update_omega(ALL_temp_0, hyper_par, initn = TRUE)
    }
    
    if (model == "Weibull") {
      # initialise k
      ALL_temp_0$k <- 1
    }
    
    ALL_temp_0 <- log_post(ALL_temp_0, hyper_par)
    
    
    # store in a large list
    ALL[[i]] <- ALL_temp_0
    
    log_posts[1, i] <- ALL_temp_0$llh + ALL_temp_0$lmp
    model_sizes[1, i] <- p_gam
    
    if (store_chains) {
      chains[[i]] <- matrix(curr, nrow = N+1, ncol = p, byrow = TRUE)
    }
  }
  
    
  
  
  # initialisation
  acc_times <- 0
  ESJD <- 0
  acc_rate <- 0
  estm_PIPs <- rep(0, p)
  
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
      
      ALL_curr <- ALL[[i]]
      
      curr <- ALL_curr$curr
      p_gam <- ALL_curr$p_gam
      
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
      ALL_prop <- ALL_curr
      ALL_prop$curr <- prop
      ALL_prop$p_gam <- sum(prop)
      
      ALL_prop <- log_post(ALL_prop, hyper_par)
      
      
      log_pi_prop <- ALL_prop$llh + ALL_prop$lmp
      log_pi_curr <- ALL_curr$llh + ALL_curr$lmp
      
      
      JD <- length(change)
      log_ratio <- log_pi_prop + log_reverse - log_pi_curr - log_prop
      acc_rate <- min(1, exp(log_ratio))
      
      
      
      if (runif(1) < acc_rate){
        
        curr <- prop
        ALL_curr <- ALL_prop
        log_pi_curr <- log_pi_prop
        
        if (iter > Nb) {
          acc_times <- acc_times + 1
        }
      }
      
      model_sizes[iter+1,i] <- ALL_curr$p_gam
      log_posts[iter+1,i] <- log_pi_curr
      
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
      
      
      if (model == "logistic" && method == "DA") {
        
        # update polya-gamma random variables
        ALL_curr <- update_omega(ALL_curr, hyper_par, initn = FALSE)
        
        # update matrices again
        ALL_curr <- log_post(ALL_curr, hyper_par)
        
        # ALL[[i]] <- ALL_curr
        
      }
      else if (model == "Weibull") {
        
        # update k
        k_lists <- weibull_update_k(ALL_curr, hyper_par, k_rw_var, iter, n_chain)
        ALL_curr <- k_lists$ALL
        ks[iter, i] <- k_lists$ALL$k
        # print(k_lists$ALL$k)
        
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
    
    
    
    if (iter == Nb) {
      end.time1 <- Sys.time()
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
              acc_rate = acc_times/c, 
              estm_PIPs = estm_PIPs/c,
              f = sum_f/c,
              ESJD = ESJD/c,
              weibull_ks = ks,
              weibull_k_acc_times = k_acc_times/(n_chain*N),
              weibull_k_rw_var = k_rw_var,
              random_gs = gs,
              random_g_acc_times = g_acc_times/(n_chain*N),
              random_g_rw_var = g_rw_var,
              CPU_time = c(difftime(end.time2, start.time1, units = "mins"),
                           difftime(end.time1, start.time1, units = "mins"),
                           difftime(end.time2, start.time2, units = "mins"))
  ))
  
}



