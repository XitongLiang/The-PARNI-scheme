# Aadd-delete-reverse with multiple chain and parallel tempering
ADR <- function(alg_par, hyper_par){
  
  # initialisation of alg_par
  N <- alg_par$N   # number of iterations
  Nb <- alg_par$Nb     # number of burn-in interations
  n_chain <- alg_par$n_chain    # number of multiple chain
  store_chains <- alg_par$store_chains
  verbose <- alg_par$verbose
  f <- alg_par$f
  
  
  #  how to initalise gamma
  if (is.null(alg_par$randon_gamma_init)){
    randon_gamma_init <- FALSE
  }
  else {
    randon_gamma_init <- alg_par$randon_gamma_init
  }
  
  
  # initialisation of hyper_para
  # n <- hyper_par$n       # number of datapoint
  p <- hyper_par$p       # number of regressors
  g <- hyper_par$g       # g
  h <- hyper_par$h       # model prior parameter
  max_p <- p*(p-1)/2
  
  # h prior type
  if (length(h) == 1) {
    # fixed h
    h_exp <- h
    hyper_par$h_odd <- function(I, h, ...) {ifelse(I, h/(1-h), (1-h)/h)}
    # h_til <- function(h, ...) {h}
    hyper_par$log_m_prior <- function(p_gam, h, p) {p_gam * (log(h) - log(1-h))}
  }
  else {
    # binomial-beta
    h_alpha <- h[1]
    h_beta <- h[2]
    h_exp <- h_alpha/(h_alpha+h_beta)
    hyper_par$h_odd <- function(I, h, p_gam, p){
      ifelse(I, (p_gam+h[1])/(p-p_gam-1+h[2]), (p-p_gam+h[2])/(p_gam-1+h[1]))}
    # h_til <- function(h, curr, p_gam) {(p_gam - curr + h[1])/(p + h[1] + h[2] - 1)}
    hyper_par$log_m_prior <- function(p_gam, h, p) {lbeta(p_gam+h[1],p-p_gam+h[2])}
  }
  
  # log likelihood type
  hyper_par$log_llh <- log_llh_DAG_table
  # hyper_par$XtX <- hyper_par$X %*% t(hyper_par$X)
  
  hyper_par$permi_pars <- H_to_permi_pars(unname(alg_par$H))
  hyper_par$tables <- marPIPs_DAG_H(hyper_par)$tables
    
  # initial model
  # currs <- list() # matrix(runif(n = p*n_chain) < h_exp, nrow = n_chain, ncol = p)# rep(0, p)
  LAs <- list()
  chains <- list()
  # llhs_curr <- rep(NA, n_chain)
  # prior_curr <- rep(NA, n_chain)
  log_posts <- matrix(NA, nrow = N+1, ncol = n_chain)
  model_sizes <- matrix(NA, nrow = N+1, ncol = n_chain)
  infs <- matrix(0, byrow = TRUE, nrow = 2, ncol = max_p)
  
  for (i in 1:n_chain) {
    
    if (randon_gamma_init) {
      gamma <- diag(p)
      while(is.DAG_adjmat(gamma)) {
        gamma <-  matrix(runif(n = p*p) < h_exp, nrow = p, ncol = p)
      }
      curr <- gamma
    } 
    else{
      curr <- matrix(0, ncol = p, nrow = p)
    }
    
    LAs_temp <- compute_LA_DAG(curr, hyper_par)
    LAs[[i]] <- LAs_temp
    
    # p_gam <- sum(curr)
    model_sizes[1, i] <- LAs_temp$p_gam
    log_posts[1, i] <- LAs_temp$log_post 
    
    # p_gam <- sum(curr)
    # model_sizes[1, i] <- p_gam
    # llhs_curr[i] <- log_llh(curr, hyper_par) 
    # prior_curr[i] <- log_m_prior(p_gam, h, max_p)
    # log_posts[1, i] <- llhs_curr[i] + prior_curr[i]
    
    if (store_chains) {
      # chains[[i]] <- matrix(as.vector(curr), nrow = N+1, ncol = p*p, byrow = TRUE)
      chains[[i]] <- list()
      chains[[i]][[1]] <- as(curr, "sparseMatrix")
    }
    
    # currs[[i]] <- curr
    
  }
  
  # transpose index for swap move
  transpose_index <- matrix(1:(p^2), nrow = p, ncol = p, byrow = TRUE)
  
  
  # initialisation
  acc_times <- 0
  ESJD <- 0
  acc_rate <- 0
  estm_PIPs <- matrix(0, nrow = p, ncol = p)
  propose_DAG <- 0
  
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
      
      LA <- LAs[[i]]
      
      curr <- LA$curr
      p_gam <- LA$p_gam
      
      if ((p_gam >= 1) && (p_gam <= (max_p-1))) {
        
        if (runif(1) < 1/3) {
          
          # reverse
          edges <- which(curr==1)
          
          if (length(edges) > 1) {
            change <- sample(edges, size = 1)
          }
          else {
            change <- edges
          }
          
          
          change <- c(change, transpose_index[change])
          
          log_prop <- 0
          log_reverse <- 0
          
        }
        else if (runif(1) < 1/2) {
          
          # add
          possible_curr <- curr + t(curr)
          diag(possible_curr) <- 1
          
          empty_edges <- which(possible_curr == 0) 
          
          if (length(empty_edges) > 1) {
            change <- sample(empty_edges, size = 1)
          }
          else {
            change <- empty_edges
          }
          
          log_prop <- - log(length(empty_edges))
          log_reverse <- - log(p_gam + 1)
          
          # if (change %in% diag(transpose_index)) {print(change)}
          
        }
        else {
          
          # delete
          connected_edges <- which(curr == 1) 
          
          if (length(connected_edges) > 1) {
            change <- sample(connected_edges, size = 1)
          }
          else {
            change <- connected_edges
          }
          
          log_prop <- - log(length(connected_edges))
          log_reverse <- - log(2) - log(max_p - (length(connected_edges)-1))
          
        }
        
      }
      else if (p_gam < 1) {
        
        # add
        possible_curr <- diag(p)
        
        empty_edges <- which(possible_curr == 0) 
        change <- sample(empty_edges, size = 1)
        
        log_prop <- - log(length(empty_edges))
        log_reverse <- - log(3)
        
      }
      else if (p_gam > (max_p-1)) {
        
        # delete
        connected_edges <- which(curr == 1) 
        change <- sample(connected_edges, size = 1)
        
        log_prop <- - log(length(connected_edges))
        log_reverse <- - log(3) - log(2)
        
      }
      
      
      prop <- curr
      prop[change] <- 1 - prop[change]
      
      
      if (is.DAG_adjmat(prop)) {
        
        propose_DAG <- propose_DAG + 1
        
        LA_prop <- compute_LA_DAG(prop, hyper_par)
        
        # cat(log_llh_prop, log_llh_2_new(prop, hyper_par), "\n")
        
        # p_gam_prop <- sum(prop)
        # log_m_prior_prop <- log_m_prior(p_gam_prop, h, max_p)
        # log_pi_prop <- log_llh_prop + log_m_prior_prop
        # log_pi_curr <- llhs_curr[i] + prior_curr[i]
        
        p_gam_prop <- LA_prop$p_gam
        log_pi_prop <- LA_prop$log_post
        log_pi_curr <- LA$log_post
        
        JD <- length(change)
        log_ratio <- log_pi_prop + log_reverse - log_pi_curr - log_prop
        acc_rate <- min(1, exp(log_ratio))
        
        
        if (runif(1) < acc_rate){
          
          curr <- prop
          LAs[[i]] <- LA_prop
          p_gam <- p_gam_prop
          log_pi_curr <- log_pi_prop
          
          # llhs_curr[i] <- log_llh_prop
          # prior_curr[i] <- log_m_prior_prop
          
          # update acceptance times after burn-in
          if (iter > Nb) {
            acc_times <- acc_times + 1
          }
        }
      }
      
      model_sizes[iter+1,i] <- p_gam
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
        chains[[i]][[iter+1]] <- as(curr, "sparseMatrix")
      } 
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
  
  
  return(list(chains = chains,
              infs = infs,
              log_post_trace = log_posts, 
              model_size_trace = model_sizes,
              acc_rate = acc_times/c, 
              propose_DAG_rate = propose_DAG/(N*n_chain),
              estm_PIPs = estm_PIPs/c,
              f = sum_f/c,
              ESJD = ESJD/c,
              CPU_time = c(difftime(end.time2, start.time1, units = "mins"),
                           difftime(end.time1, start.time1, units = "mins"),
                           difftime(end.time2, start.time2, units = "mins"))
  ))
  
}




