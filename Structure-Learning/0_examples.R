
library(rlist)
library(tidyverse)
library(pcalg)
library(BiDAG)

setwd("./algorithms")

source("ADR.R")
source("PARNI.R")
source("Compute_LA_DAG.R")
source("sample_ind_DAG.R")
source("logit_e.R")
source("is.DAG_adjmat.R")
source("marPIPs_DAG_H.R")
source("log_llh_DAG.R")
source("log_llh_DAG_update_table.R")
source("update_LA_DAG_sawp.R")
source("other_functions.R")
source("DAG_heatmap.R")

source("marPIPs_bge_H.R")
source("log_llh_bge_table.R")

source("log_llh_bge_table.R")


# example codes for PARNI

X <- gsim100
colnames(X) <- NULL
X <- scale(X)
X <- t(X)


library(pcalg)
suffStat <- list(C = cor(t(X)), n = nrow(t(X)))
varNames <- as.character(1:100)
skel.dataset <- skeleton(suffStat, indepTest = gaussCItest, labels = varNames,
                         alpha = 0.05)
skel.W <- as(skel.dataset@graph, "matrix")

DAG_heatmap(skel.W)



p <- nrow(X)

hyper_par <- list(n = ncol(X),
                  p = nrow(X),
                  g = 10,
                  h = 1/100,
                  X = X,
                  XtX = X %*% t(X))







alg_par <- list(N = 50000, 
                Nb = 10000,
                kappa = 0.01,
                n_chain = 1,
                verbose = TRUE,
                # second value of omega_par is the expected neighoburhood size
                # which will affect computational costs and the effciencies
                omega_adap = "rm", omega_par = c(-0.7, 20), use_logit_e = TRUE,
                # omega_adap = "kw", omega_init = 0.2, omega_par = c(-1, -0.5), use_logit_e = FALSE,
                eps = 1/(hyper_par$p*(hyper_par$p-1)),
                store_chains = FALSE,
                bal_fun = function(x) {pmin(1,x)},
                H = skel.W)


results_PARNI <- PARNI(alg_par, hyper_par)


results_PARNI$omega

# results$estm_PIPs
DAG_heatmap(dataset$W)
DAG_heatmap(results_PARNI$estm_PIPs, text_bound = 0.7, true_graph = dataset$W)
DAG_heatmap(results_PARNI$ad_PIPs, text_bound = 0.7, true_graph = dataset$W)

results_PARNI$ESJD
results_PARNI$acc_rate
results_PARNI$CPU_time
results_PARNI$infs
results_PARNI$k_sizes

plot(results_PARNI$log_post_trace[,1], type = "l")







# example for partitionMCMC using posterior distribution specified in the paper




usrscoreparameters <- function (initparam, usrpar = list(pctesttype = "usrCItest")) {
  
  p <- ncol(initparam$data)
  n <- nrow(initparam$data)
  
  if (is.null(usrpar$g)) {
    initparam$g <- max(n, p)
  }
  else {
    initparam$g <- usrpar$g
  }
  
  
  if (is.null(usrpar$h)) {
    initparam$h <- 10/(p*(p-1))
  }
  else {
    initparam$h <- usrpar$h
  }
  
  
  initparam$h
  
  X <- initparam$data
  initparam$XtX <- t(X) %*% X
  
  initparam$n_observations <- n
  
  initparam
}













usrDAGcorescore <- function(j,parentnodes,n,param) {
  
  
  h <- param$h
  g <- param$g
  n_observations <- param$n_observations
  XtX <- param$XtX
  
  p_pa <- length(parentnodes)
  
  
  if (p_pa == 0) {
    
    log_llh <- - n_observations*log(XtX[j,j]/2)/2
    
    
  }
  else {
    
    x_paj_xj_t <- XtX[parentnodes, j]
    
    if (p_pa == 1) {
      
      L_Vg <- sqrt(XtX[parentnodes, parentnodes] + 1/g)
      
      
      log_llh <- - log(L_Vg) - log(g)/2 - n_observations * log((XtX[j,j] -  (x_paj_xj_t/L_Vg)^2)/2)/2
      
    }
    else {
      
      diag_p_paj <- diag(p_pa)
      
      Vg <- XtX[parentnodes, parentnodes] + diag_p_paj/g
      # diag(Vg) <- diag(Vg) + 1/g
      L_Vg <- chol(Vg)
      
      log_llh <- - sum(log(diag(L_Vg))) - p_pa*log(g)/2 - n_observations * log((XtX[j,j] - sum((forwardsolve(t(L_Vg), diag_p_paj) %*%
                                                                                                  x_paj_xj_t)^2))/2)/2
      
    }
    
    
  }
  
  lmp <- p_pa * (log(h) - log(1-h))
  
  log_post <- log_llh + lmp
  
  
  
  return(log_post)
  
  
}




rlang::env_unlock(env = asNamespace('BiDAG'))
rlang::env_binding_unlock(env = asNamespace('BiDAG'))

assign('usrDAGcorescore', usrDAGcorescore, envir = asNamespace('BiDAG'))
assign('usrscoreparameters', usrscoreparameters, envir = asNamespace('BiDAG'))

rlang::env_binding_lock(env = asNamespace('BiDAG'))
rlang::env_lock(asNamespace('BiDAG'))





score_par <- scoreparameters(
  scoretype = c("usr"),
  data = t(hyper_par$X),
  usrpar = list(pctesttype = "bge",
                # XtX = hyper_par$XtX,
                h = hyper_par$h,
                g = hyper_par$g)
)




time_start = Sys.time()
results_partition <- partitionMCMC(scorepar = score_par,
                                   startspace = skel.W,
                                   scoreout = FALSE,
                                   stepsave = 1,
                                   iterations = 50000,
                                   verbose = TRUE
)

time_partition <- Sys.time() - time_start

