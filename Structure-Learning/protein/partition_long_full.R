# protein dataset


W <- as.matrix(read.csv("protein/Data/sachs_W.csv"))[, -1]
X <- as.matrix(read.csv("protein/Data/sachs_X.csv"))[, -1]

DAG_heatmap(W)


colnames(X) <- NULL
X <- scale(X)
p <- ncol(X)




hyper_par <- list(n = nrow(X),
                  p = ncol(X),
                  g = 10,
                  h = 10/110,
                  X = t(X),
                  XtX = t(X) %*% X)




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
  data = X,
  usrpar = list(pctesttype = "bge",
                # XtX = hyper_par$XtX,
                h = hyper_par$h,
                g = hyper_par$g)
)


full_W <- matrix(1, p, p)
diag(full_W) <- 0

time_start = Sys.time()
results_partition <- partitionMCMC(scorepar = score_par,
                                  # startspace = postskel_iterativeMCMC$skel,
                                  startspace = full_W,
                                  stepsave = 1, 
                                  iterations = 5000*600, 
                                  verbose = TRUE)
time_partition <- Sys.time() - time_start  


DAG_heatmap(edgep(results_partition))

# list.save(results_partition, "protein/partition_long_full.rds")


PIPs_partition <- edgep(results_partition)
DAG_heatmap(PIPs_partition)
