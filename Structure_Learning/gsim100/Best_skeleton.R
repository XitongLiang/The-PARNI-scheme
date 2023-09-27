# using iterativeMCMC to find the optimal space for gsim100

X <- gsim100
colnames(X) <- NULL
X <- scale(X)
X <- t(X)

p <- nrow(X)

hyper_par <- list(n = ncol(X),
                  p = nrow(X),
                  g = 10,
                  h = 1/100,
                  X = X,
                  XtX = X %*% t(X))



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
results_iterative = iterativeMCMC(scorepar = score_par, 
                          MAP = TRUE, 
                          cpdag = FALSE,
                          scoreout = TRUE,
                          # stepsave = 1,
                          # iterations = 50000,
                          accum = TRUE,
                          verbose = TRUE, 
                          chainout = TRUE)

time_order = Sys.time() - time_start



pip_order = as.matrix(edgep(results_iterative, pdag = FALSE, burnin = 0.17))
DAG_heatmap(pip_order, text_bound = 0.7, true_graph = dataset$W)

plot(results_iterative$trace[[1]], type = "l")


DAG_heatmap(results_iterative$startspace, text_bound = 0.7, true_graph = dataset$W)
DAG_heatmap(results_iterative$endspace, text_bound = 0.7, true_graph = dataset$W)



postskel_iterativeMCMC <- list(skel = results_iterative$endspace)
list.save(postskel_iterativeMCMC, "gsim100/gsim100_postskel_iterativeMCMC.rds")


results_iterative$score

postskel_iterativeMCMC <- list.load("gsim100/gsim100_postskel_iterativeMCMC.rds")



