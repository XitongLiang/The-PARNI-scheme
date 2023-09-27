# traceplot

library(scales)
library(ggthemes)

full_W <- matrix(1, p, p)
diag(full_W) <- 0



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












# ADR

alg_par <- list(N = 80000*2*3, 
                Nb = 16000*2*3,
                n_chain = 1,
                verbose = TRUE,
                store_chains = FALSE,
                H = full_W)


results_ADR <- ADR(alg_par, hyper_par)
plot(results_ADR$log_post_trace, type = "l")
max(results_ADR$log_post_trace)

results_ADR$CPU_time

trace_ADR <- results_ADR$log_post_trace[c(1, seq(8*3, 80000*2*3, by = 8*3)),1]    


# PARNI

alg_par <- list(N = 20000*3, 
                Nb = 4000*3,
                kappa = 0.01,
                n_chain = 1,
                verbose = TRUE,
                PIPs_update = TRUE,
                # second value of omega_par is the expected neighoburhood size
                # which will affect computational costs and the effciencies
                # omega_adap = "f", omega_init = 0.1,
                omega_adap = "rm", omega_par = c(-0.7, 10), use_logit_e = TRUE,
                # omega_adap = "kw", omega_init = 0.2, omega_par = c(-1, -0.5), use_logit_e = FALSE,
                eps = 1/(hyper_par$p*(hyper_par$p-1)),
                store_chains = FALSE,
                bal_fun = function(x) {pmin(1,x)},
                H = full_W)


results_PARNI <- PARNI(alg_par, hyper_par)
plot(results_PARNI$log_post_trace, type = "l")
max(results_PARNI$log_post_trace)

results_PARNI$CPU_time

# partitionMCMC
trace_PARNI <- results_PARNI$log_post_trace[c(1, seq(3, 20000*3, by = 3)),1]


time.start <- Sys.time()
results_partition <- partitionMCMC(scorepar = score_par,
                                   startspace = full_W,
                                   scoreout = TRUE,
                                   stepsave = 1,
                                   iterations = 20000,
                                   verbose = TRUE
)
time_partition <- Sys.time() - time.start



plot(results_partition$trace, type = "l")

trace_partition <- results_partition$trace




plot(trace_PARNI, type = "l", col = "red")
lines(trace_ADR, type = "l", col = "blue")
lines(trace_partition, type = "l", col = "green")


traces <- data.frame(trace_PARNI = trace_PARNI,
                     trace_ADR = trace_ADR, 
                     trace_partition = trace_partition)




# write_csv(traces, "protein/traces.csv")

traces <- read_csv("protein/traces.csv")

traces %>%
  mutate(iterations = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("trace_"),
               names_to = c("Sampler"),
               names_pattern = "trace_(.*)",
               values_to = "log_pmp") %>%
  mutate(Sampler = ifelse(Sampler == "PARNI", "PARNI-DAG",
                          ifelse (Sampler == "partition", "partitionMCMC", Sampler))) %>%
  ggplot(aes(x = iterations, y = log_pmp)) +
  geom_line() +
  facet_grid(~ Sampler) + 
  ylab("log posterior model probability") + 
  theme_igray() +
  ylim(-25760, -25720) +
  scale_colour_colorblind() # + theme(legend.position="bottom")



traces %>%
  mutate(iterations = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("trace_"),
               names_to = c("Sampler"),
               names_pattern = "trace_(.*)",
               values_to = "log_pmp") %>%
  mutate(Sampler = ifelse(Sampler == "PARNI", "PARNI-DAG",
                          ifelse (Sampler == "partition", "partitionMCMC", Sampler))) %>%
  ggplot(aes(x = iterations, y = log_pmp)) +
  geom_line() +
  facet_wrap(~Sampler, nrow = 3) + 
  ylab("log posterior model probability") + 
  theme_minimal() +
  ylim(-25760, -25720) +
  scale_colour_colorblind() # + theme(legend.position="bottom")






traces %>%
  mutate(iterations = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("trace_"),
               names_to = c("Sampler"),
               names_pattern = "trace_(.*)",
               values_to = "log_pmp") %>%
  mutate(Sampler = ifelse(Sampler == "PARNI", "PARNI-DAG",
                          ifelse (Sampler == "partition", "partitionMCMC", Sampler))) %>%
  ggplot(aes(x = iterations, y = log_pmp)) +
  geom_line() +
  facet_wrap(~Sampler, nrow = 3) + 
  ylab("log posterior model probability") + 
  theme_minimal() +
  # ylim(-25760, -25720) +
  scale_colour_colorblind() 



