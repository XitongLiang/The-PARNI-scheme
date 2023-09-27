reps <- 20
p <- hyper_par$p

diag_idx <- diag(matrix(1:p, p, p))

PIPs <- matrix(NA, p*(p-1), reps)
CPU_times <- rep(NA, reps)

for (i in 1:reps) {
  
  cat(i, "\n")
  
  alg_par <- list(N = 187500, 
                  Nb = 37500,
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
                  H = skel.W)
  
  
  results_PARNI <- PARNI(alg_par, hyper_par)
  
  PIPs[,i] <- as.vector(results_PARNI$estm_PIPs)[-diag_idx]
  CPU_times[i] <- results_PARNI$CPU_time[1]
  
  write_csv(as.data.frame(PIPs), "gsim100/PC_PARNI.csv")
  
}
