reps <- 20
p <- hyper_par$p

diag_idx <- diag(matrix(1:p, p, p))

PIPs <- matrix(NA, p*(p-1), reps)
CPU_times <- rep(NA, reps)

for (i in 1:reps) {
  
  cat(i, "\n")
  
  alg_par <- list(N = 80000*11, 
                  Nb = 16000*11,
                  n_chain = 1,
                  verbose = TRUE,
                  store_chains = FALSE,
                  H = postskel_iterativeMCMC$skel)
  
  
  results_ADR <- ADR(alg_par, hyper_par)
  
  PIPs[,i] <- as.vector(results_ADR$estm_PIPs)[-diag_idx]
  CPU_times[i] <- results_ADR$CPU_time[1]
  
  write_csv(as.data.frame(PIPs), "gsim100/iterative_ADR.csv")
  
}
