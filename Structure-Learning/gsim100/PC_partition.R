reps <- 20
p <- hyper_par$p

diag_idx <- diag(matrix(1:p, p, p))

PIPs <- matrix(NA, p*(p-1), reps)
CPU_times <- rep(NA, reps)

for (i in 1:reps) {
  
  cat(i, "\n")
  
  time.start <- Sys.time()
  results_partition <- partitionMCMC(scorepar = score_par,
                                     startspace = skel.W,
                                     scoreout = FALSE,
                                     stepsave = 1,
                                     iterations = 20000*8,
                                     verbose = TRUE
  )
  time_partition <- Sys.time() - time.start
  
  PIPs[,i] <- as.vector(edgep(results_partition))[-diag_idx]
  CPU_times[i] <- time_partition
  
  write_csv(as.data.frame(PIPs), "gsim100/PC_partition.csv")
  
}