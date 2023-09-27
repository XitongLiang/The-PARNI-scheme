reps <- 20
p <- hyper_par$p

diag_idx <- diag(matrix(1:p, p, p))

PIPs <- matrix(NA, p*(p-1), reps)
CPU_times <- rep(NA, reps)

for (i in 1:reps) {
  
  cat(i, "\n")
  
  time.start <- Sys.time()
  results_order <- orderMCMC(scorepar = score_par,
                             MAP = FALSE,
                                     startspace = postskel_iterativeMCMC$skel,
                                     plus1 = TRUE,
                                     scoreout = FALSE,
                             chainout = TRUE,
                                     # stepsave = 1,
                                     iterations = 20000*8*45,
                                     verbose = TRUE
  )
  time_order <- Sys.time() - time.start
  
  PIPs[,i] <- as.vector(edgep(results_order))[-diag_idx]
  CPU_times[i] <- time_order
  
  write_csv(as.data.frame(PIPs), "gsim100/iterative_order.csv")
  
}