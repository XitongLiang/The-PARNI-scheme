reps <- 3

N <- c(6315165,6802485,2613720,2450730,480000,540990,
       7249485,5706825,2299260,2658495,480000,524235,
       4290600,2912670,2286885,1798260,519675,483975,
       3980985,2912670,2431710,1757760,450000,483975)

# N <- as.integer(N/30)
Nb <- as.integer(N/3)


library(tidyverse)

settings <- read_csv("simulation/simulation.csv", col_types = cols())

sampler_name <- "ADS"

for (i in 1:24) {
  
  print(i)
  
  file.name <- settings$file.name[i]
  PIPs_fn <- paste("simulation/estm_PIPs/", sampler_name, "_", settings$name[i], ".csv", sep = "")
  
  
  dataset <- read_csv(file.name, col_types = cols(), progress = FALSE)
  y <- pull(dataset %>% select(V1))
  X <- as.matrix(dataset %>% select(-V1))
  colnames(X) <- NULL
  
  g <- 9
  h <- 10/settings$p[i]
  n <- settings$n[i]
  p <- settings$p[i]
  SNR <- settings$SNR[i]
  remove(dataset)
  
  hyper_par <- list(n = n,
                    p = p,
                    g = g,
                    h = h,
                    X = X,
                    y = y,
                    # XtX = t(X) %*% X,
                    yty = sum(y^2), 
                    ytX = t(y) %*% X,
                    g_prior_type = "ind")
  
  remove(y, X)
  
  ADS_PIPs <- matrix(NA, nrow = p, ncol = reps)
  
  alg_par <- list(N = N[i],
                  Nb = Nb[i],
                  verbose = TRUE,
                  n_chain = 1,
                  store_chains = FALSE,
                  phi = -0.7,
                  # init_taus = rep(-1,3),
                  target_swap = 0.237,
                  n_temp = 1
  )
  
  for (rep in 1:reps) {
    
    
    results <- mul_ADS_PT(alg_par, hyper_par)
    ADS_PIPs[,rep] <- results$estm_PIPs
    print(results$CPU_time[1])
    
  }
  
  write_csv(data.frame(ADS_PIPs), PIPs_fn)
  
}











