library(tidyverse)

ns <- c(500, 1000)
ps <- c(500, 5000, 50000)
SNRs <- c(0.5, 1, 2, 3)

source('simulation/linear regression random sample generator.R')

settings <- list(n = ns,
                 p = ps,
                 SNR = SNRs) %>% 
  cross_d() %>%
  mutate(name = paste(n, p, SNR, sep = "_"),
         file.name = paste("simulation/dataset/", paste(n, p, SNR, sep = "_"), ".csv", sep = ""))

write_csv(settings, "simulation/simulation.csv")



pb <- txtProgressBar(min = 0, max = nrow(settings), style = 3)

for (i in 1:nrow(settings)) {
  n <- settings$n[i]
  p <- settings$p[i]
  SNR <- settings$SNR[i]
  
  dataset <- lrrsg(n = n, p = p, rho = 0.6,SNR = SNR, seed = NULL)
  X <- scale(dataset$X)
  y <- scale(dataset$y, scale = FALSE)
  
  dataset <- as_tibble(cbind(y, X))
  
  write_csv(dataset, file = settings$file.name[i])
  
  setTxtProgressBar(pb, i)
}

close(pb)





