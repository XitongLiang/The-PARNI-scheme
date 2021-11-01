extract_llh <- function(LAs, n_temp, n_chain) {
  llhs <- matrix(NA, nrow = n_chain, ncol = n_temp)
  
  for (i in 1:n_chain) {
    for (t in 1:n_temp) {
      llhs[i, t] <- LAs[[t]][[i]]$llh
    }
  }
  return(llhs)
}

