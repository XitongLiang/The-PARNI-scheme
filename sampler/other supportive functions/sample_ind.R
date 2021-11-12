
# sampling independently
sample_ind <- function(whe_sam, probs, n_sam = NULL, samples = NULL, log = FALSE){
  
  n <- length(probs)
  
  if (whe_sam) {
    samples <- runif(n) < probs
  }
  
  if (log) {
    prob <- sum(log(probs[samples]))
  }
  else {
    prob <- prod(probs[samples])
  }
  
  return(list(prob = prob,
              sample = samples))
  
}
