
# sampling independently
sample_ind_DAG <- function(whe_sam, probs, samples = NULL, log = FALSE){
  
  d <- nrow(probs)
  
  if (whe_sam) {
    samples <- which(matrix(runif(d*d), d, d) < probs)
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