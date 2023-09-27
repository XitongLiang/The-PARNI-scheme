is.DAG_adjmat <- function(W) {
  
  num_precedent <- colSums(W)
  
  is.DAG <- TRUE
  
  d <- nrow(W)
  
  while (d > 1) {
    
    # print(num_precedent)
    
    if (all(num_precedent > 0)) {
      is.DAG <- FALSE
      break
    }
    
    leaf <- which(num_precedent == 0)
    W <- W[-leaf,-leaf]
    d <- d - length(leaf)
    # print(d)
    if (d > 1){
      num_precedent <- colSums(W)
    }
    
  }
  
  return(is.DAG)
  
}