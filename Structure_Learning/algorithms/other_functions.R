H_to_permi_pars <- function(H) {
  
  p <- nrow(H)
  permi_pars <- list()
  
  
  for (j in 1:p) {
    
    permi_pars[[j]] <- which(H[,j] == 1)
    
  }
  
  return(permi_pars)
  
}




model_encoding <- function(curr_part) {
  
  p <- length(curr_part)
  m <- sum(curr_part * 2^(0:(p-1))) + 1
  
  return(m)
  
}





model_decoding <- function(p, m) {
  
  curr_part <- convert_to_binary(m-1)
  
  curr_part <- rev(c(rep(0, p-length(curr_part)), curr_part))
  
  return(curr_part)
  
}




