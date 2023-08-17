logit <- function(x) {
  log(x)-log(1-x)
}


inv_logit <- function(x) {
  1/(1+exp(-x))
}