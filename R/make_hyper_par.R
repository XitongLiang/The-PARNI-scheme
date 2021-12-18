make_hyper_par <- function(y, X, g, h, Z = NULL, prior_type = 1){
  
  
  X <- as.matrix(X)
  y <- as.vector(y)
  colnames(X) <- NULL
  colnames(y) <- NULL
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(Z)) {
    X <- scale(X)
    y <- scale(y, scale = FALSE)
  }
  else {
    # remove the fixed linear effects
    Z <- as.matrix(cbind(1, Z))
    y <- (diag(n) - Z %*% solve(t(Z) %*% Z) %*% t(Z)) %*% y
    X <- (diag(n) - Z %*% solve(t(Z) %*% Z) %*% t(Z)) %*% X
  }
  
  
  
  if (prior_type == 1) {
    g_prior_type <- "ind"
    diag_V <- colSums(X^2) + 1/g
  }
  else if (prior_type == 2) {
    g_prior_type <- "g"
    diag_V <- colSums(X^2)
  }
  
  hyper_par <- list(n = n,
                    p = p,
                    g = g,
                    h = h,
                    X = X,
                    y = y,
                    yty = sum(y^2), 
                    ytX = t(y) %*% X,
                    diag_V = diag_V,
                    g_prior_type = g_prior_type)
  
  return(hyper_par)
}
