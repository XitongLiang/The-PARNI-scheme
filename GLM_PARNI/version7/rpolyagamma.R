# polya-gamma random variable generator

# random sampling from PG(n, z) where n is an integer and z is a positive number

rpolyagamma <- function(n, z, t = pi/2) {
  
  if (n == 1) {
    
    z <- abs(z)/2
    K <- pi^2/8 + z^2/2
    
    # normalising constants
    p <- 2*exp(-z)*pigaussian(t, mu = 1/z)
    q <- pi*exp(-K*t)/(2*K)
    
    # rejection algorithm
    while(TRUE) {
      
      # two pronged approach
      if ((runif(1)*(p+q)) < p) {
        # truncated inverse gaussian
        X <- rigaussian(mu = 1/z, t=t)
      }
      else {
        # truncated exponential
        X <- rexp(1, rate = 1)
        X <- t + X/K
      }
      
      S <- a_sequence(0, X, z, t)
      n <- 0
      V <- runif(1, min = 0, max = S)
      
      while (TRUE) {
        
        n <- n + 1
        
        if ((n %% 2) == 0) {
          
          # n is even
          S <- S + a_sequence(n, X, z, t)
          
          if (V > S) {
            break
          }
          
        }
        else {
          
          # n is odd
          S <- S - a_sequence(n, X, z, t)
          
          if (V < S) {
            return(X/4)
          }
          
        }
      }
    }
  }
  else {
    ransams <- rep(NA, n)
    for (i in 1:n) {
      ransams[i] <- rployagamma(1, z, t = t)
    }
    return(sum(ransams))
  }
  
  
}


# cdf of inverse gaussion random variable with parameter mu and lambda = 1
pigaussian <- function(t, mu) {
  
  P <- pnorm(sqrt(1/t)*(t/mu-1)) + exp(2/mu)*pnorm(-sqrt(1/t)*(t/mu+1))
  
  return(P)
}


# random variable generator for truncated (at t) inverse gaussian distribution 
# with parameter mu and lambda = 1
rigaussian <- function(mu, t = 2/pi) {
  
  success <- TRUE
  
  if (mu > t) {
    
    while (success) {
      
      success.2 <- TRUE
      
      while (success.2) {
        E <- rexp(n=1)
        E. <- rexp(n=1)
        if (E^2 <= (2*E./t)) {
          success.2 <- FALSE
        }
      }
      
      X <- t/(1+t*E)^2
      alpha <- exp(-mu^(-2)*X/2)
      
      if (runif(1) <= alpha) {
        success <- FALSE
      }
    }
  }
  else {
    while (success) {
      
      Y <- rchisq(n = 1, df = 1)
      X <- mu + 0.5*mu^2*Y-0.5*mu*sqrt(4*mu*Y+(mu*Y)^2)
      
      if (((mu + X)*runif(1)) > mu) {
        X <- mu^2/X
      }
      
      if (X <= t) {
        success <- FALSE
      }
    }
  }
  
  return(X)
}


a_sequence <- function(n, X, z, t = 2/pi) {
  
  if (X <= t) {
    # L
    a0 <- pi * (n+1/2) * (2/(pi*X))^{3/2} * exp(-2*(n+1/2)^2/X)
    a <- cosh(z)*exp(-X*z^2/2) * a0
  }
  else {
    # R
    a0 <- pi * (n+1/2) * exp(-(n+1/2)^2*pi^2*X/2)
    a <- cosh(z)*exp(-X*z^2/2) * a0
  }
  return(a)
  
}


