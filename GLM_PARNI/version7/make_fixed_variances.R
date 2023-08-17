make_fixed_variances <- function(hyper_par) {
  
  var_intercept <- hyper_par$var_intercept
  var_fixed_covariate <- hyper_par$var_fixed_covariate
  
  p_z <- hyper_par$p_z
  
  if (hyper_par$model == "Cox") {
    if (!is.null(p_z)) {
      
      if (p_z > 0) {
        # hyper_par$inv_variances <- rep(1/var_fixed_covariate, p_z)
        hyper_par$inv_variances <- 1/var_fixed_covariate
      }
      
    }
    else {
      hyper_par$p_z <- 0
    }
  }
  else {
    if (!is.null(p_z)) {
      inv_variances <- rep(0, p_z)
      
      if (!is.null(var_intercept)) {
        inv_variances[1] <- 1/var_intercept
      }
      
      if (hyper_par$p_z == 1) {
        hyper_par$var_fixed_covariate <- NULL
      }
      else if (!is.null(var_fixed_covariate)) {
        inv_variances[-1] <- 1/var_fixed_covariate
      }
      
      hyper_par$inv_variances <- inv_variances
    }
  }
  
  
  
  
  
  return(hyper_par)
  
}