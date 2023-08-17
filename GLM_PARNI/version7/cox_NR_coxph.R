library(survival)


cox_NR_coxph <- function(t, Z, X, d, p_gam, p_z, g, inv_variances, max_iter = 100){
  print(g)
  if (p_z > 0) {

    if (p_gam > 0) {

      J <- cbind(Z, X)
      full_inv_variances <- c(rep(inv_variances, p_z), rep(1/g, p_gam))
      cox_formula <- formula(Surv(t,d) ~ ridge(Z, theta = inv_variances, scale = FALSE) +
        ridge(X, theta = 1/g, scale = FALSE))

    }
    else {
      J <- Z
      full_inv_variances <- inv_variances
      cox_formula <- formula(Surv(t,d) ~ ridge(Z, theta = inv_variances, scale = FALSE))
    }


  }
  else {
    J <- X
    full_inv_variances <- 1/g
    cox_formula <- formula(Surv(t,d) ~ ridge(X, theta = 1/g, scale = FALSE))
  }

  # calucation of startval  betahat
  betahat <- rep(0, p_z + p_gam)

  # initialisation
  diff_beta <- 10
  iter <- 0

  while((diff_beta > 1e-6) & (iter <= max_iter)) {

    # print(iter)
    # print(betahat[1:10])
    cox_fit <- coxph(cox_formula, control = coxph.control(iter.max = 0),
                     init = betahat)


    inv_RXWX <- cox_fit$var

    # calculation of linear predictors and means
    eta <- as.vector(cox_fit$linear.predictors)
    # eta <- as.vector(X %*% betahat)
    lambda <- exp(eta)


    tilda_y <- cox_pseudores(lambda, t, d)
    # print(tilda_y)

    # U <- XW %*% (z - eta)

    # diag(XWX) <- diag(XWX) + 1/g
    # inv_RXWX <- solve(XWX)

    # update betahat and go back
    betahat_new <- betahat + inv_RXWX %*% (t(J) %*% tilda_y - full_inv_variances*betahat)
    print(inv_RXWX %*% t(J) %*% tilda_y )

    diff_beta <- max(abs(betahat_new - betahat))

    betahat <- betahat_new

    # update the counter
    iter <- iter + 1

  }

  # inv_Fisher <- inv_RXWX

  return(list(betahat = as.vector(betahat),
              inv_hessian = inv_RXWX,
              iter = iter,
              eta = as.vector(eta),
              J = J))

}






