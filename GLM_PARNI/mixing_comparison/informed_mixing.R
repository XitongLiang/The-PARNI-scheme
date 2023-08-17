# logistic regression

# simulate a dataset
set.seed(8533)

p <- 500
n <- 500
SNR <- 1 # sqrt(log(p)/n)*10
rho <- 0.6


X <- matrix(rnorm(n*p), ncol = p, nrow = n)

for (j in 2:p) {
  X[,j] <- rho*X[,j-1] + sqrt(1-rho^2)*X[,j]
}

beta <- SNR*c(2,-3, 2, 2,-3, 3,-2, 3,-2, 3, rep(0, p-10))

X <- scale(X)
eta <- X %*% beta
mu <- exp(eta)/(1+exp(eta))

y <- as.vector(runif(n) < mu)
mean(y)

Z <- as.matrix(rep(1, n), nrow = n)

h <- 10/p
g <- 1

hyper_par <- list(n = n,
                  p = p,
                  p_z = 1,
                  # g = g,
                  var_intercept = 100,
                  h = h,
                  X = X,
                  Z = Z,
                  y = y, 
                  model = "logistic"
)

list.save(hyper_par, "study6/logistic_hyper_par.rds")

set.seed(NULL)

# PARNI-LA
alg_par <- list(N = 10000,  # number of iterations
                Nb = 2000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                n_chain = 1,
                # omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                omega_adap = "rm",omega_par = c(-0.7, 0.6), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                method = "CPM",
                n_particles = 5, paricle_corr = 0.99,
                use_ALA2 = FALSE,
                use_ALA = FALSE,
                iter_approx = 1
)




results <- GLM_PARNI(hyper_par = hyper_par, alg_par = alg_par)
list.save(results, "study6/logistic_LA.rds")


barplot(results$estm_PIPs, ylim = c(0,1))
plot(results$log_post_trace[,1], type = "l")
results$CPU_time
results$ESJD


# barplot(results_CPM$estm_PIPs, ylim = c(0,1))
# barplot(results_CPM$estm_PIPs[1:50], ylim = c(0,1))
# barplot(results_CPM$approx_PIPs, ylim = c(0,1))
# barplot(results_CPM$ad_PIPs, ylim = c(0,1))
# plot(results_CPM$log_post_trace[,1], type = "l")
# max(results_CPM$log_post_trace)
# plot(results_CPM$model_size_trace[,1], type = "l")
# plot(results_CPM$omegas, type = "l")
# results_CPM$CPU_time
# results_CPM$acc_rate
# results_CPM$ESJD
# results_CPM$infs







# PARNI-adaptiveALA
alg_par <- list(N = 10000,  # number of iterations
                Nb = 2000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                n_chain = 1,
                # omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                omega_adap = "rm",omega_par = c(-0.7, 0.6), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                method = "CPM",
                n_particles = 5, paricle_corr = 0.99,
                use_ALA2 = TRUE,
                use_ALA = FALSE,
                iter_approx = 1
)




results <- GLM_PARNI(hyper_par = hyper_par, alg_par = alg_par)
list.save(results, "study6/logistic_aALA.rds")


barplot(results$estm_PIPs, ylim = c(0,1))
plot(results$log_post_trace[,1], type = "l")
results$CPU_time
results$ESJD




# PARNI-ALA
alg_par <- list(N = 10000,  # number of iterations
                Nb = 2000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                n_chain = 1,
                # omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                omega_adap = "rm",omega_par = c(-0.7, 0.6), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                method = "CPM",
                n_particles = 5, paricle_corr = 0.99,
                use_ALA2 = FALSE,
                use_ALA = TRUE,
                iter_approx = 1
)




results <- GLM_PARNI(hyper_par = hyper_par, alg_par = alg_par)
list.save(results, "study6/logistic_ALA.rds")


barplot(results$estm_PIPs, ylim = c(0,1))
plot(results$log_post_trace[,1], type = "l")
results$CPU_time
results$ESJD




# survival dataset
set.seed(3564)



library("flexsurv")


n <- 500
p <- 500
SNR <- 1
rho <- 0.6


X <- matrix(rnorm(n*p), ncol = p, nrow = n)

for (j in 2:p) {
  X[,j] <- rho*X[,j-1] + sqrt(1-rho^2)*X[,j]
}

beta <- SNR*c(2,-3, 2, 2,-3, 3,-2, 3,-2, 3, rep(0, p-10))
alpha <- 0

beta[1:10]

X <- scale(X)
mu <- alpha + X %*% beta
Z <- as.matrix(rep(1, n), nrow = n)

t <- rgengamma(n, -mu, sigma = 0.8, Q = -2)

boxplot(t)


sum(t > 500)/n


d <- t <= 500
t <- ifelse(t <= 500, t, 500)


weibull_hyper_par <- list(n = n,
                  p = p,
                  p_z = 1,
                  g = 1,
                  var_intercept = 100,
                  log_k_var = 10^4,
                  h = 10/p,
                  X = X,
                  Z = Z,
                  t = t,
                  d = d,
                  model = "Weibull"
)

list.save(weibull_hyper_par, "study6/weibull_hyper_par.rds")




cox_hyper_par <- list(n = n,
                  p = p,
                  p_z = 0,
                  g = 10,
                  h = 10/p,
                  X = X,
                  t = t,
                  d = d,
                  model = "Cox"
)


list.save(cox_hyper_par, "study6/weibull_hyper_par.rds")



set.seed(NULL)



# weibull
# PARNI-LA
alg_par <- list(N = 10000,  # number of iterations
                Nb = 2000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                n_chain = 1,
                # omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                omega_adap = "rm",omega_par = c(-0.7, 0.6), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                method = "CPM",
                n_particles = 5, paricle_corr = 0.99,
                use_ALA2 = FALSE,
                use_ALA = FALSE,
                iter_approx = 1
)




results <- GLM_PARNI(hyper_par = weibull_hyper_par, alg_par = alg_par)
list.save(results, "study6/weibull_LA.rds")

barplot(results$estm_PIPs, ylim = c(0,1))
plot(results$log_post_trace[,1], type = "l")
results$CPU_time
results$ESJD






# PARNI-adaptiveALA
alg_par <- list(N = 10000,  # number of iterations
                Nb = 2000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                n_chain = 1,
                # omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                omega_adap = "rm",omega_par = c(-0.7, 0.6), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                method = "CPM",
                n_particles = 5, paricle_corr = 0.99,
                use_ALA2 = TRUE,
                use_ALA = FALSE,
                iter_approx = 1
)




results <- GLM_PARNI(hyper_par = weibull_hyper_par, alg_par = alg_par)
list.save(results, "study6/weibull_aALA.rds")


barplot(results$estm_PIPs, ylim = c(0,1))
plot(results$log_post_trace[,1], type = "l")
results$CPU_time
results$ESJD


# PARNI-ALA
alg_par <- list(N = 10000,  # number of iterations
                Nb = 2000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                n_chain = 1,
                # omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                omega_adap = "rm",omega_par = c(-0.7, 0.6), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                method = "CPM",
                n_particles = 5, paricle_corr = 0.99,
                use_ALA2 = FALSE,
                use_ALA = TRUE,
                iter_approx = 1
)




results <- GLM_PARNI(hyper_par = weibull_hyper_par, alg_par = alg_par)
list.save(results, "study6/weibull_ALA.rds")


barplot(results$estm_PIPs, ylim = c(0,1))
plot(results$log_post_trace[,1], type = "l")
results$CPU_time
results$ESJD











# cox
# PARNI-LA
alg_par <- list(N = 10000,  # number of iterations
                Nb = 2000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                n_chain = 1,
                # omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                omega_adap = "rm",omega_par = c(-0.7, 0.6), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                method = "CPM",
                n_particles = 5, paricle_corr = 0.99,
                use_ALA2 = FALSE,
                use_ALA = FALSE,
                iter_approx = 1
)




results <- GLM_PARNI(hyper_par = cox_hyper_par, alg_par = alg_par)
list.save(results, "study6/cox_LA.rds")


barplot(results$estm_PIPs, ylim = c(0,1))
plot(results$log_post_trace[,1], type = "l")
results$CPU_time
results$ESJD







# PARNI-adaptiveALA
alg_par <- list(N = 10000,  # number of iterations
                Nb = 2000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                n_chain = 1,
                # omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                omega_adap = "rm",omega_par = c(-0.7, 0.6), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                method = "CPM",
                n_particles = 5, paricle_corr = 0.99,
                use_ALA2 = TRUE,
                use_ALA = FALSE,
                iter_approx = 1
)




results <- GLM_PARNI(hyper_par = cox_hyper_par, alg_par = alg_par)
list.save(results, "study6/cox_aALA.rds")


barplot(results$estm_PIPs, ylim = c(0,1))
plot(results$log_post_trace[,1], type = "l")
results$CPU_time
results$ESJD




# PARNI-ALA
alg_par <- list(N = 10000,  # number of iterations
                Nb = 2000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                n_chain = 1,
                # omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                omega_adap = "rm",omega_par = c(-0.7, 0.6), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                method = "CPM",
                n_particles = 5, paricle_corr = 0.99,
                use_ALA2 = FALSE,
                use_ALA = TRUE,
                iter_approx = 1
)




results <- GLM_PARNI(hyper_par = cox_hyper_par, alg_par = alg_par)
list.save(results, "study6/cox_ALA.rds")



barplot(results$estm_PIPs, ylim = c(0,1))
plot(results$log_post_trace[,1], type = "l")
results$CPU_time
results$ESJD










