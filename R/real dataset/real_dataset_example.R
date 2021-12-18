library(tidyverse)
library(rlist)

settings <- read_csv("real dataset/real_datasets.csv", col_types = cols())

i <- 1 # choose which dataset you would like to study

# load dataset
hyper_par <- list.load(settings$file[i])


# PARNI
alg_par <- list(N = 3000,  # number of iterations
                Nb = 1000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                target_swap = 0.25, # Parallel tempering target swapping acceptance rate
                n_chain = 25,
                n_temp = 1,
                omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                # omega_adap = "rm",omega_par = c(-0.7, 0.65), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                # bal_fun = sqrt,
                use_rb = TRUE  # whether to use Raoo-Blackwellised PIPs
)

results <- PARNI(hyper_par = hyper_par, alg_par = alg_par)


# Add-delete-swap
alg_par <- list(N = 180000,
                Nb = 60000,
                verbose = TRUE,
                n_chain = 1,
                store_chains = FALSE,
                phi = -0.7,
                # init_taus = rep(-1,3),
                target_swap = 0.237,
                n_temp = 1
)

results <- ADS(alg_par, hyper_par)






# ASI
alg_par <- list(N = 6000,
                Nb = 2000,
                full_adap = FALSE,
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                phi = c(-0.7, -0.7),
                target_prop = 0.2,
                target_swap = 0.25,
                n_chain = 25,
                init_zeta = 0.5,
                n_temp = 1,
                store_chains = FALSE)

results <- ASI(alg_par, hyper_par)
barplot(results$estm_PIPs, ylim = c(0,1))
