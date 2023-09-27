
postskel_iterativeMCMC <- list.load("gsim100/gsim100_postskel_iterativeMCMC.rds")


alg_par <- list(N = 30000*12, 
                Nb = 6000*12,
                kappa = 0.01,
                n_chain = 20,
                verbose = TRUE,
                PIPs_update = TRUE,
                # second value of omega_par is the expected neighoburhood size
                # which will affect computational costs and the effciencies
                # omega_adap = "f", omega_init = 0.1,
                omega_adap = "rm", omega_par = c(-0.7, 10), use_logit_e = TRUE,
                # omega_adap = "kw", omega_init = 0.2, omega_par = c(-1, -0.5), use_logit_e = FALSE,
                eps = 1/(hyper_par$p*(hyper_par$p-1)),
                store_chains = FALSE,
                bal_fun = function(x) {pmin(1,x)},
                H = postskel_iterativeMCMC$skel)


results_PARNI <- PARNI(alg_par, hyper_par)

results_PARNI$CPU_time

list.save(results_PARNI, "gsim100/PARNI_long.rds")

