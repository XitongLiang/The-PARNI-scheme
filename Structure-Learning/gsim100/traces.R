# traceplot


X <- gsim100
colnames(X) <- NULL
X <- scale(X)
X <- t(X)


library(pcalg)
suffStat <- list(C = cor(t(X)), n = nrow(t(X)))
varNames <- as.character(1:100)
skel.dataset <- skeleton(suffStat, indepTest = gaussCItest, labels = varNames,
                         alpha = 0.05)
skel.W <- as(skel.dataset@graph, "matrix")

DAG_heatmap(skel.W)



p <- nrow(X)

hyper_par <- list(n = ncol(X),
                  p = nrow(X),
                  g = 10,
                  h = 1/100,
                  X = X,
                  XtX = X %*% t(X))


library(scales)
library(ggthemes)

# skeleton from PC 0.05


# ADR

alg_par <- list(N = 80000, 
                Nb = 16000,
                n_chain = 1,
                verbose = TRUE,
                store_chains = FALSE,
                H = skel.W)


results_ADR <- ADR(alg_par, hyper_par)
plot(results_ADR$log_post_trace, type = "l")
max(results_ADR$log_post_trace)

results_ADR$CPU_time

trace_ADR <- results_ADR$log_post_trace[c(1, seq(4, 80000, by = 4)),1]    


# PARNI

alg_par <- list(N = 20000, 
                Nb = 4000,
                kappa = 0.01,
                n_chain = 1,
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
                H = skel.W)


results_PARNI <- PARNI(alg_par, hyper_par)
plot(results_PARNI$log_post_trace, type = "l")
max(results_PARNI$log_post_trace)

results_PARNI$CPU_time

# partitionMCMC
trace_PARNI <- results_PARNI$log_post_trace[,1]


time.start <- Sys.time()
results_partition <- partitionMCMC(scorepar = score_par,
                                   startspace = skel.W,
                                   scoreout = TRUE,
                                   stepsave = 1,
                                   iterations = 20000,
                                   verbose = TRUE
                                   )
time_partition <- Sys.time() - time.start



plot(results_partition$trace, type = "l")

trace_partition <- results_partition$trace




plot(trace_PARNI, type = "l", col = "red")
lines(trace_ADR, type = "l", col = "blue")
lines(trace_partition, type = "l", col = "green")


traces <- data.frame(trace_PC_PARNI = trace_PARNI,
                    trace_PC_ADR = trace_ADR, 
                    trace_PC_partition = trace_partition)


# skeleton from iterativeMCMC


# ADR

alg_par <- list(N = 80000, 
                Nb = 16000,
                n_chain = 1,
                verbose = TRUE,
                store_chains = FALSE,
                H = postskel_iterativeMCMC$skel)


results_ADR <- ADR(alg_par, hyper_par)
plot(results_ADR$log_post_trace, type = "l")
max(results_ADR$log_post_trace)

results_ADR$CPU_time

trace_ADR <- results_ADR$log_post_trace[c(1, seq(4, 80000, by = 4)),1]    


# PARNI

alg_par <- list(N = 20000, 
                Nb = 4000,
                kappa = 0.01,
                n_chain = 1,
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
plot(results_PARNI$log_post_trace, type = "l")
max(results_PARNI$log_post_trace)

results_PARNI$CPU_time

# partitionMCMC
trace_PARNI <- results_PARNI$log_post_trace[,1]


time.start <- Sys.time()
results_partition <- partitionMCMC(scorepar = score_par,
                                   startspace = postskel_iterativeMCMC$skel,
                                   scoreout = TRUE,
                                   stepsave = 1,
                                   iterations = 20000,
                                   verbose = TRUE
)
time_partition <- Sys.time() - time.start



plot(results_partition$trace, type = "l")

trace_partition <- results_partition$trace




plot(trace_PARNI, type = "l", col = "red")
lines(trace_ADR, type = "l", col = "blue")
lines(trace_partition, type = "l", col = "green")


traces <- traces %>%
  mutate(trace_iterative_PARNI = trace_PARNI,
                    trace_iterative_ADR = trace_ADR, 
                    trace_iterative_partition = trace_partition)


# write_csv(traces, "gsim100/traces.csv")

traces %>%
  mutate(iterations = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("trace_"),
               names_to = c("Skeleton", "Sampler"),
               names_pattern = "trace_(.*)_(.*)",
               values_to = "log_pmp") %>%
  mutate(Skeleton = factor(Skeleton, levels = c("PC", "iterative")),
         Sampler = ifelse(Sampler == "PARNI", "PARNI-DAG",
                          ifelse (Sampler == "partition", "partitionMCMC", Sampler))) %>%
  ggplot(aes(x = iterations, y = log_pmp, col = Sampler)) +
  geom_line() +
  facet_wrap(~ Skeleton) +
  ylab("log posterior model probability") + 
  theme_igray() +
  scale_colour_colorblind() # + theme(legend.position="bottom")



