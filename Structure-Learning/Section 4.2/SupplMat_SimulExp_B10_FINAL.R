###################################
# 1 - Static BN Simulated Example #
###################################

# Notes:
# "Directed" version of the models is learning the DAG directly, while metrics on
# the "Undirected" learn the respective skeleton (different than equivalent class and CPDAG)
# Definitions: https://kevinbinz.com/tag/cpdag/ or/and 
# https://www.hiit.fi/wp-content/uploads/2018/04/Learning-Markov-Equivalence-S29.pdf
#
#
# While Accuracy, Brier and SHD are affected by DAG or CPDAG, TPR and FPR are not!


rm(list = ls())

# retrieving script path 
library("rstudioapi") 
setwd(dirname(getSourceEditorContext()$path))


#   #   0   Packages & function loading -------------------------------------------

library(rlist)
library(tidyverse)
library(pcalg)

# ADR + PARNI functions

# Restricted
source("utils/PARNI.R")
source("utils/ADR.R")

source("utils/Compute_LA_DAG.R")
# source("utils/update_LA_DAG.R")
source("utils/sample_ind_DAG.R")
source("utils/logit_e.R")
source("utils/is.DAG_adjmat.R")
source("utils/DAG_heatmap.R")
source("utils/marPIPs_DAG_H.R")
source("utils/log_llh_DAG.R")
source("utils/log_llh_DAG_update_table.R")
source("utils/update_LA_DAG_sawp.R")
source("utils/other_functions.R")

source("utils/marPIPs_bge_H.R")
source("utils/log_llh_bge_table.R")


# # Install BiDAG (RUN ONLY IF NEEDED)
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install("graph")
# BiocManager::install("Rgraphviz")
# BiocManager::install("RBGL")
# 
# install.packages("BiDAG")

# BIDAGs algorithms
library(BiDAG)


# Brier Score & accuracy
brier_score = function(true_W, est_W) { return(mean((true_W - est_W)^2)) }
accuracy = function(true_W, est_W) { sum(true_W == ifelse(est_W >= thresh_prob, 1, 0)) / (nrow(true_W)*ncol(true_W)) }
MC_se <- function(x, B) qt(0.975, B-1)*sd(x, na.rm = T)/sqrt(B)


# Load data generating process functions
source("utils/datageneration.R")


#### Options
N = c(50, 100) # sample size
d = 50 # number of variables
true_pedge = 0.1 # true probability of an edge b/w variables

B = 10 # Replications


#   #   1   Graph Generation -----------------------------------------------
W = DAG_adjmat_generation(d = d, p = true_pedge)
num_edges = sum(W)


#### Storing
models_compared = c("PC", "GES", "LINGAM", "Order", "Iterative", 
                    "Partition", "ADR", "PARNI")
metrics = c("Brier", "TPR", "FPR", "SHD") 


metrics_list = as.list(rep(c(NA), length(metrics)))
metrics_list = setNames(metrics_list, metrics)

results_list = list()
final_list = list()
for (mod in models_compared) { results_list[[mod]] = metrics_list }
for (n in N) { final_list[[paste(n)]] = results_list }

rm(results_list, metrics_list)


time_mark = Sys.time()

# Sample size N
for (n in N) {
  
  
  #### Run Simulations
  for (b in 1:B) {
    
    
    cat("\n\n\n++++ Sample N:", n, "\n\n###### Iteration ", b, "\n\n\n\n")
    set.seed(b + b*10)
    
    
    #   #   1.1   Data generation from graph
    X = t(BCD_generation(W = W, n = n, const = 0.5))
    X = scale(X)
    
    
    p = nrow(W)
    
    
    hyper_par <- list(n = nrow(X),
                      p = ncol(X),
                      g = 10,
                      h = 1/100,
                      X = t(X),
                      XtX = t(X) %*% X)
    
    
    
    #   #   #   2. Run Models  -----------------------------------------------
    
    #### 1) PC algorithm
    suffStat <- list(C = cor(X), n = nrow(X))
    pc.gmG <- pc(suffStat, indepTest = gaussCItest,
                 p = ncol(X), alpha = 0.01)
    
    
    est_dirW_PC = as.matrix(as(pc.gmG, "amat"))
    
    # Store directed W results
    for (metr in c("TPR", "FPR", "SHD")) {
      final_list[[paste(n)]][["PC"]][[metr]][b] = compareDAGs(est_dirW_PC, W)[metr]
    }
    
    
    
    
    
    ####
    #### WARNING  Not able to converge with p = n, so discard the first prolly
    ####
    
    #### 2) GES algorithm
    GES_score <- new("GaussL0penObsScore", data = X,
                     lambda = floor(2*log(nrow(X))))
    
    GES_ = ges(GES_score)
    
    
    str(GES_, max=2)
    est_dirW_GES = as.matrix( as(as(GES_$essgraph,"graphNEL"),"Matrix") )
    
    
    # Store directed W results
    for (metr in c("TPR", "FPR", "SHD")) {
      final_list[[paste(n)]][["GES"]][[metr]][b] = compareDAGs(est_dirW_GES, W)[metr]
    }   
    
    
    
    
    
    ####
    # # LINGAM AS WELL MIGHT NOT WORK AS SINGULARITY AND TOO HIGH-DIMENSIONAL PROBLEM
    ####
    
    
    #### 3) LINGAM algorithm
    try({
      
      lingam.fit = lingam(X)
      
      est_dirW_LINGAM = ifelse(lingam.fit$Bpruned == 0, 0, 1)
      est_undirW_LINGAM = as.matrix(igraph::as_adjacency_matrix(igraph::graph_from_adjacency_matrix(est_dirW_LINGAM, "undirected")))
      
      # Store directed W results
      for (metr in c("TPR", "FPR", "SHD")) {
        final_list[[paste(n)]][["LINGAM"]][[metr]][b] = compareDAGs(est_dirW_LINGAM, W)[metr]
      }
      
    })
    
    
    
    
    
    #### 4) (Restricted) Order MCMC 
    
    # OrderMCMC starts by restricting the search space with the PC algorithm 
    # skeleton - the non restricted version is feasible only for low number of nodes)
    order_score = scoreparameters(scoretype = "bge",
                                  data = X, 
                                  bgepar = list(edgepf = floor(2*log(p))))
    
    time_start = Sys.time()
    results_order = orderMCMC(scorepar = order_score, 
                              MAP = FALSE, 
                              plus1 = FALSE,
                              iterations = 200000, chainout = T)
    time_order = Sys.time() - time_start
    
    
    pip_order = as.matrix(edgep(results_order, pdag = FALSE, burnin = 0.4))
    est_dirW_order = as.matrix(results_order$DAG)
    
    
    # Store directed W results
    final_list[[paste(n)]][["Order"]][["Brier"]][b] = brier_score(W, pip_order)
    
    for (metr in c("TPR", "FPR", "SHD")) {
      final_list[[paste(n)]][["Order"]][[metr]][b] = compareDAGs(est_dirW_order, W)[metr]
    }
    
    
    
    #### 5) Iterative(/order) MCMC 
    
    # IterativeMCMC implements the orderMCMC scheme on the restricted starting
    # space defined by the skeleton of the PC algorithm - BUT implements it for
    # more iterations (approx 20), where at each iteration the initial search space 
    # is expanded by allowing each node to have up to one additional parent each,
    # not included in the previous search space
    try({
      
      iterative_score = scoreparameters(scoretype = "bge",
                                        data = X,
                                        bgepar = list(edgepf = floor(2*log(p))))
      
      time_start = Sys.time()
      results_iterative = iterativeMCMC(scorepar = iterative_score, 
                                        MAP = FALSE, stepsave = 100,
                                        iterations = 200000, chainout = T, plus1it = 10)
      time_iterative = Sys.time() - time_start
      
      
      pip_Iterative = as.matrix(edgep(results_iterative, pdag = FALSE, burnin = 0.4))
      est_dirW_Iterative = as.matrix(results_iterative$DAG)
      
      # Store directed W results
      final_list[[paste(n)]][["Iterative"]][["Brier"]][b] = brier_score(W, pip_Iterative)
      
      for (metr in c("TPR", "FPR", "SHD")) {
        final_list[[paste(n)]][["Iterative"]][[metr]][b] = compareDAGs(est_dirW_Iterative, W)[metr]
      }
      
    })
    
    
    
    # #### 6) (Restricted) Partition MCMC 
    # 
    # # PartitionMCMC starts by  restricting the space with the PC algorithm 
    # # skeleton - the non restricted version is feasible only for low number of nodes)
    try({
      
      partition_score = scoreparameters(scoretype = "bge",
                                        data = X, 
                                        bgepar = list(edgepf = floor(2*log(p))))
      
      time_start = Sys.time()
      results_partition = partitionMCMC(scorepar = partition_score,
                                        verbose = T,
                                        iterations = 200000, 
                                        startspace = as.matrix(results_iterative$DAG) + t(as.matrix(results_iterative$DAG)))
      time_partition = Sys.time() - time_start
      
      pip_partition = as.matrix(edgep(results_partition, pdag = FALSE, burnin = 0.4))
      est_dirW_partition = as.matrix(results_partition$DAG)
      
      
      # Store directed W results
      final_list[[paste(n)]][["Partition"]][["Brier"]][b] = brier_score(W, pip_partition)
      
      for (metr in c("TPR", "FPR", "SHD")) {
        final_list[[paste(n)]][["Partition"]][[metr]][b] = compareDAGs(est_dirW_partition, W)[metr]
      }  
      
    })
    
    
    
    
    #### 7) Restricted ADR
    skel.W = try(as.matrix(results_iterative$DAG) + t(as.matrix(results_iterative$DAG)))
    
    if ('try-error' %in% class(skel.W)) {
      
      # Skeleton PC
      p <- ncol(W)
      suffStat <- list(C = cor(X), n = nrow(X))
      varNames <- as.character(1:p)
      skel.dataset <- skeleton(suffStat, indepTest = gaussCItest, labels = varNames,
                               alpha = 0.05)
      skel.W <- as(skel.dataset@graph, "matrix")
      
    }
    
    
    #
    alg_par <- list(N = 20000, 
                    Nb = 10000,
                    n_chain = 10,
                    verbose = TRUE,
                    store_chains = FALSE,
                    H = skel.W)
    
    
    time_start = Sys.time()
    results_ADR <- ADR(alg_par, hyper_par)
    time_ADR = Sys.time() - time_start
    
    
    est_dirW_ADR = ifelse(results_ADR$estm_PIPs >= 0.5, 1, 0)
    
    # Store directed W results
    final_list[[paste(n)]][["ADR"]][["Brier"]][b] = brier_score(W, results_ADR$estm_PIPs)
    
    for (metr in c("TPR", "FPR", "SHD")) {
      final_list[[paste(n)]][["ADR"]][[metr]][b] = compareDAGs(est_dirW_ADR, W)[metr]
    }
    
    
    
    
    #### 8) Restricted PARNI
    hyper_par$use_bge <- FALSE
    
    alg_par <- list(N = 6000, 
                    Nb = 4000,
                    kappa = 0.01,
                    n_chain = 10,
                    verbose = TRUE,
                    random_gamma_init = FALSE,
                    PIPs_update = TRUE,
                    omega_adap = "rm", omega_par = c(-0.7, 10), use_logit_e = TRUE,
                    eps = 1/(p*(p-1)),
                    store_chains = FALSE,
                    bal_fun = function(x) {pmin(1,x)},
                    H = skel.W)
    
    
    time_start = Sys.time()
    results_PARNI <- PARNI(alg_par, hyper_par)
    time_PARNI = Sys.time() - time_start
    
    est_dirW_PARNI = ifelse(results_PARNI$estm_PIPs >= 0.5, 1, 0)
    
    # Store directed W results
    final_list[[paste(n)]][["PARNI"]][["Brier"]][b] = brier_score(W, results_PARNI$estm_PIPs)
    
    for (metr in c("TPR", "FPR", "SHD")) {
      final_list[[paste(n)]][["PARNI"]][[metr]][b] = compareDAGs(est_dirW_PARNI, W)[metr]
    }
    
    
  } # B replication
  
} # N sample size


time_expired = Sys.time() - time_mark

print(time_expired)


df_results = expand.grid(N = N,
                         model = models_compared,
                         metric = metrics,
                         ave = NA,
                         SE = NA)


for (n in N) {
  for (mod in models_compared) {
    for (metric in metrics) {
      
      cond = df_results$N == n & df_results$model == mod & df_results$metric == metric
      
      df_results[cond, "ave"] = mean(final_list[[paste(n)]][[mod]][[metric]], na.rm = T)
      df_results[cond, "SE"] = MC_se(final_list[[paste(n)]][[mod]][[metric]], B)
      
    }
  }
}


df_results %>%
  filter(metric == "SHD" & N == 50)

sum(is.na(df_results$ave))

write.csv(df_results,
          paste0("../Results/SupplMat_50_100_SimExp_B", B,".csv"))

# 
# ggplot(df_results, aes(x = N)) + geom_line(aes(y = ave, color = model), size = 1) +
#   geom_errorbar(aes(ymin = (ave - SE), ymax = (ave + SE), color = model)) +
#   ggh4x::facet_grid2(c("edge_marks", "metric"), 
#                      scales = "free", independent = "all", ) +
#   theme_minimal() + ylab("") + labs(color = "Model:") +
#   scale_color_manual(values=hcl.colors(n=length(models_compared), palette = "Dynamic", rev = T))


final_list = final_list

save(final_list, file = "../Results/final_list_B10.RData")
load("../Results/final_list_B10.RData")


df_results = expand.grid(B = 1:B,
                                N = N,
                                model = models_compared,
                                SHD = NA)


for (n in N) {
  for (mod in models_compared) {
    for (b in 1:B) {
      
      cond = df_results$B == b & df_results$N == n & df_results$model == mod
      df_results[cond, 'SHD'] = final_list[[paste(n)]][[mod]][['SHD']][[b]]
      
    }
  }
}


df_results = dplyr::bind_rows(df_results, df_results_150) %>%
  select(-B) %>%
  mutate(N = factor(paste0('N = ', N), levels = c('N = 50', 'N = 100', 'N = 150'))) %>%
  rename(Model = model)


p = ggplot(df_results, aes(x = Model, y = SHD, fill = Model)) +
  geom_violin() + 
  facet_wrap(~ N, scales = 'free_y') +
  theme_minimal()


ggsave("../Results/Violin.pdf", plot = p, width = 14.5, height = 4.5, dpi = 800)
write.csv(df_results,
          paste0("../Results/SupplMat_Result_SimExp_B", B,".csv"))


