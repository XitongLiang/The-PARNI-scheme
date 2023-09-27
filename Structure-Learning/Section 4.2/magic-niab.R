############################################
# magic-niab Simulated Example 
# Ref: https://pubmed.ncbi.nlm.nih.gov/16646851/
############################################

rm(list = ls())

# retrieving script path 
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

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



##### Install BiDAG (Uncomment if needed)
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install("graph")
# BiocManager::install("Rgraphviz")
# BiocManager::install("RBGL")
# 
# install.packages("BiDAG")

# Load packages
library(BiDAG)
library(igraph)
library(rlist)
library(tidyverse)


# Brier Score & accuracy
brier_score = function(true_W, est_W) { return(mean((true_W - est_W)^2)) }
accuracy = function(true_W, est_W) { sum(true_W == ifelse(est_W >= thresh_prob, 
                                                          1, 0)) / (nrow(true_W)*ncol(true_W)) }
MC_se <- function(x, B) qt(0.975, B-1)*sd(x, na.rm = T)/sqrt(B)



# Load data generating process functions
source("utils/datageneration.R")

#### Load data (magic-niab - https://www.degruyter.com/document/doi/10.2202/1544-6115.1175/html)
magic_niab = readRDS(paste0(dirname(getwd()), '/Data/magic-niab.rds'))
magic_niab_W = bnlearn::as.igraph(magic_niab)

W = as.matrix(as_adjacency_matrix(magic_niab_W))
num_edges = sum(W)


#### Options
N = 100 # sample size
B = 20 # Replications


#### Storing
models_compared = c("PC", "GES", "LINGAM", "Order", "Iterative", 
                    "Partition", "ADR", "PARNI")
metrics = c("Brier", "TPR", "FPR", "SHD") 

metrics_list = as.list(rep(c(NA), length(metrics)))
metrics_list  = setNames(metrics_list, metrics)

final_list = list()
for (mod in models_compared) { final_list[[mod]] = metrics_list }

rm(metrics_list)


time_mark = Sys.time()


#### Run Simulations
for (b in 1:B) {
  
  
  cat("\n\n\n\n###### Iteration ", b, "\n\n\n\n")
  set.seed(b + b*20)
  
  
  
  #   #   1.1   Data generation from graph
  X = as.matrix( bnlearn::rbn(magic_niab, N) )
  X = scale(X)
  
  # X = t(BCD_generation_laplace(W = W, n = N, const = 0.5, mu = 0, beta_l = 1))
  # X = scale(X)
  
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
    final_list[["PC"]][[metr]][b] = compareDAGs(est_dirW_PC, W)[metr]
  }
  
  
  
  ####
  #### WARNING  Not able to converge with p = n, so discard the first prolly
  ####
  
  #### 2) GES algorithm
  GES_score <- new("GaussL0penObsScore", data = X,
                   lambda = 0.5*log(nrow(X)))
  
  GES_ = ges(GES_score)
  
  
  str(GES_, max=2)
  est_dirW_GES = as.matrix( as(as(GES_$essgraph,"graphNEL"),"Matrix") )
  
  
  # Store directed W results
  for (metr in c("TPR", "FPR", "SHD")) {
    final_list[["GES"]][[metr]][b] = compareDAGs(est_dirW_GES, W)[metr]
  }
  
  
  
  
  ####
  # # LINGAM AS WELL MIGHT NOT WORK AS SINGULARITY AND TOO HIGH-DIMENSIONAL PROBLEM
  ####
  
  
  #### 3) LINGAM algorithm
  try({
    
    lingam.fit = lingam(X)
    est_dirW_LINGAM = ifelse(lingam.fit$Bpruned == 0, 0, 1)
    
    # Store directed W results
    for (metr in c("TPR", "FPR", "SHD")) {
      final_list[["LINGAM"]][[metr]][b] = compareDAGs(est_dirW_LINGAM, W)[metr]
    }
    
  })
  
  
  
  
  
  #### 4) (Restricted) Order MCMC 
  
  # OrderMCMC starts by restricting the search space with the PC algorithm 
  # skeleton - the non restricted version is feasible only for low number of nodes)
  order_score = scoreparameters(scoretype = "bge",
                                data = X)
  
  time_start = Sys.time()
  results_order = orderMCMC(scorepar = order_score, 
                            MAP = FALSE, 
                            plus1 = FALSE,
                            iterations = 200000, chainout = T)
  time_order = Sys.time() - time_start
  
  
  pip_order = as.matrix(edgep(results_order, pdag = FALSE, burnin = 0.4))
  est_dirW_order = as.matrix(results_order$DAG)
  
  
  # Store directed W results
  final_list[["Order"]][["Brier"]][b] = brier_score(W, pip_order)
  
  for (metr in c("TPR", "FPR", "SHD")) {
    final_list[["Order"]][[metr]][b] = compareDAGs(est_dirW_order, W)[metr]
  }
  
  
  
  #### 5) Iterative(/order) MCMC 
  
  # IterativeMCMC implements the orderMCMC scheme on the restricted starting
  # space defined by the skeleton of the PC algorithm - BUT implements it for
  # more iterations (approx 20), where at each iteration the initial search space 
  # is expanded by allowing each node to have up to one additional parent each,
  # not included in the previous search space
  try({
    
    iterative_score = scoreparameters(scoretype = "bge",
                                      data = X)
    
    time_start = Sys.time()
    results_iterative = iterativeMCMC(scorepar = iterative_score, 
                                      MAP = FALSE, stepsave = 100,
                                      iterations = 200000, chainout = T, plus1it = 10)
    time_iterative = Sys.time() - time_start
    
    
    pip_Iterative = as.matrix(edgep(results_iterative, pdag = FALSE, burnin = 0.4))
    est_dirW_Iterative = as.matrix(results_iterative$DAG)
    
    # Store directed W results
    final_list[["Iterative"]][["Brier"]][b] = brier_score(W, pip_Iterative)
    
    for (metr in c("TPR", "FPR", "SHD")) {
      final_list[["Iterative"]][[metr]][b] = compareDAGs(est_dirW_Iterative, W)[metr]
    }
    
  })
  
  
  
  
  # #### 6) (Restricted) Partition MCMC 
  # 
  # # PartitionMCMC starts by  restricting the space with the PC algorithm 
  # # skeleton - the non restricted version is feasible only for low number of nodes)
  try({
    
    partition_score = scoreparameters(scoretype = "bge",
                                      data = X)
    
    time_start = Sys.time()
    results_partition = partitionMCMC(scorepar = partition_score,
                                      verbose = T,
                                      iterations = 200000, 
                                      startspace = as.matrix(results_iterative$DAG) + t(as.matrix(results_iterative$DAG)))
    time_partition = Sys.time() - time_start
    
    pip_partition = as.matrix(edgep(results_partition, pdag = FALSE, burnin = 0.4))
    est_dirW_partition = as.matrix(results_partition$DAG)
    
    
    # Store directed W results
    final_list[["Partition"]][["Brier"]][b] = brier_score(W, pip_partition)
    
    for (metr in c("TPR", "FPR", "SHD")) {
      final_list[["Partition"]][[metr]][b] = compareDAGs(est_dirW_partition, W)[metr]
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
  final_list[["ADR"]][["Brier"]][b] = brier_score(W, results_ADR$estm_PIPs)
  
  for (metr in c("TPR", "FPR", "SHD")) {
    final_list[["ADR"]][[metr]][b] = compareDAGs(est_dirW_ADR, W)[metr]
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
  final_list[["PARNI"]][["Brier"]][b] = brier_score(W, results_PARNI$estm_PIPs)
  
  for (metr in c("TPR", "FPR", "SHD")) {
    final_list[["PARNI"]][[metr]][b] = compareDAGs(est_dirW_PARNI, W)[metr]
  }
  
  
} # B replication


time_expired = Sys.time() - time_mark

print(time_expired)


df_results = expand.grid(model = models_compared,
                         metric = metrics,
                         ave = NA,
                         SE = NA)


for (mod in models_compared) {
  for (metric in metrics) {
    
    cond = df_results$model == mod & df_results$metric == metric
    
    df_results[cond, "ave"] = mean(final_list[[mod]][[metric]], na.rm = T)
    df_results[cond, "SE"] = MC_se(final_list[[mod]][[metric]], B)
  }
}


df_results %>%
  filter(metric == "SHD")


write.csv(df_results,
          paste0("../Results/magic-niab_B", B,".csv"))

