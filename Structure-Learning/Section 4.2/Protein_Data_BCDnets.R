#################################
# 1 - ecoli70 Simulated Example #
#################################

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
source("utils/update_LA_DAG.R")
source("utils/sample_ind_DAG.R")
source("utils/logit_e.R")
source("utils/is.DAG_adjmat.R")
source("utils/DAG_heatmap.R")
source("utils/marPIPs_DAG_H.R")
source("utils/log_llh_DAG.R")
source("utils/update_LA_DAG_sawp.R")
source("utils/other_functions.R")


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
library(igraph)

# Brier Score & accuracy
brier_score = function(true_W, est_W) { return(mean((true_W - est_W)^2)) }
accuracy = function(true_W, est_W) { sum(true_W == ifelse(est_W >= thresh_prob, 
                                                          1, 0)) / (nrow(true_W)*ncol(true_W)) }
MC_se <- function(x, B) qt(0.975, B-1)*sd(x, na.rm = T)/sqrt(B)



# Load data generating process functions
source("utils/datageneration.R")

#### Load data (Protein data - https://pubmed.ncbi.nlm.nih.gov/15845847/)
W = as.matrix(read.csv(paste0(dirname(getwd()), '/Data/sachs_W.csv')))[, -1]
X = as.matrix(read.csv(paste0(dirname(getwd()), '/Data/sachs_X.csv')))[, -1]
X = scale(X)

# W_graph = graph_from_adjacency_matrix(W, "directed")
# 
# plot.igraph(W_graph, layout = layout.reingold.tilford(W_graph),
#             vertex.color = rgb(0.9,0.3,0.3,0.3), vertex.frame.color="red",
#             vertex.size=5, vertex.shape="circle", vertex.label.color="black",
#             vertex.label=parse(text= paste("italic(X[", 1:nrow(W), "])", sep="")), 
#             edge.color="black", edge.width=1,
#             edge.arrow.size=0.3, edge.curved=0,
#             vertex.label.cex = 0.6,
#             asp = 0.4,
#             margin = -0.1)



#### Storing
models_compared = c("PC", "GES", "LINGAM", "Order", 
                    "Iterative", "Restr_ADR", "Restr_PARNI")
metrics = c("Brier", "TPR", "FPR", "SHD") 

metrics_list = as.list(rep(c(NA), length(metrics)))
metrics_list  = setNames(metrics_list, metrics)

final_list = list()
for (mod in models_compared) { final_list[[mod]] = metrics_list }

rm(metrics_list)

time_mark = Sys.time()


#### Run Simulations
for (b in 1:B) {
  
  set.seed(b + b*10)
  
  library(rlist)
  dataset <- list.load("utils/dataset.rds")
  default_hyper_par <- dataset$hyper_par
  
  
  hyper_par <- list(n = as.integer(nrow(X)), 
                    p = ncol(W),
                    g = ceiling(ncol(W)/2),
                    h = default_hyper_par$h,
                    X = X,
                    log_llh_type = default_hyper_par$log_llh_type)
  
  
  
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
                            iterations = 15000, chainout = T)
  time_order = Sys.time() - time_start
  
  
  pip_order = as.matrix(edgep(results_order, pdag = FALSE, burnin = 0.34))
  
  est_dirW_order = ifelse(pip_order >= 0.5, 1, 0)
  
  
  # Store directed W results
  final_list[["Order"]][["Brier"]][b] = brier_score(W, pip_order)
  
  for (metr in c("TPR", "FPR", "SHD")) {
    final_list[["Order"]][[metr]][b] = compareDAGs(est_dirW_order, W)[metr]
  }
  
  
  
  
  ####
  #### PARTITION MCMC (ALREADY PAIRED WITH ITERATIVE PROCEDURE) REALLY STRUGGLE TO RUN ON 50 NODES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ####
  
  
  # #### 5) (Restricted) Partition MCMC 
  # 
  # # PartitionMCMC starts by  restricting the space with the PC algorithm 
  # # skeleton - the non restricted version is feasible only for low number of nodes)
  # partition_score = scoreparameters(scoretype = "bge",
  #                                   data = X)
  
  
  
  
  #### 6) Iterative(/order) MCMC 
  
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
                                      MAP = FALSE, 
                                      iterations = 15000, chainout = T, plus1it = 15)
    time_iterative = Sys.time() - time_start
    
    
    pip_Iterative = as.matrix(edgep(results_iterative, pdag = FALSE, burnin = 0.34))
    est_dirW_Iterative = ifelse(pip_Iterative >= 0.5, 1, 0)
    
    # Store directed W results
    final_list[["Iterative"]][["Brier"]][b] = brier_score(W, pip_Iterative)
    
    for (metr in c("TPR", "FPR", "SHD")) {
      final_list[["Iterative"]][[metr]][b] = compareDAGs(est_dirW_Iterative, W)[metr]
    }
    
  })
  
  
  
  #### 7) Restricted ADR
  # Skeleton PC
  p <- ncol(W)
  suffStat <- list(C = cor(X), n = nrow(X))
  varNames <- as.character(1:p)
  skel.dataset <- skeleton(suffStat, indepTest = gaussCItest, labels = varNames, 
                           alpha = 0.01)
  skel.W <- as(skel.dataset@graph, "matrix")
  
  
  alg_par <- list(N = 15000, 
                  Nb = 5000,
                  n_chain = 20,
                  verbose = TRUE,
                  store_chains = FALSE,
                  H = skel.W)
  
  
  time_start = Sys.time()
  results_ADR <- ADR(alg_par, hyper_par)
  time_ADR = Sys.time() - time_start
  
  
  est_dirW_ADR = ifelse(results_ADR$estm_PIPs >= 0.25, 1, 0)
  
  # Store directed W results
  final_list[["Restr_ADR"]][["Brier"]][b] = brier_score(W, results_ADR$estm_PIPs)
  
  for (metr in c("TPR", "FPR", "SHD")) {
    final_list[["Restr_ADR"]][[metr]][b] = compareDAGs(est_dirW_ADR, W)[metr]
  }
  
  
  
  
  #### 8) Restricted PARNI
  alg_par <- list(N = 4000, 
                  Nb = 2000,
                  kappa = 0.001,
                  n_chain = 20,
                  verbose = TRUE,
                  omega_adap = "rm", omega_par = c(-0.7, 0.5),
                  eps = 1/(p*(p-1)),
                  store_chains = FALSE,
                  bal_fun = function(x) {pmin(1,x)},
                  H = skel.W)
  
  
  time_start = Sys.time()
  results_PARNI <- PARNI(alg_par, hyper_par)
  time_PARNI = Sys.time() - time_start
  
  
  est_dirW_PARNI = ifelse(results_PARNI$estm_PIPs >= 0.25, 1, 0)
  
  
  # Store directed W results
  final_list[["Restr_PARNI"]][["Brier"]][b] = brier_score(W, results_PARNI$estm_PIPs)
  
  for (metr in c("TPR", "FPR", "SHD")) {
    final_list[["Restr_PARNI"]][[metr]][b] = compareDAGs(est_dirW_PARNI, W)[metr]
  }
  
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
  filter(metric == "Brier" & N == 50)


# write.csv(df_results,
#           paste0("./Results/SimExp_Run6_Sparse_B", B,".csv"))


levels(df_results$model)[match("Restr_ADR", levels(df_results$model))] = "ADR"
levels(df_results$model)[match("Restr_PARNI", levels(df_results$model))] = "PARNI-DAG"


ggplot(df_results, aes(x = N)) + geom_line(aes(y = ave, color = model), size = 1) +
  # geom_errorbar(aes(ymin = (ave - SE), ymax = (ave + SE), color = model)) +
  # geom_ribbon(aes(ymin = (ave - SE), ymax = (ave + SE), color = model), fill = "grey", alpha = 0.2) +
  ggh4x::facet_grid2(c("edge_marks", "metric"), 
                     scales = "free", independent = "all") +
  theme_minimal() + ylab("") + labs(color = "Model:") +
  scale_color_manual(values=hcl.colors(n=length(models_compared), palette = "Dynamic", rev = T))


