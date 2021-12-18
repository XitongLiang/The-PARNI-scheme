library(tidyverse)

# load dataset information
settings <- read_csv("simulation/simulation.csv", col_types = cols())

# choose a dataset between 1 and 24
i <- 1

# define file name of the dataset
file.name <- settings$file.name[i]

# import dataset
dataset <- read_csv(file.name, col_
                    types = cols(), progress = FALSE)
y <- pull(dataset %>% select(V1))
X <- as.matrix(dataset %>% select(-V1))
colnames(X) <- NULL

# setting hyperparameters
g <- 9
h <- 10/settings$p[i]
n <- settings$n[i]
p <- settings$p[i]
SNR <- settings$SNR[i]
remove(dataset)

# define model hyperparameters
hyper_par <- list(n = n,
                  p = p,
                  g = g,
                  h = h,
                  X = X,
                  y = y,
                  yty = sum(y^2), 
                    ytX = t(y) %*% X,
                  g_prior_type = "ind")

# remove used variables
remove(dataset, y, X)

