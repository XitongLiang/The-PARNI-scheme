PCR_X <- read.table("real dataset/dataset/PCR/gene_expression_data.txt")
PCR_y <- read.table("real dataset/dataset/PCR/gene_expression_target.txt")
PCR_sex <- read.table("real dataset/dataset/PCR/Gender_num.txt")

Z <- as.matrix(cbind(1, PCR_sex))
X <- as.matrix(PCR_X[,-1])
ys <- as.matrix(PCR_y[,-1])
colnames(X) <- NULL
colnames(ys) <- NULL

n <- nrow(X)
p <- ncol(X)

ys <- (diag(n) - Z %*% solve(t(Z) %*% Z) %*% t(Z)) %*% ys
X <- (diag(n) - Z %*% solve(t(Z) %*% Z) %*% t(Z)) %*% X


rm(PCR_X, PCR_y,  PCR_sex, Z)

n <- nrow(X)
p <- ncol(X)
X <- scale(X)
y <- scale(ys[,1], scale = FALSE)

h_alpha <- 1
h_beta <- (p-5)/5
# g <- max(n, p)
g <- 1/2
# h <- 5/p
h <- c(h_alpha, h_beta)

hyper_par <- list(n = nrow(X),
                  p = ncol(X),
                  g = g,
                  h = h,
                  X = X,
                  y = y,
                  # XtX = t(X) %*% X,
                  yty = sum(y^2), 
                  ytX = t(y) %*% X,
                  g_prior_type = "ind")

list.save(hyper_par, settings$file[5])

remove(hyper_par)





n <- nrow(X)
p <- ncol(X)
X <- scale(X)
y <- scale(ys[,2], scale = FALSE)

h_alpha <- 1
h_beta <- (p-5)/5
# g <- max(n, p)
g <- 1/4
# h <- 5/p
h <- c(h_alpha, h_beta)

hyper_par <- list(n = nrow(X),
                  p = ncol(X),
                  g = g,
                  h = h,
                  X = X,
                  y = y,
                  # XtX = t(X) %*% X,
                  yty = sum(y^2), 
                  ytX = t(y) %*% X,
                  g_prior_type = "ind")

list.save(hyper_par, settings$file[6])

remove(hyper_par)




n <- nrow(X)
p <- ncol(X)
X <- scale(X)
y <- scale(ys[,3], scale = FALSE)

h_alpha <- 1
h_beta <- (p-5)/5
# g <- max(n, p)
g <- 1/4
# h <- 5/p
h <- c(h_alpha, h_beta)

hyper_par <- list(n = nrow(X),
                  p = ncol(X),
                  g = g,
                  h = h,
                  X = X,
                  y = y,
                  # XtX = t(X) %*% X,
                  yty = sum(y^2), 
                  ytX = t(y) %*% X,
                  g_prior_type = "ind")

list.save(hyper_par, settings$file[7])

remove(hyper_par)
