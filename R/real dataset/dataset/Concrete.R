concrete <- read.csv("real dataset/dataset/Concrete_Data/Concrete_Data.csv", header = TRUE)

y <- as.vector(concrete[,9])
X <- as.matrix(concrete[,-9])
X <- cbind(X, log(X[,c(1,4,6,7,8)]))
colnames(X) <- NULL
X <- model.matrix(lm(y~(.)^2-1, data = data.frame(X)))

n <- nrow(X)
p <- ncol(X)
X <- scale(X)
y <- scale(y, scale = TRUE)
h <- 5/100 #  #
g <- 100

hyper_par <- list(n = nrow(X),
                  p = ncol(X),
                  g = g,
                  h = h,
                  X = X,
                  y = y,
                  yty = sum(y^2), 
                  ytX = t(y) %*% X,
                  g_prior_type = "ind",
                  diag_V = colSums(X^2) + 1/g)

list.save(hyper_par, settings$file[2])

remove(hyper_par)