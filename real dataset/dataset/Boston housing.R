library(mlbench)
data(BostonHousing)

y <- log(BostonHousing$medv)
# y <- BostonHousing$medv
# X <- BostonHousing[,c(-4, -14)]

X <- model.matrix(lm(y~(.)^2-1, data = BostonHousing[,c(-14)]))[,-4]
X <- cbind(X, BostonHousing[,c(-4, -14)]^2)

colnames(X) <- NULL
# model.matrix(~(colnames(X[,-4]))^2, X)
# X <- as.numeric(X)
# X <- as.ma3trix(X)

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

list.save(hyper_par, settings$file[3])

remove(hyper_par)

