library(fda.usc)

data(tecator)


X_fda <- tecator$absorp.fdata
# View(X_fda)

# X_fda$names
X_data <- X_fda$data

y_t <- tecator$y

# View(cor(X_data))

X <- X_data
y <- y_t$Fat
y <- y-mean(y)
p <- ncol(X)
n <- nrow(X)
X <- t(t(X)-colMeans(X))

hyper_par <- list(y = y,
               X = X,
               p = p,
               n = n,
               yty = sum(y^2),
               ytX = t(y) %*% X,
               h = 5/100,
               g = 100,
               g_prior_type = "ind",
               diag_V = colSums(X^2) + 1/g)

list.save(hyper_par, settings$file[1])

remove(hyper_par)

