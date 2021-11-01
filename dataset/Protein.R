# Protein activity data

# install.packages("BAS")
library(BAS)

data(protein)
# protein

y <- protein$prot.act1

protein$buf <- factor(protein$buf)
protein$ra <- factor(protein$ra)
protein$det <- factor(protein$det)

X <- model.matrix(lm(y~.-1, data = protein[,1:8]))
X <- model.matrix(lm(y~(.)^2-1, data = data.frame(X)))

X <- X[, which(colSums(X)!=0)]
colnames(X) <- NULL

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
                  # XtX = t(X) %*% X,
                  yty = sum(y^2), 
                  ytX = t(y) %*% X,
                  g_prior_type = "ind")

list.save(hyper_par, settings$file[4])

remove(hyper_par)