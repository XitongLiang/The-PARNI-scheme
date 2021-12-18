load("real dataset/dataset/SNR/cfw.RData")



colnames(pheno)
sacwt <- pheno[,1]

pheno_new <- pheno[!is.na(sacwt),]
geno_new <- geno[!is.na(sacwt),]

testis <- pheno_new[,6]
pheno_new <- pheno_new[!is.na(testis),]
geno_new <- geno_new[!is.na(testis),]


X1 <- pheno_new[,1]
Z <- cbind(1, X1)
X <- geno_new

n <- nrow(X)
p <- ncol(X)

y <- (diag(n) - Z %*% solve(t(Z) %*% Z) %*% t(Z)) %*% pheno_new[,6]
X <- (diag(n) - Z %*% solve(t(Z) %*% Z) %*% t(Z)) %*% X
# X
X <- scale(X)
colnames(X) <- NULL
colnames(y) <- NULL
y <- scale(y, scale = FALSE)

rm(Z, geno, geno_new, gwscan.bvsr, gwscan.gemma, map, sacwt, testis, X1, pheno, pheno_new)




# h_alpha <- 1
# h_beta <- (p-5)/5
# g <- max(n, p)
g <- 1/4
h <- 5/p

hyper_par <- list(n = nrow(X),
                  p = ncol(X),
                  g = g,
                  # h = h,
                  # h_alpha = h_alpha,
                  # h_beta = h_beta,
                  # h = c(h_alpha, h_beta),
                  h = h,
                  X = X,
                  y = y,
                  yty = sum(y^2), 
                  ytX = t(y) %*% X,
                  g_prior_type = "ind",
                  diag_V = colSums(X^2) + 1/g)


list.save(hyper_par, settings$file[8])

remove(hyper_par)
# hyper_par <- list.load("real dataset/dataset/SNP.RData")

