#### Synthetic data experiment (DGPs) ####
X <- matrix(c(rnorm(100,1,0.05),
              c(rep(1,50),rep(0,50)),
              sin(seq(0.1,10,0.1))) , ncol = 3)
s <- matrix(c(seq(0.1,10,0.1),
              runif(100,0,1)),ncol = 2)
beta <- matrix(c(0.5,1,2,
                 1,2,1,
                 2,0.5,0), ncol = 3) # K=3
eta <- X %*% beta # + ( s[,1]^(1/5) - s[,2]^(1/3) + 0.5 * sin(s[,1]) + 0.5 * cos(s[,2]) )  # spacial covariates
pi <- mclust::softmax(eta)
theta <- matrix(runif(30,0,10),ncol = 3) # d = 10
z <- apply(pi , 1, function(prob) sample.int(3, size = 1, prob = prob) )
Y <- t(sapply(z, function(k) rdirichlet(1, theta[, k])))

#### Real GRUMP data experiment ####
# load(file = "./relative_abund.Rdata")
# 
# df_wide_mat <- as.matrix(df_wide[,-1])



#### Fitting the model ####
# On GRUMP data
# res <- EM_dirichlet(df_wide_mat+1e-10,nclust = 2,max_iter = 100)
# res3 <- EM_dirichlet(df_wide_mat+1e-10,nclust = 3,max_iter = 100)

# On simulated data
res <- EM_dirichlet(Y,nclust = 3,max_iter = 100, regress = FALSE) # without regression
t(theta)
res$theta # theta estimated well

res_regress_X <- EM_dirichlet(Y = Y, X = X, nclust = 3, max_iter = 100, regress = TRUE)
t(theta)
res_regress_X$theta # theta estimated well again.
plot(pi[,3],type = "l")
lines(res_regress_X$pi[,2],col="red")
plot(pi[,1],type = "l")
lines(res_regress_X$pi[,1],col="red")



res_regress <- EM_dirichlet(Y = Y, X = X, nclust = 3, max_iter = 500, regress = TRUE)
res_regress$beta
beta
plot(pi[,3],type = "l")
lines(res_regress$pi[,3],col="red")

