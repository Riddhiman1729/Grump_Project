#### Synthetic data experiment (DGPs) ####
# Only X
X <- matrix(c(cos(seq(0,2*3.14,0.05)), # Represent -temperature
              -cos(seq(0,4*3.14,0.1))), # Represent salinity
            ncol = 2)
beta <- matrix(c(2,0, # prototype 1 sensitive to temperature
                 0,2, # prototype 2 sensitive to salinity
                 2,2), # prototype 3 sensitive to both
               ncol = 3) # K=3
eta <- X %*% beta
pi <- mclust::softmax(eta)
d = 10
theta <- matrix(runif(d * ncol(beta), 0 , d),ncol = ncol(beta)) # d = 10
z <- apply(pi , 1, function(prob) sample.int(length(prob), size = 1, prob = prob) )
Y <- t(sapply(z, function(k) rdirichlet(1, theta[, k])))

# Both s and X
s <- matrix(c(seq(0.1,10,0.1),
              runif(100,0,1)),ncol = 2)
g1 <-  0.1*sin(s[,1]) + 0.8*s[,2]^2 - 0.3*s[,1]*s[,2]
g2 <- -0.7*cos(s[,2]) + 0.2*s[,2]^3
g3 <-  0.2*sin(s[,1] + 2*s[,2]) + exp(-s[,1]^2)
eta_full <- X %*% beta + cbind(g1,g2,g3)  # spacial covariates
pi_full <- mclust::softmax(eta_full)
theta_full <- matrix(runif(30,0,10),ncol = 3) # d = 10
z_full <- apply(pi_full , 1, function(prob) sample.int(3, size = 1, prob = prob) )
Y_full <- t(sapply(z_full, function(k) rdirichlet(1, theta_full[, k])))



#### Real GRUMP data experiment ####
# load(file = "./relative_abund.Rdata")
# 
# df_wide_mat <- as.matrix(df_wide[,-1])




#### Fitting the model ####
# On GRUMP data
# res <- EM_dirichlet(df_wide_mat+1e-10,nclust = 2,max_iter = 100)
# res3 <- EM_dirichlet(df_wide_mat+1e-10,nclust = 3,max_iter = 100)

# On simulated data
# without regression
res <- EM_dirichlet(Y,nclust = 3,max_iter = 100, regress = FALSE) 
t(theta)
res$theta # theta are well estimated 

# regression with X
res_regress <- EM_dirichlet(Y = Y, X = X, nclust = 3, max_iter = 500, regress = TRUE)
res_regress$beta
beta
for(i in 1:3){
  plot(pi[,i],type = "l",ylim = c(0,1))
  lines(res_regress$pi[,i],col="red") # pi are well estimated but not the beta
}
res_regress$theta 
theta # theta are also well estimated

# regression with X and s
res_regress_full <- EM_dirichlet(Y = Y_full, X = X, s = s, nclust = 3, max_iter = 500, regress = TRUE)
plot(pi_full[,3],type = "l",ylim = c(0,1))
lines(res_regress_full$pi[,2],col="blue") # The pi estimated well but not the beta again
res_regress_full$theta
theta_full # theta are well estimated
