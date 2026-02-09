# Generate a mixture of dirichlet dataset
library(MCMCpack)
library(fda)

# Helper: Log dirichlet density function(stable)
log_dirichlet <- function(X, alpha) {
  log_dirichlet_1d <- function(x , alpha){
    lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
  }
  n <- nrow(X) # Check is matrix or vector
  if(is.null(n)){
    return(log_dirichlet_1d(X,alpha))
  }else{
    return(apply(X , 1 , log_dirichlet_1d , alpha = alpha))
  }
}


# Helper: Estimate the dirichlet parameters (with weights)
estimate_dirichlet <- function(X, w, maxiter = 100, tol = 1e-6) {
  n <- nrow(X)
  d <- ncol(X)
  Nk <- sum(w)
  
  log_ybar <- colSums(w * log(X))/Nk 
  alpha <- rep(1, d) 
  
  for (t in 1:maxiter) {
    alpha_old <- alpha
    s <- sum(alpha_old)
    alpha <- digamma_inv( digamma(s) + log_ybar )
    # stop if converged
    if (max(abs(alpha - alpha_old)) < tol) break
  }
  return(alpha)
}

# Helper: invert digamma(x) - Ref: Minka TP (2012)
digamma_inv <- function(x) {
  y <- ifelse(x >= -2.22, exp(x) + 0.5, -1 / (x + digamma(1)))
  for (i in 1:5) y <- y - (digamma(y) - x) / trigamma(y)
  
  return(y)
}


# Helper: solve multinomial glm
solve_multinom <- function(y, X, lambda = 0){
  fit <- glmnet::glmnet(x = X,
                        y = y,
                        family = "multinomial",
                        intercept = TRUE)
  coefs = glmnet::coef.glmnet(fit, s = lambda)
  return(as.matrix(do.call(cbind, coefs)))
}


# Helper: base spline
solve_spline <- function(s){
  if(ncol(s) == 1){
    basis <- create.bspline.basis(range(s), nbasis = 4, norder = 3)
    B <- eval.basis(s,basis)
    B <- scale(B, center=TRUE, scale=FALSE)
    
  }else if(ncol(s) == 2){
    s1 <- s[,1]
    s2 <- s[,2]
    basis1 <- create.bspline.basis(range(s1), nbasis=4, norder=3)
    basis2 <- create.bspline.basis(range(s2), nbasis=4, norder=3)
    
    B1 <- eval.basis(s1, basis1)
    B2 <- eval.basis(s2, basis2)
    # center columns to avoid intercept-like duplication
    B1 <- scale(B1, center=TRUE, scale=FALSE)
    B2 <- scale(B2, center=TRUE, scale=FALSE)
    
    B <- as.matrix(do.call(cbind, lapply(1:ncol(B1), function(a) B1[,a] * B2)))
  }
  return(B)
}



##### glmnet test (To be deleted) #####
# X = matrix(c(1,2,3,5,2,1,3,1),ncol = 2)
# y = matrix(c(0.3,0.6,0.1,0.2,0.6,0.2,0.4,0.2,0.4,0.3,0.5,0.2),ncol = 3,byrow = TRUE)
# fit <- glmnet(x = X,
#        y = y,
#        family = "multinomial",
#        intercept = TRUE)
# coefs = coef(fit,s = fit$lambda[10])
# as.matrix(do.call(cbind, coefs)) # Get coefficents
# predict(fit, newx = X, s = fit$lambda[10], type = "response") # Make prediction
##### end of test (To be deleted) #####


##### spline version glmnet test (To be deleted) #####
# library(splines) # Choice one
# library(fda) # Choice two
# s <- c(100,101,102,103) # spacial coeff
# # Using splines package
# B_sp <- splines::bs(s, df = 3, degree = 3, intercept = FALSE)
# # Using fda package
# basis <- fda::create.bspline.basis(range(s), nbasis = 6, norder = 4)  # cubic: norder=4
# B_fda <- fda::eval.basis(s, basis)
# B_fda <- scale(B_fda, center = TRUE, scale = FALSE)
##### end of test (To be deleted) #####






# Estep: Calculate the responsibilities
e_step <- function(Y, pi , theta){
  n <- nrow(Y)
  K <- ncol(pi)
  
  log_w <- matrix(nrow = n,ncol = K)
  
  for(k in 1:K){
    log_fk <- log_dirichlet(Y, theta[[k]])
    log_w[,k] <- log(pi[,k]) + log_fk
  }
  
  m <- apply(log_w,1,max)
  log_denom <- m + log(rowSums(exp(log_w-m))) # replace log_denom <- log(rowSums(exp(log_w))) to aviod -log0
  resp <- exp(log_w - log_denom)
  
  return(resp)
}


# Mstep: update the parameters
m_step <- function(Y, X = NULL, s = NULL, resp, theta, regress = FALSE){
  n <- nrow(Y)
  K <- ncol(resp)
  
  # Update pi and beta
  if(!regress){
    Nk <- colSums(resp)
    pi_new <- matrix(rep(Nk/n,n),ncol = K,byrow = TRUE)
    beta = NULL
  }else{
    if(is.null(X)){
      print("X is missing")
    }
    
    if(is.null(s)){
      X_full = X
    }else if(ncol(s) %in% c(1,2)){
      B <- solve_spline(s)
      X_full = cbind(X,B)
    }else{
      print("Spatial covariates dimension error")
    }
    
    # only X corvariates
    beta = solve_multinom(y = resp, X = X_full)
    pi_new_mat = as.matrix(exp(cbind(1,X_full) %*% beta))
    pi_new <- pi_new_mat / rowSums(pi_new_mat)
    # X + f(s)
    # X + f(s1,s2)
  }

  # Update theta
  theta_new <- vector("list", K)
  for(k in 1:K){
    theta_new[[k]] <- estimate_dirichlet(Y, 
                                         resp[,k],
                                         maxiter = 1000,
                                         tol = 1e-6)
  }
  
  return(list(pi = pi_new, beta = beta, theta = theta_new))
}


# EM function
EM_dirichlet <- function(Y,
                         X = NULL,
                         s = NULL,
                         nclust = 2,
                         max_iter = 1000,
                         tol = 1e-6,
                         regress = FALSE,
                         verbose = FALSE){
  n <- nrow(Y)
  d <- ncol(Y)
  K <- nclust
  
  # Initialization
  pi <- matrix(rep(1/K , K * n), nrow = n)
  theta <- vector("list", K)
  for(k in 1:K){
    samp <- sample.int(nrow(Y) , nrow(Y) %/% K )
    theta[[k]] <- colSums(Y[samp,])   # need to think about the initialization
  }
  # loglik_trace <- rep(-Inf, max_iter)
  
  # EM algorithm
  for(iter in 1:max_iter){
    if(iter %% 10 == 0){print(iter)} # print iter(TO BE DELETED)
    # E step
    resp <- e_step(Y, pi, theta)
    
    # M step
    mstep_res <- m_step(Y = Y, X = X, s = s, resp , theta, regress = regress)
    pi_new <- mstep_res$pi
    beta_new <- mstep_res$beta
    theta_new <- mstep_res$theta
    
    # Convergence check
    # if(iter > 1 && abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol){
    #   loglik_trace <- loglik_trace[1:iter]
    #   break
    # }
    
    # Update parameters
    pi <- pi_new
    beta <- beta_new
    theta <- theta_new
  }
  
  return( list(
    pi = pi,
    beta = beta,
    theta = theta,
    resp = resp
    #log_likeli <- log_likeli_trace
  ))
}




# Test
# theta1 <- c(2,5,3,4)
# theta2 <- c(8,2,1,3)
# Y <- rbind(rdirichlet(1000,theta1),
#            rdirichlet(500,theta2))
# 
# res <- EM_dirichlet(Y = Y, nclust = 2, max_iter = 100, regress = FALSE)
# res$pi
# res$theta
