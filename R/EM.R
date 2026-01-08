# Generate a mixture of dirichlet dataset
library(MCMCpack)

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
  alpha
}

# Helper: invert digamma(x) - Ref: Minka TP (2012)
digamma_inv <- function(x) {
  y <- ifelse(x >= -2.22, exp(x) + 0.5, -1 / (x + digamma(1)))
  for (i in 1:5) y <- y - (digamma(y) - x) / trigamma(y)
  
  return(y)
}

# Estep: Calculate the responsibilities
e_step <- function(X, pi , theta){
  n <- nrow(X)
  K <- length(pi)
  
  log_w <- matrix(nrow = n,ncol = K)
  
  for(k in 1:K){
    log_fk <- log_dirichlet(X, theta[[k]])
    log_w[,k] <- log(pi[k]) + log_fk
  }
  log_denom <- log(rowSums(exp(log_w)))
  resp <- exp(log_w - log_denom)
  
  return(resp)
}


# Mstep: update the parameters
m_step <- function(X, resp, theta){
  n <- nrow(X)
  K <- ncol(resp)
  
  # Update p
  Nk <- colSums(resp)
  pi_new <- Nk/n
  
  # Update theta
  theta_new <- vector("list", K)
  for(k in 1:K){
    theta_new[[k]] <- estimate_dirichlet(X, 
                                         resp[,k],
                                         maxiter = 1000,
                                         tol = 1e-6)
  }
  
  return(list(pi = pi_new, theta = theta_new))
}


# EM function
EM_dirichlet <- function(X,
                         nclust = 2,
                         max_iter = 1000,
                         tol = 1e-6,
                         verbose = FALSE){
  n <- nrow(X)
  d <- ncol(X)
  K <- nclust
  
  # Initialization
  pi <- rep(1/K , K)
  theta <- vector("list", K)
  for(k in 1:K){
    theta[[k]] <- sample(1:10,d)
  }
  # loglik_trace <- rep(-Inf, max_iter)
  
  # EM algorithm
  for(iter in 1:max_iter){
    # E step
    resp <- e_step(X, pi, theta)
    
    # M step
    mstep_res <- m_step(X, resp , theta)
    pi_new <- mstep_res$pi
    theta_new <- mstep_res$theta
    
    # Convergence check
    # if(iter > 1 && abs(loglik_trace[iter] - loglik_trace[iter - 1]) < tol){
    #   loglik_trace <- loglik_trace[1:iter]
    #   break
    # }
    
    # Update parameters
    pi <- pi_new
    theta <- theta_new
  }
  
  return( list(
    pi = pi,
    theta = theta,
    resp = resp
    #log_likeli <- log_likeli_trace
  ))
}




# Test
theta1 <- c(2,5,3,4)
theta2 <- c(8,2,1,3)
X <- rbind(rdirichlet(1000,theta1),
           rdirichlet(500,theta2))

res <- EM_dirichlet(X, nclust = 2, max_iter = 100)
res$pi
res$theta
