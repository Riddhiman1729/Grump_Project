#install.packages('gtools')
library('gtools')

#'@param dimdat.y the number of dimensions of the compostions 
#'@param dimdat.x the number of dimensions of the covariates 
#'@param nsamples the number of samples to be generated 
#'@param nclusters the number of clusters
#'@param left_tol tolerance below which an observation is deemed zero.
#'@return a list of matrices containing the response, the covariates, the true theta, beta and cluster probabilities



data_gen_fn <- function(dimdat.y, dimdat.x, nsamples, nclusters, left_tol){
  theta.mat = matrix(runif(dimdat.y * nclusters,2,3),nrow = dimdat.y, ncol = nclusters)
  beta.mat = matrix(rnorm((dimdat.x + 1)*nclusters,2,1),(dimdat.x + 1), nclusters)
  beta.mat[1,] = beta.mat[1,]+1
  beta.mat[,1] = 0
  X=matrix(NA,nsamples,dimdat.x)
  X[,1] = sin(seq(from=1,to=100,length=nsamples))
  X[,2] = c(rep(1,floor(nsamples/2)),rep(0,nsamples-floor(nsamples/2)))
  
  if(dimdat.x>2){
    for(dd in 3:dimdat.x){
      X[,dd]=rnorm(nsamples,0,1)
    }
    
  }
  
  Xa=cbind(1,X)
  
  ymat=matrix(NA,nsamples,dimdat.y)
  
  for(ii in 1:nsamples){
    cluster_id = sample(1:nclusters,1,F)
    ymat[ii, ] = rdirichlet(1, theta.mat[,cluster_id])
  }
  
  ymat=ymat/rowSums(ymat)
  
  ymat=ymat*(ymat>left_tol)
  
  ymat=ymat/rowSums(ymat)
  
  cluster.probs = exp(Xa%*%beta.mat)/rowSums(exp(Xa%*%beta.mat))
  
  return(list(ymat = ymat,
              X = X,
              cluster_probs = cluster.probs,
              beta_mat = beta.mat, 
              theta_mat = theta.mat
            ))  
}



