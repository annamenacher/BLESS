##############
### Script ###
##############

# Run Gibbs sampler for BLESS model (= Bayesian Lesion Estimation for a Structured Spike-and-Slab Prior) with MCAR prior over sparsity parameters.
# (Code for Simulation Study)

#################
### Libraries ###
#################

library(LaplacesDemon)
library(brglm2)
library(mvtnorm)
library(truncnorm)
library(hutils)

#################
### Constants ###
#################

# path: data folder
path_data = ''
# path: folder with output from parameter estimation with BLESS (backwards DPE results)
path_BLESS = ''
# path: output folder
path = ''
# path: folder with true parameter estimates (if initialization with true coefficients is chosen)
path_truth = ''
# path: folder with Firth parameter estimates (if initialization with Firth estimated coefficients is chosen)
path_Firth = ''
# path: output folder
path = ''

# Load true values of coefficients. 
# Base rate intensity (lambda = {1,2,3}).
lambda = 3
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

# Data 
# Sample size (N = {500,1000,5000}).
N = 1000
# Number of dataset to perform analysis on.
dataset = 1
# Number of voxels.
M = 2500
# Dimension of 2D slice.
dimension =  dim1 = dim2 = sqrt(M)
# Number of covariates.
P = 2
# Initialization of variational infrence algorithm (options: Firth, truth).
init = 'truth'

# Hyperparameter: spike variance (plot marginal posterior of gamma list under prior v0=0 & identify plateau in marginal distribtuion /
# stablization of parameters)
# Index of best spike variance
v0_idx = 12
# Range of spike variances over which DPE was run (in simulation studies: range of 15 spike variances over exp(-20) to exp(-1)). 
#v0 = as.numeric(unlist(read.csv(paste0(path_BLESS, 'v0.csv'))[,2]))
v0 = exp(seq(-20, -1, length.out = 15))
v0 = v0[v0_idx]
# Hyperparameter: slab variance
v1 = 10
# Hyperparameter: degrees of freedom for Wishart distribution prior over precision matrix of MCAR prior over sparsity parameters theta.
v = P
# Set hyperparameter standard deviation of spatially-varying intercept to large value.
sigma_beta0 = 10
# Start of Gibbs sample after initialization.
n_start = 2

# Data
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)
Y = data.matrix(read.csv(paste0(path_data, "Y", dataset, ".csv"), header = T)[,2:(M+1)])
# Number of lesion occurences over sample size at voxel s_j.
N1 = apply(Y, 2, sum)
# Number of no lesions occurrences over sample size at voxel s_j.
N0 = N - N1

# Initialization
if(init == 'Firth'){
  params0 = list()
  x = as.numeric(unlist(read.csv(paste0(path_Firth, 'beta0.csv'))[dataset,2:(M+1)]))
  x[is.na(x)] = 0
  params0$beta0 = x
  x = rbind2(as.numeric(unlist(read.csv(paste0(path_Firth, 'Beta1.csv'))[dataset,2:(M+1)])),
             as.numeric(unlist(read.csv(paste0(path_Firth, 'Beta2.csv'))[dataset,2:(M+1)])))
  x[is.na(x)] = 0
  params0$beta = x
  img = matrix(0, dim1, dim2)
  img[1:dim1, (dim2/2 + 1):dim2] = 1
  params0$gamma = matrix(0, P, M)
  params0$gamma[1,] = matrix(as.vector(img))
  img = matrix(0,	dim1, dim2)
  img[(dim1/2 + 1):dim1, 1:(dim2/2)] = 1
  params0$gamma[2,] = matrix(as.vector(img))
  params0$theta = matrix(0, P, M)
  params0$theta[1,] = rep(0.5, M)
  params0$theta[2,] = rep(0.25, M)
  params0$Sigma_Inv = diag(P)
  params0$xi = matrix(params0$theta^2, nrow = P, ncol = M)
  params0$v0 = v0
  params0$v1 = v1
  params0$z = (cbind(rep(1, N), X) %*% rbind(params0$beta0, params0$beta))
  params0$sigma_beta0 = sigma_beta0
}

if(init=='truth'){
  params0 = list()
  params0$beta0 = as.vector(truth$beta0)
  params0$beta = matrix(0, P, M)
  params0$beta[1,] = as.vector(truth$beta1)
  params0$beta[2,] = as.vector(truth$beta2)
  img = matrix(0, dim1, dim2)
  img[1:dim1, (dim2/2 + 1):dim2] = 1
  params0$gamma = matrix(0, P, M)
  params0$gamma[1,] = matrix(as.vector(img))
  img = matrix(0,	dim1, dim2)
  img[(dim1/2 + 1):dim1, 1:(dim2/2)] = 1
  params0$gamma[2,] = matrix(as.vector(img))
  params0$theta = matrix(0, P, M)
  params0$theta[1,] = rep(0.5, M)
  params0$theta[2,] = rep(0.25, M)
  params0$Sigma_Inv = diag(P)
  params0$xi = matrix(params0$theta^2, nrow = P, ncol = M)
  params0$v0 = v0
  params0$v1 = v1
  params0$z = (cbind(rep(1, N), X) %*% rbind(params0$beta0, params0$beta))
  params0$sigma_beta0 = sigma_beta0
}

# Number of iterations of Gibbs sampler & burn-in
n_iter = 15000
burn_in = 5000

#################
### Functions ###
#################

# Logit function.
logit = function(x){
  x = log(x/(1-x))
  return(x)
}

# Function for Gibbs sampler
gibbs_sampler_probit_regression_structured_ss = function(X, y, dim1, dim2, params0, n_iter, burn_in){

  # Inputs:
  # - X: input data (scalar covariates)
  # - y: output data (image -> lesion masks)
  # - dim1: dimension1 of 3D MRI scan
  # - dim2: dimension2 of 3D MRI scan
  # - params0: list containing intial parameter values
  # - n_iter: number of MCMC iterations
  # - burn_in: burn-in iterations
  
  # Logistic function (this specification avoids numerical issues)
  logistic = function(x){
    x = ifelse(x >= 0, 1 / ( 1 + exp(-x) ), exp(x) / ( 1 + exp(x) ))
    return(x)
  }
  
  # Function to calculate lambda(xi) within the Taylor approximation for acquiring updates of 
  # the sparsity parameter theta. 
  lambda_xi_func = function(x){
    x = (1/(2*x))*(logistic(x) - (1/2))
    return(x)
  }
  
  # Function for deriving indices of adjacent neighbors of every single voxel location for 
  # spatial MCAR prior (a voxel location is considered adjacent if they share the same face).
  adjacency_matrix = function(dim1, dim2, dim3){

    # Input:
    # - dim1: dimension1 of 3D MRI scan
    # - dim2: dimension2 of 3D MRI scan
    # - dim3: dimension3 of 3D MRI scan
    
    if(missing(dim3)){
      A = data.frame(x = integer(), y = integer())
      ind = 1:(dim1*dim2)
      conv = as.vector(matrix(1:(dim1*dim2), dim1, dim2, byrow = T))
      
      for (i in 1:(dim1 * dim2)){
        up = i - dim2
        down = i + dim2
        left = i - 1
        right = i + 1
        if (up > 0){
          A = rbind(A, c(i, up))
        }
        if (down <= (dim1 * dim2)){
          A = rbind(A, c(i, down))
        }
        if (left %% dim2 != 0){
          A = rbind(A, c(i, left))
        }
        if (i %% dim2 != 0){
          A = rbind(A, c(i, right))
        }
      }
      colnames(A) = c('x', 'y')
      Ax = numeric(length(A$x))
      Ay = numeric(length(A$y))
      for(i in 1:length(A$x)){
        Ax[i] = ind[which(conv == A$x[i], arr.ind = T)]
        Ay[i] = ind[which(conv == A$y[i], arr.ind = T)]
      }
      A$x = Ax
      A$y = Ay
    } else{
      A_2D = data.frame(x = integer(), y = integer())
      ind = 1:(dim1*dim2*dim3)
      conv = as.vector(aperm(array(1:(dim1*dim2*dim3), dim = c(dim2, dim1, dim3)), perm = c(2, 1, 3)))
      
      for (i in 1:(dim1 * dim2)){
        up = i - dim2
        down = i + dim2
        left = i - 1
        right = i + 1
        if (up > 0){
          A_2D = rbind(A_2D, c(i, up))
        }
        if (down <= (dim1 * dim2)){
          A_2D = rbind(A_2D, c(i, down))
        }
        if (left %% dim2 != 0){
          A_2D = rbind(A_2D, c(i, left))
        }
        if (i %% dim2 != 0){
          A_2D = rbind(A_2D, c(i, right))
        }
      }
      colnames(A_2D) = c('x', 'y')
      A = data.frame(x = integer(), y = integer())
      for (k in 0:(dim3-1)) {
        A = rbind(A, (A_2D + (k*dim1*dim2)))
      }
      for(i in 1:(dim1*dim2*dim3)){
        bottom = i - dim1*dim2
        top = i + dim1*dim2
        if(bottom > 0){
          A = rbind(A, c(i, bottom))
        }
        if(top <= (dim1*dim2*dim3)){
          A = rbind(A, c(i, top))
        }
      }
      Ax = numeric(length(A$x))
      Ay = numeric(length(A$y))
      for(i in 1:length(A$x)){
        Ax[i] = ind[conv == A$x[i]]
        Ay[i] = ind[conv == A$y[i]]
      }
      A$x = Ax
      A$y = Ay
    }
    return(A)
  }
  
  # Indices of adjacency matrix of 2D lattice.
  A = adjacency_matrix(dim1, dim2)
  
  # Function for deriving number of neighbors of every single voxel location for spatial MCAR prior.
  n_neighbors = function(dim1, dim2, dim3){

    # Input:
    # - dim1: dimension1 of 3D MRI scan
    # - dim2: dimension2 of 3D MRI scan
    # - dim3: dimension3 of 3D MRI scan
    
    if(missing(dim3)){
      if (dim1 < 3 | dim2 < 3){ 
        stop("Image dimensions need to be greater than 2!")
      }
      n_sj = matrix(4, nrow = dim1, ncol = dim2)
      n_sj[1,1] = n_sj[1,dim2] = n_sj[dim1,1] = n_sj[dim1,dim2] = 2
      n_sj[2:(dim1-1),1] = n_sj[2:(dim1-1),dim2] = n_sj[1,2:(dim2-1)] = n_sj[dim1,2:(dim2-1)] = 3
      n_sj = as.vector(n_sj)
    } else{
      if (dim1 < 3 | dim2 < 3 | dim3 < 3){ 
        stop("Image dimensions need to be greater than 2!")
      }
      n_sj = array(6, c(dim1, dim2, dim3))
      n_sj[1,1,1] = n_sj[1,dim2,1] = n_sj[dim1,1,1] = n_sj[dim1,dim2,1] = n_sj[1,1,dim3] = n_sj[1,dim2,dim3] = n_sj[dim1,1,dim3] = n_sj[dim1,dim2,dim3] = 3
      n_sj[2:(dim1-1),1,1] = n_sj[2:(dim1-1),dim2,1] = n_sj[dim1,2:(dim2-1),1] = n_sj[1,2:(dim2-1),1]  = 4
      n_sj[2:(dim1-1),1,dim3] = n_sj[2:(dim1-1),dim2,dim3] = n_sj[dim1,2:(dim2-1),dim3] = n_sj[1,2:(dim2-1),dim3] = 4
      n_sj[1,1,2:(dim3-1)] = n_sj[1,dim2,2:(dim3-1)] = n_sj[dim1,1,2:(dim3-1)] = n_sj[dim1,dim2,2:(dim3-1)] = 4
      n_sj[2:(dim1-1), 2:(dim2-1),1] = n_sj[2:(dim1-1), 2:(dim2-1),dim3] = n_sj[1, 2:(dim2-1),2:(dim3-1)] = 5
      n_sj[dim1, 2:(dim2-1),2:(dim3-1)] = n_sj[2:(dim1-1), 1,2:(dim3-1)] = n_sj[2:(dim1-1), dim2,2:(dim3-1)] = 5
      n_sj = as.vector(n_sj)
    }
    return(n_sj)
    
  }
 
  # Number of neighbors of 2D lattice.
  n_sj = n_neighbors(dim1, dim2)
  
  # Function to acquire indices of upper triangular of adjacency matrix.
  upper_triangular = function(A, M){

    # Input: 
    # - A: spatial adjacency matrix
    # - M: total number of voxels
    
    A = A[order(A$x),]
    UT = data.frame(x = integer(), y = integer())
    K = dim(A)[1]
    for(k in 1:K){
      idx = A[k,]
      if(idx[1] < idx[2]){
        UT = rbind(UT, idx)
      }
    }
    return(UT)
    
  }
  
  # Get indices of upper triangular of adjacency matrix for 2D lattice. 
  A = upper_triangular(A, dim1*dim2)
  
  # Make a list with every voxel location (vector: 1 to M) and its neighbor indices 
  # (for example: voxel 1 in a 5x5 image will have neighbors at voxels 2 & 6, 
  # and has therefore 2 neighbors.
  list_ind = list()
  for(j in 1:M){
    list_ind[[j]] = c(A$y[which(A$x == j, arr.ind = T)], A$x[which(A$y == j, arr.ind = T)])
  }
  
  # Initialize parameters & hyperparameters.
  sigma_beta0 = params0$sigma_beta0
  v0 = params0$v0
  v1 = params0$v1
  gamma = params0$gamma
  theta = params0$theta
  Sigma_Inv = params0$Sigma_Inv
  beta = params0$beta
  beta0 = params0$beta0
  z = params0$z
  xi = params0$xi
  
  # Matrix storing samples of parameters
  beta0_chain = matrix(NA, nrow = (n_iter-burn_in), ncol = M)
  beta1_chain = matrix(NA, nrow = (n_iter-burn_in), ncol = M)
  beta2_chain = matrix(NA, nrow = (n_iter-burn_in), ncol = M)
  theta1_chain = matrix(NA, nrow = (n_iter-burn_in), ncol = M)
  theta2_chain = matrix(NA, nrow = (n_iter-burn_in), ncol = M)
  Sigma_Inv_chain = matrix(NA, nrow = (n_iter-burn_in), ncol = P*P)
  gamma1_chain = matrix(NA, nrow = (n_iter-burn_in), ncol = M)
  gamma2_chain = matrix(NA, nrow = (n_iter-burn_in), ncol = M)
  xi1_chain = matrix(NA, nrow = (n_iter-burn_in), ncol = M)
  xi2_chain = matrix(NA, nrow = (n_iter-burn_in), ncol = M)
 
  # Run chain for n_iter samples.
  for (t in n_start:n_iter) {
      
      sum_si_sj = matrix(NA, nrow = P, ncol = M)
    
      for(j in 1:M){
      # Draw latent variable z from its full conditional: z | beta, y, X
      mu_z_j = X %*% matrix(beta[,j],P,1) + c(beta0[j])
      z_j = z[,j]
      y_j = Y[,j]
      if(N1[j] == 0){
        z_j[y_j == 0] = rtruncnorm(N0[j], mean = mu_z_j[y_j == 0], sd = 1, a = -Inf, b = 0)
      } else{
        z_j[y_j == 0] = rtruncnorm(N0[j], mean = mu_z_j[y_j == 0], sd = 1, a = -Inf, b = 0)
        z_j[y_j == 1] = rtruncnorm(N1[j], mean = mu_z_j[y_j == 1], sd = 1, a = 0, b = Inf)
      }
      z[,j] = z_j
      
      # Compute posterior variance of beta
      precision_beta = diag(c(gamma[,j]*(1/v1) + (1-gamma[,j])*(1/v0)), P)
      V_beta = solve(precision_beta + crossprod(X, X))
      # Compute posterior mean of beta.
      mean_beta = V_beta %*%crossprod(X, (z[,j] - beta0[j]))
      # Draw variable beta from its full conditional: beta | z, X
      beta[,j] = c(rmvnorm(1, mean_beta, V_beta))
      
      # Draw variable beta0 from its full conditional: beta0 | z, X
      V_beta0 =  (1/(N + 1/sigma_beta0^2))
      mean_beta0 = V_beta0*(sum(z[,j] - X%*%beta[,j]))
      beta0[j] = rnorm(1, mean_beta0, sqrt(V_beta0))
      
      # Calculate posterior mean of Bernoulli.
      l0 = dnorm(c(beta[,j]), 0, sqrt(v0))*logistic(1-c(theta[,j]))
      l1 = dnorm(c(beta[,j]), 0, sqrt(v1))*logistic(c(theta[,j]))
      
      # Draw variable gamma from its full conditional: gamma | theta, beta
      gamma[,j] = rbinom((l1/(l0 + l1)), 1, (l1/(l0 + l1)))
      
      # Draw variable gamma from its full conditional: theta | Sigma_Inv, xi
      ind = c(A$y[which(A$x == j, arr.ind = T)], A$x[which(A$y == j, arr.ind = T)])
      sum_si_sj[,j] = matrix(apply(matrix(theta[,ind], nrow = P), 1, sum), nrow = P)
      var_theta = solve(n_sj[j]*Sigma_Inv + 2*diag(c(lambda_xi_func(xi[,j])), P))
      mu_theta = var_theta%*%(Sigma_Inv%*%matrix(sum_si_sj[,j], P, 1) + gamma[,j] - (1/2))
      theta[,j] = c(rmvnorm(1, mu_theta, var_theta))
      xi[,j] = sqrt(diag(var_theta) + theta[,j]^2)
      }
    
    # Draw variable gamma from its full conditional: Sigma_Inv | theta
    term = matrix((theta[,A$x] - theta[,A$y]), P)%*%t(matrix((theta[,A$x] - theta[,A$y]), P))
      
    v = M - 1 
    I = solve(diag(P) + term)
    Sigma_Inv = matrix(rWishart(1, v, I), P, P)
    
    # Don't save burn-in samples (lessen memory cost)
    if(t > (burn_in)){
        # Store the beta draws.
        beta0_chain[(t-burn_in), ] = as.vector(beta0)
        beta1_chain[(t-burn_in), ] = as.vector(beta[1,])
        beta2_chain[(t-burn_in), ] = as.vector(beta[2,])
        theta1_chain[(t-burn_in), ] = as.vector(theta[1,])
        theta2_chain[(t-burn_in), ] = as.vector(theta[2,])
        Sigma_Inv_chain[(t-burn_in), ] = as.vector(Sigma_Inv)
        xi1_chain[(t-burn_in), ] = as.vector(xi[1,])
        xi2_chain[(t-burn_in), ] = as.vector(xi[2,])
        gamma1_chain[(t-burn_in), ] = as.vector(gamma[1,])
        gamma2_chain[(t-burn_in), ] = as.vector(gamma[2,])
    }
    # Save parameters every 100 samples.
    if(t%%100 == 0){
    write.csv(beta0_chain, paste0(path, "beta0_chain.csv"))
    write.csv(beta1_chain, paste0(path, "beta1_chain.csv"))
    write.csv(beta2_chain, paste0(path, "beta2_chain.csv"))
    write.csv(gamma1_chain, paste0(path, "gamma1_chain.csv"))
    write.csv(gamma2_chain, paste0(path, "gamma2_chain.csv"))
    write.csv(Sigma_Inv_chain, paste0(path, "Sigma_Inv_chain.csv"))
    write.csv(theta1_chain, paste0(path, "theta1_chain.csv"))
    write.csv(theta2_chain, paste0(path, "theta2_chain.csv"))
    write.csv(xi1_chain, paste0(path, "xi1_chain.csv"))
    write.csv(xi2_chain, paste0(path, "xi2_chain.csv"))
    }
  }
  
  # Save final chain. 
  result = list()
  result$beta0_chain = beta0_chain
  result$beta1_chain = beta1_chain
  result$beta2_chain = beta2_chain
  result$theta1_chain = theta1_chain
  result$theta2_chain = theta2_chain
  result$gamma1_chain = gamma1_chain
  result$gamma2_chain = gamma2_chain
  result$xi1_chain = xi1_chain
  result$xi2_chain = xi2_chain
  result$Sigma_Inv_chain = Sigma_Inv_chain
  
  write.csv(result$beta0_chain, paste0(path, "beta0_chain.csv"))
  write.csv(result$beta1_chain, paste0(path, "beta1_chain.csv"))
  write.csv(result$beta2_chain, paste0(path, "beta2_chain.csv"))
  write.csv(result$gamma1_chain, paste0(path, "gamma1_chain.csv"))
  write.csv(result$gamma2_chain, paste0(path, "gamma2_chain.csv"))
  write.csv(result$Sigma_Inv_chain, paste0(path, "Sigma_Inv_chain.csv"))
  write.csv(result$theta1_chain, paste0(path, "theta1_chain.csv"))
  write.csv(result$theta2_chain, paste0(path, "theta2_chain.csv"))
  write.csv(result$xi1_chain, paste0(path, "xi1_chain.csv"))
  write.csv(result$xi2_chain, paste0(path, "xi2_chain.csv"))
  return(result)
}

# Run Gibbs sampler for BLESS.
result = gibbs_sampler_probit_regression_structured_ss(X, Y, dim1, dim2, params0, n_iter, burn_in)

