##############
### Script ###
##############

# Run BSGLMM (=Bayesian Generalized Linear Mixed Model) (Code for Simulation Study)

#################
### Libraries ###
#################

library(hutils)
library(brglm2)
library(MASS)
library(logOfGamma)

#################
### Constants ###
#################

# Number of dataset to perform analysis on.
sim = 1
# Number of voxels.
M = 2500
# Dimension of 2D slice.
dim1 = dim2 = dimension = 50
# Number of covariates.
P = 2
# Convergence threshold.
eps = 0.0001
# Sample size (N = {500,1000,5000}).
N = 1000
# Base rate intensity (lambda = {1,2,3}).
lambda = 3
# Initialisation of Parameters (Options: Firth, random)
init = 'Firth'

# Path: data
path_data = ""
# Path: output
path = ""
# Path: initialization (if Firth is chosen as option for initial values)
path_Firth = ""

#######################
### BSGLMM function ###
#######################

estimate_BSGLMM=function(X, Y, params0, eps){

  # Inputs:
  # - X: input data (scalar covariates)
  # - Y: output data (image -> lesion masks)
  # - params0: list containing intial parameter values
  # - eps: optimization threshold = difference between log-likelihoods which determines when to stop optimization
  
  
  # Function to calculate sum of adjacent parameter values that are neighbors as 
  # determined by list of neighboring values.
  sum_si_sj_function = function(j){
    x = matrix(apply(matrix(Beta[,list_ind[[j]]], nrow = P), 1, sum), nrow = P)
    return(x)
  }
  
  # Function to update spatially-varying coefficients.
  Beta_function = function(j){
    x = matrix(Sigma_Beta[,j], P, P) %*% (t(X) %*% (expected_Z[,j] - beta0[j]) + Sigma_Inv %*% sum_si_sj[,j]) 
    return(x)
  }
  
  # Function to update variance of parameter estimates.
  Sigma_Beta_function = function(n_sj){
    x = as.vector(solve(t(X)%*%X + Sigma_Inv * c(n_sj)))
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
  # and has therefore 2 neighbors)
  list_ind = list()
  for(j in 1:M){
    list_ind[[j]] = c(A$y[which(A$x == j, arr.ind = T)], A$x[which(A$y == j, arr.ind = T)])
  }
  
  # Time the function call.
  time.start = Sys.time()
  # Initialize counter for number of iterations of optimization.
  counter = 0
  
  # Initialise difference between parameter values to large value.
  diff = 100
  # Set parameters & hyperparameters to initial values.
  sigma_beta0 = params0$sigma_beta0
  beta0 = params0$beta0
  Beta=params0$Beta
  Sigma_Inv = params0$Sigma_Inv

  # Run optimization until convergence is reached.
  while(diff > eps){
    
    # Increase number of iterations by 1.
    counter = counter + 1
    # Keep old parameters to check for difference in parameter estimates.
    Beta_old = Beta
    
    # Calculate expected value of latents. 
    Eta = (cbind(rep(1, N), X) %*% rbind(beta0, Beta))
    expected_Z = structure(hutils::if_else(Y == 1, Eta  + (structure(dnorm(-Eta , mean = 0, sd = 1), dim = c(N, M)) / (structure(pnorm(Eta, mean = 0, sd = 1), dim = c(N, M)))) , Eta - (structure(dnorm(-Eta, mean = 0, sd = 1), dim = c(N, M)) / (1 - structure(pnorm(Eta, mean = 0, sd = 1), dim = c(N, M)))) ), dim = c(N, M))
    
    # Calculate sum of neighboring parameter values. 
    sum_si_sj = apply(matrix(1:M, nrow = M), 1, sum_si_sj_function)
    
    # Update spatially-varying intercept.
    beta0 = (1/(N + 1/sigma_beta0^2))*(colSums(expected_Z-X%*%Beta))
    
    # Update variance of Beta.
    Sigma_Beta  = apply(matrix(n_sj, M, 1), 1, FUN = Sigma_Beta_function)
    
    # Update spatially-varying coefficients. 
    Beta=apply(matrix(1:M, nrow = M), 1, Beta_function)
    
    # Update smoothing covariance matrix. 
    term = matrix((Beta[,A$x] - Beta[,A$y]), P)%*%t(matrix((Beta[,A$x] - Beta[,A$y]), P)) 
    Sigma_Inv = (M - 1) * solve(diag(P) + term)
  
    # Calculate difference between new and old parameters.
    diff = max(abs(Beta - Beta_old))
  }
  
  # Calculate time needed to execute function.
  time.end = Sys.time()
  time = time.end - time.start
  
  # Save parameters in list.
  params = list()
  params$beta0 = beta0
  params$Beta = Beta
  params$Sigma_Inv = Sigma_Inv
  params$counter = counter
  params$time = time 
  params$Sigma_Beta = Sigma_Beta
  return(params)
}

########################
### Simulation Study ###
########################

X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

Y = data.matrix(read.csv(paste0(path_data, "Y", sim, ".csv"), header = T)[,2:(M+1)])

# Initialize parameters 
if(init == 'Firth'){
  params0 = list()
  x = as.numeric(unlist(read.csv(sprintf("%sbeta0%03d.csv",path_Firth, sim))[,2]))
  x[is.na(x)] = 0
  params0$beta0 = x
  x = rbind2(as.numeric(unlist(read.csv(sprintf("%sBeta1%03d.csv",path_Firth, sim))[,2])),
             as.numeric(unlist(read.csv(sprintf("%sBeta2%03d.csv",path_Firth, sim))[,2]))
  )
  x[is.na(x)] = 0
  params0$Beta = x
  params0$Sigma_Inv = 1*diag(P)
  params0$sigma_beta0 = 10
}

if(init == 'random'){
params0 = list()
params0$beta0 = rnorm(M, 0, 1)
params0$Beta = matrix(rnorm(M*P, 0, 1), P, M)
params0$Sigma_Inv = 1*diag(P)
params0$sigma_beta0 = 10

}

# Parameter estimation
params = estimate_BSGLMM(X, Y, params0, eps)

# Get test statistics t = beta / sd_beta
t = params$Beta / sqrt(params$Sigma_Beta[c(1,4),])

# Save output.
write.csv(c(params$beta0), sprintf("%sbeta0%03d.csv", path, sim))
write.csv(params$Beta, sprintf("%sBeta%03d.csv", path, sim))
write.csv(params$Sigma_Beta, sprintf("%sSigma_Beta%03d.csv", path, sim))
write.csv(t, sprintf("%stest_statistic%03d.csv", path, sim))
write.csv(params$Sigma_Inv, sprintf("%sSigma_Inv%03d.csv", path, sim))
write.csv(params$Q, sprintf("%sQ%03d.csv", path, sim))
write.csv(params$counter, sprintf("%scounter%03d.csv", path, sim))
write.csv(params$time, sprintf("%stime%03d.csv", path, sim))

