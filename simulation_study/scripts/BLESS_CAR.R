##############
### Script ###
##############

# Run BLESS (= Bayesian Lesion Estimation for a Structured Spike-and-Slab Prior) with shared CAR prior over sparsity parameters 
# over all covariates (Code for Simulation Study).

#################
### Libraries ###
#################

library(hutils)
library(brglm2)
library(logOfGamma)

#################
### Constants ###
#################

# Sample size (N = {500,1000,5000}).
N = 1000
# Base rate intensity (lambda = {1,2,3}).
lambda = 3
# Number of dataset to perform analysis on.
sim = 1
# Number of voxels.
M = 2500
# Dimension of 2D slice.
dimension =  dim1 = dim2 = sqrt(M)
# Number of covariates.
P = 2
# Initialization of variational infrence algorithm (options: 'Firth_no_saved_results' (no stored files from 
# a previous parameter estimation with Firth regression), 'Firth_saved_results' (saved files stored - need
# to specify file path to stored results), 'random').
init = 'Firth_no_saved_results'
# Convergence threshold of difference in ELBO.
eps_DPE = 0.1

# Range of spike variance to evaluate in Dynamic Posterior Exploration (DPE).
v0 = exp(seq(-20, -1, length.out = 15))
n_sim_v0 = length(v0)
# Fix slab variance to large value.
v1 = 10
# Set hyperparameter standard deviation of spatially-varying intercept to large value.
sigma_beta0 = 10

# Path: output folder
path = ""
# Path: data folder
path_data = ""
# Path: folder with inital values (if Firth initialization is chosen)
path_Firth = ""

#################
### Functions ###
#################

# Logit function
logit = function(x){
  x = log(x/(1-x))
  return(x)
}

# Function to estimate parameters in BLESS.
estimate_BLESS = function(X, Y, params0, eps){

  # Inputs:
  # - X: input data (scalar covariates)
  # - Y: output data (image -> lesion masks)
  # - params0: list containing intial parameter values
  # - eps: optimization threshold = difference between ELBOs which determines when to stop optimization
  
  # Time function call.
  time.start = Sys.time()
  
  # Hyperparameters
  # spike & slab variance parameter
  v0 = params0$v0
  v1 = params0$v1
  # hyperparameter for CAR prior on precision parameter
  v = params0$v

  # Functions & other quantities
  # Exp-normalize trick to avoid numerical issues with large exponential values.
  exp_normalize = function(x){
    b = max(x)
    y = exp(x-b)
    return(y / sum(y))
  }
  
  # Logistic function
  logistic = function(x){
    x = exp(x) / (1 + exp(x))
    return(x)
  }
  
  # Function to calculate lambda(xi) within the Taylor approximation for acquiring updates of 
  # the sparsity parameter theta. 
  lambda_xi_func = function(x){
    x = (1/(2*x))*(logistic(x) - (1/2))
    return(x)
  }
  
  # Function for deriving indices of adjacent neighbors of every single voxel location for 
  # spatial CAR prior (a voxel location is considered adjacent if they share the same face).
  adjacency_matrix = function(dim1, dim2, dim3){
    
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
      colnames(A) = c('x','y')
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
      conv = as.vector(aperm(array(1:(dim1*dim2*dim3), dim = c(dim2,dim1,dim3)), perm = c(2, 1, 3)))
      
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
      colnames(A_2D) = c('x','y')
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
  
  # Function to update variance of Beta.
  var_beta_function = function(expected_gamma2){
    x = solve(t(X)%*%X + diag(expected_gamma2, P))
    return(x)
  }
  
  # Function to update spatially-varying parameters.
  Beta_function = function(j){
    x = matrix(var_Beta[,j], P, P) %*% t(X)%*%(expected_Z[,j] - beta0[j])
    return(x)
  }
  
  # Function to update gamma values (probability of inclusion). 
  gamma_function = function(j){
  
    # Exp-normalize trick
    active = (-0.5*log(v1) - (matrix(diag(matrix(var_Beta[,j], P, P) + Beta[,j]%*%t(Beta[,j])), nrow = P, ncol = 1))/(2*v1) + theta[j])
    not_active = (-0.5*log(v0) - (matrix(diag(matrix(var_Beta[,j], P, P) + Beta[,j]%*%t(Beta[,j])), nrow = P, ncol = 1))/(2*v0))
    x = apply(matrix(1:P, nrow=P), 1, function(p) exp_normalize(c(active[p], not_active[p]))[1])
    
    return(x)
  }
  
  # Function to sum neighboring values of 2D/3D lattice together.
  sum_si_sj_function = function(j){
    x = sum(theta[list_ind[[j]]])
    return(x)
  }
  
  # Initialize counter for number of iterations of optimization.
  counter = 0
  # Initialize difference in ELBO values to large value.
  diff = 100
  # Set parameters & hyperparameters to initial values.
  sigma_beta0 = 10
  beta0 = params0$beta0
  Beta = params0$Beta
  Sigma_Inv = params0$Sigma_Inv
  theta = params0$theta
  xi = params0$xi
  expected_gamma = params0$expected_gamma
  expected_gamma2 = (expected_gamma/v1) + ((1-expected_gamma)/v0)
  ELBO = NA
  Q = -300000000000000000000000000000000

  # Run optimization until convergence criteria is reached.
  while(diff > eps){
    
    # Increase number of iterations by 1.
    counter = counter + 1
    # Save old quantities. 
    Q_old = Q
 
    # Update Z.
    expected_Z = structure(hutils::if_else(Y == 1, 
                                           (cbind(rep(1, N), X) %*% rbind(beta0, Beta))  + 
                                             (structure(dnorm(-(cbind(rep(1, N), X) %*% rbind(beta0, Beta)) , mean = 0, sd = 1), dim = c(N, M)) / 
                                                (structure(pnorm((cbind(rep(1, N), X) %*% rbind(beta0, Beta)), mean = 0, sd = 1), dim = c(N, M)))) , 
                                           (cbind(rep(1, N), X) %*% rbind(beta0, Beta)) - (structure(dnorm(-(cbind(rep(1, N), X) %*% rbind(beta0, Beta)), mean = 0, sd = 1), dim = c(N, M)) / 
                                                                                            (1 - structure(pnorm((cbind(rep(1, N), X) %*% rbind(beta0, Beta)), mean = 0, sd = 1), dim = c(N, M)))) ), 
                           dim = c(N, M))
    
    # Update beta0.
    beta0 = (1/(N + 1/sigma_beta0^2))*(colSums(expected_Z) - colSums(X%*%Beta))
    
    # Update variance of Beta.
    var_Beta  = apply(expected_gamma2, 2, FUN = var_beta_function)
    
    # Update Beta.
    Beta=apply(matrix(1:M, nrow = M), 1, Beta_function)
    
    # Update gamma.
    expected_gamma = apply(matrix(1:M, nrow = M), 1, gamma_function)
    # Update expectation: E[gamma*(1/v1) + (1-gamma)*(1/v0)]
    expected_gamma2 = (expected_gamma/v1) + ((1-expected_gamma)/v0)
    
    # Update theta
    sum_si_sj = as.vector(apply(matrix(1:M, nrow = M), 1, sum_si_sj_function))
    Sigma_theta = as.vector(apply(matrix(1:M, nrow = M), 1, function(j) as.vector(solve(n_sj[j]*Sigma_Inv + 2*P*lambda_xi_func(xi[j])))))
    theta = as.vector(apply(matrix(1:M, nrow = M), 1, function(j) Sigma_theta[j]%*%(Sigma_Inv%*%sum_si_sj[j] + sum(expected_gamma[,j]) - P/2)))
    
    # Update Sigma Inverse.
    term = sum((theta[A$x] - theta[A$y])^2) + sum(Sigma_theta[A$y])
    Sigma_Inv = (M + v)/(1 + term)
    
    # Update xi.
    xi = sqrt(Sigma_theta + theta^2)
    
    # Calculate ELBO. 
    log_p_y_z = 0
    log_p_z_beta_beta0 = (-0.5)*sum(apply(matrix(1:M, nrow = M), 1, function(j) sum(diag(t(X)%*%X%*%(matrix(var_Beta[,j], P, P) + Beta[,j]%*%t(Beta[,j])))))) - (N/2)*sum((1/(N + 1/sigma_beta0^2)) + beta0^2) -
      sum(apply(matrix(1:M, nrow = M), 1, function(j) t(Beta[,j])%*%t(X)%*%matrix(rep(beta0[j], N), N, 1)))
    log_p_beta_gamma = (-0.5)*sum(expected_gamma*log(v1) + (1-expected_gamma)*log(v0)) - 
      0.5*sum(apply(matrix(1:M, nrow = M), 1, function(j) sum(diag(diag(expected_gamma2[,j], P)%*%(matrix(var_Beta[,j], P, P) + Beta[,j]%*%t(Beta[,j]))))))
    log_p_beta0 = -M*log(sigma_beta0^2) - sum((1/(2*sigma_beta0^2))*((1/(N + 1/sigma_beta0^2)) + beta0^2))
    log_p_gamma_theta = sum(sweep(expected_gamma, 2, theta, FUN = `*`)) + P*sum(log(logistic(xi) + 10^(-10))) - 0.5*P*sum(theta + xi) - 
       P*sum(lambda_xi_func(xi)*(Sigma_theta + theta^2 - xi^2))
    log_p_theta_Sigma_Inv = (-0.5)*Sigma_Inv*term
    log_p_Sigma_Inv = (-0.5)*Sigma_Inv - (1/2)*log(2) - gammaln(1/2) 
    Eta = cbind(rep(1,N), X) %*% rbind(beta0, Beta) 
    log_q_z = (-0.5)*sum(apply(matrix(1:M, nrow = M), 1, function(j) t(Eta[,j])%*%Eta[,j])) - 
      sum(Y*log(structure(pnorm(Eta, mean = 0, sd = 1), dim = dim(Eta)) + 10^(-10)) + 
            (1-Y)*log(1-structure(pnorm(Eta, mean = 0, sd = 1),dim = dim(Eta)) + 10^(-10)))
    log_q_beta = (-0.5)*sum(apply(matrix(1:M, nrow = M), 1, function(j) log(det(matrix(var_Beta[,j], P, P)) + 10^(-10)))) - (P*M)/2
    log_q_beta0 = - (M/2)*log((1/(N + 1/sigma_beta0^2))) -(M/2) 
    log_q_gamma = sum(expected_gamma*log(expected_gamma + 10^(-10)) + (1-expected_gamma)*log(1-expected_gamma + 10^(-10)))
    log_q_theta = (-0.5)*sum(log(Sigma_theta + 10^(-10))) - (M)/2 
    log_q_Sigma_Inv = -(M + 1)/2 - (((M+1))/2)*log(2) + ((M+1)/2)*log((1 + term)) - gammaln((M+1)/2)
    
    Q = log_p_y_z + log_p_z_beta_beta0 + log_p_beta_gamma + log_p_beta0 + log_p_gamma_theta + log_p_theta_Sigma_Inv + log_p_Sigma_Inv -
      log_q_z - log_q_beta - log_q_beta0 - log_q_gamma - log_q_theta - log_q_Sigma_Inv
    
    diff = (Q - Q_old)
    Eta = 0
    ELBO= append(ELBO, Q)
    
  }

  # Evaluate time needed to execute function call.
  time.end = Sys.time()
  time = time.end - time.start
  
  # Save final parameters. 
  params = list()
  params$beta0 = beta0
  params$Beta = Beta
  params$Sigma_Inv = Sigma_Inv
  params$theta = theta
  params$expected_gamma = expected_gamma
  params$expected_gamma2 = expected_gamma2
  params$xi = xi
  params$counter = counter
  params$Q = Q
  params$time = time
  params$ELBO = ELBO
  params$term = term
  return(params)
}

# Function to calculate marginal posterior of gamma. 
calculate_marginal = function(result, params0, X, y){
  
  # Set constants & hyperparameters.
  v0 = params0$v0
  v1 = params0$v1
  v = params0$v
  sigma_beta0 = params0$sigma_beta0
  N = dim(X)[1]
  P = dim(X)[2]
  M = dim(params$Beta)[2]
  dimension = dim1 = dim2 = sqrt(M) 
  
  # Functions
  # Logistic funtion.
  logistic = function(x){
    x = exp(x)/(1+exp(x))
    return(x)
  }
  
  # Function to calculate lambda(xi) within the Taylor approximation for acquiring updates of 
  # the sparsity parameter theta. 
  lambda_xi_func = function(x){
    x = (1/(2*x))*(logistic(x) - (1/2))
    return(x)
  }
  
  # Function for deriving indices of adjacent neighbors of every single voxel location for 
  # spatial CAR prior (a voxel location is considered adjacent if they share the same face).
  adjacency_matrix = function(dim1, dim2, dim3){
    
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
    
    if(missing(dim3)){
      if (dim1<3 | dim2<3){ 
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
  
  # Function to calculate sum of adjacent parameter values that are neighbors as 
  # determined by list of neighboring values.
  sum_si_sj_function = function(j){
    x = sum(theta[list_ind[[j]]])
    return(x)
  }
  
  # Set values to posterior quantities evaluated with BLESS for evaluated spike variance v0.
  beta = structure(params$Beta, dim = c(P, M))
  beta0 = params$beta0
  theta = as.vector(params$theta)
  xi = as.vector(params$xi)
  Sigma_Inv = as.vector(params$Sigma_Inv)
  exp_gamma = structure(params$expected_gamma, dim = c(P, M))
  # Acquire binary values after thresholding posterior gamma values.
  gamma = exp_gamma
  gamma[exp_gamma > 0.5] = 1
  gamma[exp_gamma <= 0.5] = 0
  
  # Sum neighboring sparsity parameters theta together ( = sum_si_sj E[theta(s_i)]).
  sum_si_sj = as.vector(apply(matrix(1:M, nrow = M), 1, sum_si_sj_function))
  
  # Acquire E[sum_si_sj [[theta(s_i) - theta(s_j)][theta(s_i) - theta(s_j)]^T].
  term = sum((theta[A$x] - theta[A$y])^2)
  
  # Evaluate marginal posterior of gamma under prior v0 = 0.
  marginal = 0
  for(j in 1:M){
    # X_gamma 
    X_gamma = X*matrix(rep(gamma[,j], N), nrow = N, ncol = P, byrow = T)
    X_gamma = X_gamma[,apply(X_gamma, 2, function(x) !all(x == 0))]
    Q = sum(gamma[,j])
    
    # Q=0 (all gamma = 0 for voxel s_j) vs Q>0 (at least one covariate has gamma = 1 f or voxel s_j)
    if(Q==0){
      eta = rep(beta0[j], N)
      Phi = structure(pnorm(-eta, mean = 0, sd = 1),dim = dim(eta))
      marginal = marginal + sum(y[,j]*log(1-Phi) + (1-y[,j])*log(Phi)) +
        0.5*t(eta)%*%eta
    } else{
      X_gamma = matrix(X_gamma, N, Q)
      beta_gamma = beta[,j]*gamma[,j]
      beta_gamma = beta_gamma[beta_gamma != 0]
      expected_gamma = exp_gamma[,j]*gamma[,j]
      expected_gamma = expected_gamma[expected_gamma != 0]
      var_beta = solve(t(X_gamma)%*%X_gamma + diag(c(expected_gamma*(1/v1)), Q))
      eta = X_gamma%*%beta_gamma + beta0[j]
      Phi = structure(pnorm(-eta, mean = 0, sd = 1), dim = dim(eta))
      
      marginal = marginal - (0.5)*sum(diag(t(X_gamma)%*%X_gamma)%*%(var_beta + beta_gamma%*%t(beta_gamma))) - 
        sum(beta0[j]*X_gamma%*%beta_gamma) - (Q/2)*log(v1) -
        (1/2)*sum(diag(rep((1/v1), Q), Q)%*%(var_beta + beta_gamma%*%t(beta_gamma))) +
        0.5*t(eta)%*%eta + sum(y[,j]*log(1-Phi) + (1-y[,j])*log(Phi )) + 
        0.5*log(det(var_beta)) + (Q/2) 
    }
  }
  
  marginal = marginal - (N/2)*sum(((1/(N + 1/sigma_beta0^2)) + beta0^2))  - 
    sum((1/(2*sigma_beta0^2))*((1/(N + 1/sigma_beta0^2)) + beta0^2)) +
    sum(sweep(gamma, 2, log(logistic(theta)), FUN = `*`) + sweep((1-gamma), 2, log(1 - logistic(theta) ), FUN = `*`)) - ((M)/2)*log(2*pi) +
    (M/2)*log(Sigma_Inv) - 0.5*Sigma_Inv*term + 
    ((-1)/2)*log(Sigma_Inv) - ((1)/2)*log(2)  - gammaln((1)/2) - 
    0.5*Sigma_Inv + (M/2) + 0.5*sum(log((1/(N + 1/sigma_beta0^2)) + beta0^2)) - 
    (M/2)*log(sigma_beta0^2)
  return(marginal)
}

########################
### Simulation Study ###
########################

ELBO_data = list()

# Data
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)
Y = data.matrix(read.csv(paste0(path_data, "Y", sim, ".csv"), header = T)[,2:(M+1)])

# Initialization
if(init == 'random'){
  params0 = list()
  params0$beta0 = rnorm(M, 0, 1)
  params0$Beta = matrix(rnorm(P*M, 0, 1), nrow = P, ncol = M)
  params0$Sigma_Inv = 1
  params0$theta = rnorm(M, 0, 1)
  params0$v0 = max(v0)
  params0$v1 = v1 
  params0$v = 1
  params0$xi = rnorm(M, 0, 1)
  params0$sigma_beta0 = sigma_beta0
  params0$expected_gamma = matrix(0.5, P, M)
  params0$Q = - 30000000000000000000000000
}

if(init == 'Firth_saved_results'){
  params0 = list()
  x = as.numeric(unlist(read.csv(sprintf("%sbeta0%03d.csv", path_Firth, sim))[,2]))
  x[is.na(x)] = 0
  params0$beta0 = x
  x = rbind2(as.numeric(unlist(read.csv(sprintf("%sBeta1%03d.csv", path_Firth, sim))[,2])),
             as.numeric(unlist(read.csv(sprintf("%sBeta2%03d.csv", path_Firth, sim))[,2]))
  )
  x[is.na(x)] = 0
  params0$Beta = x
  params0$Sigma_Inv = 1
  x = rbind2(as.numeric(unlist(read.csv(sprintf("%sstd_error_Beta1%03d.csv", path_Firth, sim))[,2])),
             as.numeric(unlist(read.csv(sprintf("%sstd_error_Beta2%03d.csv", path_Firth, sim))[,2]))
  )
  t = params0$Beta / x
  t[abs(t) > 1.96] = 1
  t[is.na(t)] = 0
  t[t != 1] = 0
  params0$theta = logit(rep(mean(t), M))
  params0$xi = sqrt(params0$theta^2)
  params0$v = 1
  params0$sigma_beta0 = sigma_beta0
  params0$v0 = max(v0)
  params0$v1 = v1
  params0$expected_gamma = matrix(0.5, P, M)
  params0$Q = - 30000000000000000000000000
}

if(init == 'Firth_no_saved_results'){
  params0 = list()
  # Run independent Firth regression for every voxel s_j.
  params_Firth = matrix(NA, nrow = P+1, ncol = M)
  var_Firth = matrix(NA, nrow = P+1, ncol = M)
  y_pred_Firth = matrix(NA, nrow = 4, ncol = M)
  for (j in 1:M) {
    data = data.frame(cbind(Y[,j], X))
    colnames(data) = c('y', 'x1', 'x2')
    tryCatch({
      model = glm(y ~ 1 + x1 + x2, family = binomial(probit), data = data, method = "brglmFit", type = 'AS_mean', maxit = 10000, epsilon = 1e-05, slowit = 0.5)
      y_pred_Firth[,j] = fitted(model)[c(1, (N*0.25 + 1), (N*0.5 + 1),(N*0.75 + 1))]
      model = summary(model)
      params_Firth[,j] = matrix(model$coefficients[,1], nrow = (P+1), ncol = 1)
      var_Firth[,j] = matrix(model$coefficients[,2], nrow = (P+1), ncol = 1)
    }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
  }
  params0$beta0 = params_Firth[1,]
  params0$beta0[is.na(params0$beta0)] = 0
  params0$Beta = params_Firth[2:3,]
  params0$Beta[is.na(params0$Beta)] = 0
  params0$Sigma_Inv = 1
  t = params0$Beta / var_Firth[2:3,]
  t[abs(t) > 1.96] = 1
  t[is.na(t)] = 0
  t[t != 1] = 0
  params0$theta = logit(rep(mean(t), M))
  params0$xi = sqrt(params0$theta^2)
  params0$v = 1
  params0$sigma_beta0 = sigma_beta0
  params0$v0 = max(v0)
  params0$v1 = v1
  params0$expected_gamma = matrix(0.5, P, M)
  params0$Q = - 30000000000000000000000000
}

# Matrix setup for results
gamma1_list = matrix(NA, nrow = n_sim_v0, ncol = M)
gamma2_list = matrix(NA, nrow = n_sim_v0, ncol = M)
gamma21_list = matrix(NA, nrow = n_sim_v0, ncol = M)
gamma22_list = matrix(NA, nrow = n_sim_v0, ncol = M)
Beta1_list = matrix(NA, nrow = n_sim_v0, ncol = M)
Beta2_list = matrix(NA, nrow = n_sim_v0, ncol = M)
time_list = matrix(NA, nrow = n_sim_v0, ncol = 1)
theta1_list = matrix(NA, nrow = n_sim_v0, ncol = M)
counter_list = matrix(NA, nrow = n_sim_v0, ncol = 1)
Sigma_Inv_list = matrix(NA, nrow = n_sim_v0, ncol = 1)
xi1_list = matrix(NA, nrow = n_sim_v0, ncol = M)
beta0_list = matrix(NA, nrow = n_sim_v0, ncol = M)
ELBO_list = matrix(NA, nrow = n_sim_v0, ncol = 1)
marginal_list = numeric(n_sim_v0)

# Backwards DPE: hence start from the largest v0 value and go backwards in the spike variance range
# taking the rest of the estimated parameter values as a warm-start. 
for(sim_v0 in n_sim_v0:1){
  
  # Print stage of DPE.
  print(sim_v0)

  # If not at start of DPE, set inital values params0 to  previously estimated parameter estimates 
  # and set spike variance to a lower value in the spike variance range. 
  if(sim_v0 != n_sim_v0){
    params0 = list()
    params0$beta0 = beta0_list[(sim_v0+1),]
    params0$beta0 = c(params0$beta0)
    params0$Beta = matrix(NA, nrow = P, ncol = M)
    params0$Beta[1,] = Beta1_list[(sim_v0+1),]
    params0$Beta[2,] = Beta2_list[(sim_v0+1),]
    params0$expected_gamma = matrix(NA, nrow = P, ncol = M)
    params0$expected_gamma[1,] = gamma1_list[(sim_v0+1),]
    params0$expected_gamma[2,] = gamma2_list[(sim_v0+1),]
    params0$Sigma_Inv = c(Sigma_Inv_list[(sim_v0+1),])
    params0$theta = as.vector(theta1_list[(sim_v0+1),])
    params0$xi = as.vector(xi1_list[(sim_v0+1),])
    params0$v = 1
    params0$sigma_beta0 = sigma_beta0
    
    params0$v0 = v0[sim_v0]
    params0$v1 = v1
  }

 	# Parameter estimation
  params = estimate_BLESS(X, Y, params0, eps_DPE)

  # Storage of simulated parameters
  ELBO_data[[sim_v0]] = params$ELBO
  gamma1_list[sim_v0, ] = as.vector(params$expected_gamma[1,])
  gamma2_list[sim_v0, ] = as.vector(params$expected_gamma[2,])
  gamma21_list[sim_v0, ] = as.vector(params$expected_gamma2[1,])
  gamma22_list[sim_v0, ] = as.vector(params$expected_gamma2[2,])
  Beta1_list[sim_v0, ] = as.vector(params$Beta[1,])
  Beta2_list[sim_v0, ] = as.vector(params$Beta[2,])
  time_list[sim_v0, ] = as.vector(params$time)
  counter_list[sim_v0, ] = as.vector(params$counter)
  Sigma_Inv_list[sim_v0, ] = as.vector(params$Sigma_Inv)
  theta1_list[sim_v0, ] = as.vector(params$theta)
  xi1_list[sim_v0, ] = as.vector(params$xi)
  beta0_list[sim_v0, ] = as.vector(params$beta0)
  ELBO_list[sim_v0,] = as.vector(params$Q)
  marginal_list[sim_v0] = c(calculate_marginal(params, params0, X, Y))
  
  # Find lowest marginal posterior of gamma value and hence optimal spike variance over the range 
  # of evaluated spike variances (Optimal value will be the last evaluated quantitiy closest to 
  # v0 = 0 as we evaluate the marginal posterior of gamma under v0 = 0!)
  # However, good sanity check to plot marginal values and check for plateu of values towards 0
  # which implies a stabilization of the parameter estimates! 
  opt = marginal_list
  opt[marginal_list == 0] = NA
  optimal_v0_idx = which.max(opt)
  optimal_v0 = v0[optimal_v0_idx]
  
  # Save parameter estimates.
  write.csv(marginal_list, sprintf("%smarginal%03d.csv", path, sim))
  write.csv(optimal_v0, sprintf("%soptimal_v0%03d.csv", path, sim))
  write.csv(optimal_v0_idx, sprintf("%soptimal_v0_idx%03d.csv", path, sim))
  write.csv(beta0_list, sprintf("%sbeta0%03d.csv", path, sim))
  write.csv(Beta1_list, sprintf("%sBeta1%03d.csv", path, sim))
  write.csv(Beta2_list, sprintf("%sBeta2%03d.csv", path, sim))
  write.csv(gamma1_list, sprintf("%sgamma1%03d.csv", path, sim))
  write.csv(gamma2_list, sprintf("%sgamma2%03d.csv", path, sim))
  write.csv(gamma21_list, sprintf("%sgamma21%03d.csv", path, sim))
  write.csv(gamma22_list, sprintf("%sgamma22%03d.csv", path, sim))
  write.csv(xi1_list, sprintf("%sxi1%03d.csv", path, sim))
  write.csv(theta1_list, sprintf("%stheta1%03d.csv", path, sim))
  write.csv(Sigma_Inv_list, sprintf("%sSigma_Inv%03d.csv", path, sim))
  write.csv(counter_list, sprintf("%scounter%03d.csv", path, sim))
  write.csv(time_list, sprintf("%stime%03d.csv", path, sim))
  write.csv(ELBO_list, sprintf("%sELBO%03d.csv", path, sim))
  save(ELBO_data, file = sprintf("%sELBO_data%03d.RData", path, sim))  
}

write.csv(v0, sprintf("%sv0.csv", path))
