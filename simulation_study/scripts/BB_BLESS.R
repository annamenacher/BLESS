##############
### Script ###
##############

# Run BB-BLESS (= Bayesian Bootstrap - Bayesian Lesion Estimation for a Structured Spike-and-Slab Prior) with MCAR prior over sparsity parameters.
# (Code for Simulation Study)

#################
### Libraries ###
#################

library(LaplacesDemon)
library(brglm2)
library(logOfGamma)

#################
### Constants ###
#################

# path: data folder
path_data = ''
# path: folder with output from parameter estimation with BLESS (backwards DPE results)
path_BLESS = ''
# path: output folder
path = ''

# Sample size (N = {500,1000,5000}).
N = 1000
# Base rate intensity (lambda = {1,2,3}).
lambda = 3
# Number of dataset to perform analysis on.
dataset = 1
# Number of voxels.
M = 2500
# Dimension of 2D slice.
dimension =  dim1 = dim2 = sqrt(M)
# Number of covariates.
P = 2
# Initialization of variational infrence algorithm (options: Firth, random).
init = 'DPE'
# Convergence threshold of difference in parameter estimates negligible small.
eps = 0.001
# alpha: concentration parameter of Dirichlet distribution
dir_alpha = 1
# Set hyperparameter standard deviation of spatially-varying intercept to large value.
sigma_beta0 = 10
# Number of bootstrap replicates
B = 1500
# Hyperparameter: last value of backwards DPE procedure
v0_idx = 1
# Range of spike variances over which DPE was run (in simulation studies: range of 15 spike variances over exp(-20) to exp(-1)). 
#v0 = as.numeric(unlist(read.csv(paste0(path_BLESS, 'v0.csv'))[,2]))
v0 = exp(seq(-20, -1, length.out = 15))
v0 = v0[v0_idx]
# Hyperparameter: slab variance
v1 = 10

# Data
X = matrix(cbind(c(rep(1, N/2),rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)
Y = data.matrix(read.csv(paste0(path_data, "Y", dataset, ".csv"), header = T)[,2:(M+1)])

#################
### Functions ###
#################

# Logit function
logit = function(x){
  x = log(x/(1-x))
  return(x)
}

# Function to estimate BB-BLESS model.
bb_bless = function(X, Y, params0, eps){

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
  # degrees of freedom of Wishart prior on covariance matrix of MCAR prior
  v = P
  # weights drawn from Dirichlet prior distribution in order to re-weight the likelihood
  w = N * c(rdirichlet(1, alpha))
  # perturbation of prior mean (='jitter') drawn from spike distribution
  mu = rnorm(M*P, 0, sqrt(v0))
  mu = matrix(mu, P, M)
  
  # Functions & other quantities
  # Logistic function
  logistic = function(x){
    #x = exp(x) / (1 + exp(x))
    #return(x)
    x = ifelse(x >= 0, 1 / ( 1 + exp(-x) ), exp(x) / ( 1 + exp(x) ))
    return(x)
  }
  
  # Exp-normalize trick to avoid numerical issues with large exponential values.
  exp_normalize = function(x){
    b = max(x)
    y = exp(x-b)
    return(y / sum(y))
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
    
    # Inputs:
    # - dim1: dimension1
    # - dim2: dimension2
    # - dim3: dimension3
   
    if(missing(dim3)){
      A = data.frame(x = integer(),y = integer())
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

    # Inputs:
    # - dim1: dimension1
    # - dim2: dimension2
    # - dim3: dimension3
    
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

    # Inputs:
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
  
  # Function to update variance of Beta.
  var_beta_function = function(expected_gamma2){
    x = solve(XWX + diag(expected_gamma2, P))
    return(x)
  }
  
  # Function to update spatially-varying parameters.
  Beta_function = function(j){
    x = matrix(var_Beta[,j], P, P) %*%(XW%*%(expected_Z[,j] - beta0[j])+ diag(expected_gamma2[,j], P)%*%matrix(mu[,j], P, 1))
    return(x)
  }
  
  # Function to update gamma values (probability of inclusion). 
  gamma_function = function(j){
    active = (-0.5*log(v1) - (matrix(diag(matrix(var_Beta[,j], P, P)), nrow = P, ncol = 1) + Beta[,j]^2 - 2*Beta[,j]*mu[,j] + mu[,j]^2)/(2*v1) + theta[,j])
    not_active = (-0.5*log(v0) - (matrix(diag(matrix(var_Beta[,j], P, P)), nrow = P, ncol = 1) + Beta[,j]^2 - 2*Beta[,j]*mu[,j] + mu[,j]^2)/(2*v0))
    x = apply(matrix(1:P, nrow = P), 1, function(p) exp_normalize(c(active[p],not_active[p]))[1])
    return(x)
  }
  
  # Function to sum neighboring values of 2D/3D lattice together.
  sum_si_sj_function = function(j){
    x = matrix(apply(matrix(theta[,list_ind[[j]]], nrow = P), 1, sum), nrow = P)
    return(x)
  }
  
  # Function to update variational parameter xi.
  xi_function = function(j){
    x = sqrt((diag(matrix(Sigma_theta[,j], P, P) + theta[,j]%*%t(theta[,j]))))
    return(x)
  }
  
  # Initialize counter for number of iterations of optimization.
  counter = 0
  # Initialize difference in ELBO values to large value.
  diff = 100
  # Set parameters & hyperparameters to initial values.
  sigma_beta0 = params0$sigma_beta0
  beta0 = params0$beta0
  Beta = params0$Beta
  Sigma_Inv = params0$Sigma_Inv
  theta = params0$theta
  xi = params0$xi
  expected_gamma = params0$expected_gamma
  expected_gamma2 = (expected_gamma/v1) + ((1-expected_gamma)/v0)
  Q = -3000000000000000000000
  ELBO = NA
  # Calculate quantities repeatedly used in below updates.
  XW = t(sweep(X, 1, w, FUN = `*`))
  XWX = t(sweep(X, 1, sqrt(w), FUN = `*`))%*%sweep(X, 1, sqrt(w), FUN = `*`)
  WX = sweep(X, 1, w, FUN = `*`)
  WY = sweep(Y, 1, w, FUN = `*`)
  W1minusY = sweep((1-Y), 1, w, FUN = `*`)

  # Run optimization until convergence criteria is reached.
  while(diff > eps){
    
    # Increase number of iterations by 1.
    counter = counter + 1
    # Save old quantities. 
    Q_old = Q
    counter = counter + 1
    Q_old = Q
    Beta_old = Beta
    
    # Update Z.
    expected_Z = structure(hutils::if_else(Y == 1,
                                           (cbind(rep(1, N), X) %*% rbind(beta0, Beta))  +
                                             (structure(dnorm(-(cbind(rep(1, N), X) %*% rbind(beta0, Beta)) , mean = 0, sd = 1), dim = c(N, M)) /
                                                (structure(pnorm((cbind(rep(1, N), X) %*% rbind(beta0, Beta)), mean = 0, sd = 1) + 10^(-10), dim = c(N, M)))) ,
                                           (cbind(rep(1, N), X) %*% rbind(beta0, Beta)) - (structure(dnorm(-(cbind(rep(1, N), X) %*% rbind(beta0, Beta)), mean = 0, sd = 1), dim = c(N, M)) /
                                                                                            (1 - structure(pnorm((cbind(rep(1, N), X) %*% rbind(beta0, Beta)), mean = 0, sd = 1), dim = c(N, M)) + 10^(-10))) ),
                           dim=c(N, M))

    # Update beta0.
    beta0 = (1/(sum(w) + 1/sigma_beta0^2))*(as.vector(matrix(w,1,N)%*%(expected_Z - X%*%Beta)))
    
    # Update variance of Beta.
    var_Beta  = apply(expected_gamma2, 2, FUN = var_beta_function)
    
    # Update Beta.
    Beta=apply(matrix(1:M, nrow = M), 1, Beta_function)

    # Update gamma.
    expected_gamma = apply(matrix(1:M, nrow = M), 1, gamma_function)

    # Update expectation: E[gamma*(1/v1) + (1-gamma)*(1/v0)]
    expected_gamma2 = (expected_gamma/v1) + ((1-expected_gamma)/v0)
    
    # Update theta
    sum_si_sj = apply(matrix(1:M, nrow = M), 1, sum_si_sj_function)
    Sigma_theta = apply(matrix(1:M, nrow = M), 1, function(j) as.vector(solve(n_sj[j]*Sigma_Inv + 2*diag(c(lambda_xi_func(xi[,j])), P))))
    theta = apply(matrix(1:M, nrow = M), 1, function(j) matrix(Sigma_theta[,j], P, P)%*%(Sigma_Inv%*%sum_si_sj[,j] + expected_gamma[,j] - 1/2))
    
    # Update Sigma Inverse.
    term_ind = matrix(1:(P*P), P, P)
    term = matrix((theta[,A$x] - theta[,A$y]), P)%*%t(matrix((theta[,A$x] - theta[,A$y]), P)) + matrix(rowSums(Sigma_theta[,A$y]), P)
    Sigma_Inv = (M + v - P - 1)*solve(diag(P) + term)

    # Update xi.
    xi = apply(matrix(1:M, nrow = M), 1, xi_function)
    
    # ELBO
    expected_Z = 0
    log_p_y_z = 0
    log_p_z_beta_beta0 = (-0.5)*sum(apply(matrix(1:M, nrow = M), 1, function(j) sum(diag(XWX%*%(matrix(var_Beta[,j], P, P) + Beta[,j]%*%t(Beta[,j])))))) - 
      (sum(w)/2)*sum((1/(sum(w) + 1/sigma_beta0^2)) + beta0^2) -
      sum(apply(matrix(1:M, nrow = M), 1, function(j) t(Beta[,j])%*%t(WX)%*%matrix(rep(beta0[j], N), N, 1)))
    log_p_beta_gamma = (-0.5)*sum(expected_gamma*log(v1) + (1-expected_gamma)*log(v0)) - 
      0.5*sum(apply(matrix(1:M, nrow = M), 1, function(j) sum(diag(diag(expected_gamma2[,j], P)%*%(matrix(var_Beta[,j], P, P) + Beta[,j]%*%t(Beta[,j])))))) -
      0.5*sum(apply(matrix(1:M, nrow = M), 1, function(j) sum(diag(diag(expected_gamma2[,j], P)%*%matrix(mu[,j], P, 1)%*%t(matrix(mu[,j], P, 1)))))) +
      sum(apply(matrix(1:M, nrow = M), 1, function(j) sum(diag(diag(expected_gamma2[,j], P)%*%matrix(mu[,j], P, 1)%*%t(Beta[,j])))))
    log_p_beta0 = -(M/2)*log(sigma_beta0^2) - sum((1/(2*sigma_beta0^2))*((1/(N + 1/sigma_beta0^2)) + beta0^2))
    log_p_gamma_theta = sum(theta*expected_gamma) + sum(log(logistic(xi) + 10^(-10))) - 0.5*sum(theta + xi) - 
      sum(t(apply(matrix(1:P, ncol = P), 2, function(i) sum(lambda_xi_func(xi[i,])*(Sigma_theta[term_ind[i,i],] + theta[i,]^2 - xi[i,]^2)))))
    log_p_theta_Sigma_Inv = (-0.5)*sum(diag(Sigma_Inv%*%term))
    log_p_Sigma_Inv = (-0.5)*sum(diag(diag(P)%*%Sigma_Inv)) - ((P^2)/2)*log(2) - sum(apply(matrix(1:P, nrow = P), 1, function(p) gammaln((P+1-p)/2))) - (P/2)*log(P)
    Eta = cbind(rep(1,N), X) %*% rbind(beta0, Beta) 
    log_q_z = (-0.5)*sum(apply(matrix(1:M, nrow = M), 1, function(j) t(Eta[,j]*w)%*%Eta[,j])) - 
      sum(WY*log(structure(pnorm(Eta, mean = 0, sd = 1), dim = dim(Eta)) + 10^(-10)) + W1minusY*log(1-structure(pnorm(Eta, mean = 0, sd = 1),dim = dim(Eta)) + 10^(-10)))
    log_q_beta = (-0.5)*sum(apply(matrix(1:M, nrow = M), 1, function(j) log(det(matrix(var_Beta[,j], P, P)) + 10^(-10)))) - (P*M)/2
    log_q_beta0 = - (M/2)*log((1/(N + 1/sigma_beta0^2))) -(M/2) 
    log_q_gamma = sum(expected_gamma*log(expected_gamma + 10^(-10)) + (1-expected_gamma)*log(1-expected_gamma + 10^(-10)))
    log_q_theta = (-0.5)*sum(apply(matrix(1:M, nrow = M), 1, function(j) sum(log(det(matrix(Sigma_theta[,j], P, P) + 10^(-10)))))) - (P*M)/2 
    log_q_Sigma_Inv = -(M + P)/2 - (((M+P)*P)/2)*log(2) - ((M+P)/2)*log(det(solve(diag(P) + term))) - sum(apply(matrix(1:P, nrow = P),1, function(p) gammaln((M+P+1-p)/2)))
    
    Q = log_p_y_z + log_p_z_beta_beta0 + log_p_beta_gamma + log_p_beta0 + log_p_gamma_theta + log_p_theta_Sigma_Inv + log_p_Sigma_Inv -
      log_q_z - log_q_beta - log_q_beta0 - log_q_gamma - log_q_theta - log_q_Sigma_Inv
    Eta = 0
    ELBO= append(ELBO, Q)
    
    diff = max(abs(Beta - Beta_old))

    if(counter %% 100){
      
      write.csv(c(beta0), sprintf("%sbeta0%04d.csv", path, sim))
      write.csv(Beta[1,], sprintf("%sbeta1%04d.csv", path, sim))
      write.csv(Beta[2,], sprintf("%sbeta2%04d.csv", path, sim))
      write.csv(theta[1,], sprintf("%stheta1%04d.csv", path, sim))
      write.csv(theta[2,], sprintf("%stheta2%04d.csv", path, sim))
      write.csv(Sigma_Inv, sprintf("%sSigma_Inv%04d.csv", path, sim))
      write.csv(expected_gamma[1,], sprintf("%sgamma1%04d.csv", path, sim))
      write.csv(expected_gamma[2,], sprintf("%sgamma2%04d.csv", path, sim))
      write.csv(xi[1,], sprintf("%sxi1%04d.csv", path, sim))
      write.csv(xi[2,], sprintf("%sxi2%04d.csv", path, sim))
      write.csv(c(counter), sprintf("%scounter%04d.csv", path, sim))
      write.csv(c(ELBO[2:length(ELBO)]), sprintf("%sELBO%04d.csv", path, sim))
    }
  }
  
  # Save results. 
  time.end = Sys.time()
  time = time.end - time.start
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
  
  return(params)
}

# Run BB-BLESS for B number of bootstraps. 
for(sim in 1:B){
  # Print bootstrap sample run.
  print(sim)
  
  # Initialization
  if(init == 'DPE'){
    idx = v0_idx
    params0 = list()
    params0$beta0 = as.numeric(unlist(read.csv(sprintf("%sbeta0%03d.csv", path_BLESS, dataset))[idx,2:(M+1)]))
    params0$Beta = data.matrix(rbind2(as.numeric(unlist(read.csv(sprintf("%sBeta1%03d.csv", path_BLESS, dataset))[idx,2:(M+1)])), as.numeric(unlist(read.csv(sprintf("%sBeta2%03d.csv", path_BLESS, dataset))[idx,2:(M+1)]))))
    params0$theta = data.matrix(rbind2(as.numeric(unlist(read.csv(sprintf("%stheta1%03d.csv", path_BLESS, dataset))[idx,2:(M+1)])), as.numeric(unlist(read.csv(sprintf("%stheta2%03d.csv", path_BLESS, dataset))[idx,2:(M+1)]))))
    params0$Sigma_Inv = matrix(as.numeric(unlist(read.csv(sprintf("%sSigma_Inv%03d.csv", path_BLESS, dataset))[idx,2:(P+1)])), P, P)
    params0$xi = data.matrix(rbind2(as.numeric(unlist(read.csv(sprintf("%sxi1%03d.csv", path_BLESS, dataset))[idx,2:(M+1)])), as.numeric(unlist(read.csv(sprintf("%sxi2%03d.csv", path_BLESS, dataset))[idx,2:(M+1)]))))
    params0$v = P
    params0$sigma_beta0 = sigma_beta0
    params0$expected_gamma = data.matrix(rbind2(as.numeric(unlist(read.csv(sprintf("%sgamma1%03d.csv", path_BLESS, dataset))[idx,2:(M+1)])), as.numeric(unlist(read.csv(sprintf("%sgamma2%03d.csv", path_BLESS, dataset))[idx,2:(M+1)]))))
    params0$Q = - 30000000000000000000000000
  }
  params0$v0 = v0
  params0$v1 = v1
  alpha = rep(dir_alpha, N)
  
result = bb_bless(X, Y, params0, eps)

write.csv(c(result$beta0), sprintf("%sbeta0%04d.csv", path, sim))
write.csv(result$Beta[1,], sprintf("%sbeta1%04d.csv", path, sim))
write.csv(result$Beta[2,], sprintf("%sbeta2%04d.csv", path, sim))
write.csv(result$theta[1,], sprintf("%stheta1%04d.csv", path, sim))
write.csv(result$theta[2,], sprintf("%stheta2%04d.csv", path, sim))
write.csv(result$Sigma_Inv, sprintf("%sSigma_Inv%04d.csv", path, sim))
write.csv(result$expected_gamma[1,], sprintf("%sgamma1%04d.csv", path,sim))
write.csv(result$expected_gamma[2,], sprintf("%sgamma2%04d.csv", path, sim))
write.csv(result$xi[1,], sprintf("%sxi1%04d.csv", path, sim))
write.csv(result$xi[2,], sprintf("%sxi2%04d.csv", path, sim))
write.csv(c(result$counter), sprintf("%scounter%04d.csv", path, sim))
write.csv(c(result$time), sprintf("%stime%04d.csv", path, sim))
write.csv(c(result$ELBO[2:length(result$ELBO)]), sprintf("%sELBO%04d.csv", path, sim))

# Check convergence.
write.csv(1, sprintf("%sfinal%04d.csv", path, sim))

}

