##############
### Script ###
##############

# Run BLESS (= Bayesian Lesion Estimation for a Structured Spike-and-Slab Prior) with MCAR prior over sparsity parameters.
# (3D lesion masks with parameter estimation and inference for masked version of 3D scan)
# (Code for UK Biobank application)

#################
### Libraries ###
#################

library(hutils)
library(data.table)
library(brglm2)
library(logOfGamma)

#################
### Constants ###
#################

# Sample size.
N = 2000
# Dimensions of 3D MRI scan / lesion mask.
dim1 = 91  
dim2 = 109
dim3 = 91
# Total number of voxels in 3D MRI scan.
M_image = dim1 * dim2 * dim3 
# Number of covariates (age, sex, headsize scaling factor, age by sex).
P = 4
# Convergence threshold of difference in ELBO.
eps = 10
# Initialization of variational infrence algorithm (options: 'Firth', 'random').
init = 'Firth'

# Number of spike variances to evaluate over. 
n_sim_v0 = n_sim  = 5
# Backwards DPE
DPE_list = n_sim:1
# Hyperparameter: slab variance (fixed to large value)
v1 = 10
# Hyperparameter: degrees of freedom for Wishart prior over precision matrix in MCAR prior 
v = P
# Hyperparameter: standard deviation in prior over spatially-varying intercept
sigma_beta0 = 10
# Starting point of spike variances to evaluate over: set starting point to a point where you don't see 
# any activation in model and start evaluation of spike variances under marginal posterior of gamma from there. 
# Note: a good way of doing this is to evaluate a fixed number of models with independent different spike variances
# in range of 0 < v0 < v1.
v0_start = log(0.05)
# Set range of spike variance values.
v0_list = exp(seq(-10, v0_start, length.out = n_sim))

# Set paths.
# path: folder with data
path_data = ""
# path: folder with mask & spatial adjacency list & vector with number of adjacent neighbors
path_mask = ""
# path: folder with Firth parameter estimates (for initialization)
path_Firth = ""
# path: folder for results of backwards DPE 
path_DPE = paste0("")

# Read in indices of mask, adjacency indices matrix & vector containing number of adjacent neighbors.
ind_mask = as.numeric(unlist(read.csv(paste0(path_mask, "ind_mask.csv"))[,2]))
A = data.frame(read.csv(paste0(path_mask, "A_mask.csv"))[,2:3])
colnames(A) = c('x', 'y')
n_sj = as.numeric(unlist(read.csv(paste0(path_mask, "n_sj_mask.csv"))[,2]))
load(paste0(path_mask, "list_ind_mask.RData"))

# Total number of voxels after masking lesion masks.
M = length(ind_mask)

# Data
# X (standardise data)
X = data.matrix(data.table::fread(paste0(path_data, "X_standardised_N", N, ".csv"), header = T)[,2:(P+1)])

# Y (binary lesion masks: format NxM matrix where every 3D lesion mask (.nii.gz file was vectorised and masked)
Y = data.matrix(data.table::fread(paste0(path_data, "Y_mask_N", N, ".csv"), header = T)[,2:(M+1)])

#################
### Functions ###
#################

# Logit function.
logit = function(x){
  return(log(x/(1-x)))
}

# Function to estimate parameter estimates of BLESS.
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
  # degrees of freedom of Wishart prior on covariance matrix of MCAR prior
  v = params0$v
  # Term for evaluation of Wishart covariance matrix.
  term_ind = matrix(1:(P*P), P, P)
  
  # Functions & other quantities
  # Exp-normalize trick to avoid numerical issues with large exponential values.
  exp_normalize = function(x){
    b = max(x)
    y = exp(x - b)
    return(y / sum(y))
  }
  
  # Logistic function which avoids numerical instabilities.
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
    active = (-0.5*log(v1) - (matrix(diag(matrix(var_Beta[,j], P, P) + Beta[,j]%*%t(Beta[,j])), nrow = P, ncol = 1))/(2*v1) + theta[,j])
    not_active = (-0.5*log(v0) - (matrix(diag(matrix(var_Beta[,j], P, P) + Beta[,j]%*%t(Beta[,j])), nrow = P, ncol = 1))/(2*v0))
    x = apply(matrix(1:P, nrow = P), 1, function(p) exp_normalize(c(active[p], not_active[p]))[1])
    
    return(x)
  }
  
  # Function to sum neighboring values together.
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
  ELBO = NA
  Q = -3000000000000000000000000000
  
  # Run optimization until convergence criteria is reached.
  while(diff > eps){
    
    # Avoid memory issues by clearing unnecessairly loaded quantities every iteration of the optimization.
    gc()    
    
    # Increase number of iterations by 1.
    counter = counter + 1
    # Save old quantities. 
    Q_old = Q

    # Update Z.
    expected_Z = structure(hutils::if_else(Y == 1, 
                                           (cbind(rep(1, N), X) %*% rbind(beta0, Beta))  + 
                                             (structure(dnorm(-(cbind(rep(1, N), X) %*% rbind(beta0, Beta)) , mean = 0, sd = 1), dim=c(N, M)) / 
                                                (structure(pnorm((cbind(rep(1, N), X) %*% rbind(beta0, Beta)), mean = 0, sd = 1) + 10^(-10),dim = c(N, M)))) , 
                                           (cbind(rep(1, N), X) %*% rbind(beta0, Beta)) - (structure(dnorm(-(cbind(rep(1, N), X) %*% rbind(beta0, Beta)), mean = 0, sd = 1),dim = c(N, M)) / 
                                                                                            (1 + 10^(-10) - structure(pnorm((cbind(rep(1, N), X) %*% rbind(beta0, Beta)), mean = 0, sd = 1),dim = c(N, M)))) ), 
                           dim = c(N, M))
    
    # Update beta0.
    beta0 = (1/(N + 1/sigma_beta0^2))*(colSums(expected_Z) - colSums(X%*%Beta))
    
    # Update variance of Beta.
    var_Beta  = apply(expected_gamma2, 2, FUN = var_beta_function)
    
    # Update beta.
    Beta=apply(matrix(1:M, nrow = M), 1, Beta_function)
    
    # Update gamma.
    expected_gamma = apply(matrix(1:M, nrow = M), 1, gamma_function)
    
    # Update expectation: E[gamma*(1/v1) + (1-gamma)*(1/v0)]
    expected_gamma2 = (expected_gamma/v1) + ((1-expected_gamma)/v0)
    
    # Update theta.
    sum_si_sj = apply(matrix(1:M, nrow = M), 1, sum_si_sj_function)
    Sigma_theta = apply(matrix(1:M, nrow = M), 1, function(j) as.vector(solve(n_sj[j]*Sigma_Inv + 2*diag(c(lambda_xi_func(xi[,j])), P))))
    theta = apply(matrix(1:M, nrow = M), 1, function(j) matrix(Sigma_theta[,j], P, P)%*%(Sigma_Inv%*%sum_si_sj[,j] + expected_gamma[,j] - 1/2))
    
    # Update Sigma Inverse.
    term_ind = matrix(1:(P*P), P, P)
    term = matrix((theta[,A$x] - theta[,A$y]), P)%*%t(matrix((theta[,A$x] - theta[,A$y]), P)) + matrix(rowSums(Sigma_theta[,A$y]), P)
    Sigma_Inv = (M + v - P - 1)*solve(diag(P) + term)
    
    # Update xi.
    xi = apply(matrix(1:M, nrow = M), 1, xi_function)
    
    # Calculate ELBO. 
    expected_Z = 0
    log_p_y_z = 0
    log_p_z_beta_beta0 = (-0.5) * sum(apply(matrix(1:M, nrow = M), 1, function(j) sum(diag(t(X)%*%X%*%(matrix(var_Beta[,j], P, P) + Beta[,j]%*%t(Beta[,j])))))) - (N/2)*sum((1/(N + 1/sigma_beta0^2)) + beta0^2) -
      sum(apply(matrix(1:M, nrow = M), 1, function(j) t(Beta[,j])%*%t(X)%*%matrix(rep(beta0[j], N), N, 1)))
    log_p_beta_gamma = (-0.5) * sum(expected_gamma * log(v1) + (1-expected_gamma) * log(v0)) - 
      0.5*sum(apply(matrix(1:M, nrow=M), 1, function(j) sum(diag(diag(expected_gamma2[,j], P)%*%(matrix(var_Beta[,j], P, P) + Beta[,j]%*%t(Beta[,j]))))))
    log_p_beta0 = -M*log(sigma_beta0^2) - sum((1/(2*sigma_beta0^2))*((1/(N + 1/sigma_beta0^2)) + beta0^2))
    log_p_gamma_theta = sum(theta*expected_gamma) + sum(log(logistic(xi) + 10^(-10))) - 0.5*sum(theta + xi) - 
      sum(t(apply(matrix(1:P, ncol = P), 2, function(i) sum(lambda_xi_func(xi[i,])*(Sigma_theta[term_ind[i,i],] + theta[i,]^2 - xi[i,]^2)))))
    log_p_theta_Sigma_Inv = (-0.5)*sum(diag(Sigma_Inv%*%term))
    log_p_Sigma_Inv = (-0.5)*sum(diag(diag(P)%*%Sigma_Inv)) - ((P^2)/2)*log(2) - sum(apply(matrix(1:P,nrow=P), 1, function(p) gammaln((P+1-p)/2))) - (P/2)*log(P)
    Eta = cbind(rep(1, N), X) %*% rbind(beta0, Beta) 
    log_q_z = (-0.5)*sum(apply(matrix(1:M, nrow = M), 1, function(j) t(Eta[,j])%*%Eta[,j])) - 
      sum(Y*log(structure(pnorm(Eta, mean = 0, sd = 1), dim = dim(Eta)) + 10^(-10)) + 
            (1-Y)*log(1-structure(pnorm(Eta, mean = 0, sd = 1), dim = dim(Eta)) + 10^(-10)))
    log_q_beta = (-0.5)*sum(apply(matrix(1:M, nrow = M), 1, function(j) log(det(matrix(var_Beta[,j], P, P)) + 10^(-10)))) - (P*M)/2
    log_q_beta0 = - (M/2)*log((1/(N + 1/sigma_beta0^2))) -(M/2) 
    log_q_gamma = sum(expected_gamma*log(expected_gamma + 10^(-10)) + (1-expected_gamma)*log(1-expected_gamma + 10^(-10)))
    log_q_theta = (-0.5)*sum(apply(matrix(1:M, nrow = M), 1, function(j) sum(log(det(matrix(Sigma_theta[,j], P, P) + 10^(-10)))))) - (P*M)/2 
    log_q_Sigma_Inv = -(M + P)/2 - (((M+P)*P)/2)*log(2) - ((M+P)/2)*log(det(solve(diag(P) + term))) - sum(apply(matrix(1:P,nrow=P), 1, function(p) gammaln((M+P+1-p)/2)))
    Eta = 0
    
    Q = log_p_y_z + log_p_z_beta_beta0 + log_p_beta_gamma + log_p_beta0 + log_p_gamma_theta + log_p_theta_Sigma_Inv + log_p_Sigma_Inv -
      log_q_z - log_q_beta - log_q_beta0 - log_q_gamma - log_q_theta - log_q_Sigma_Inv
    
    diff = (Q - Q_old)
    ELBO = append(ELBO, Q)
    
    # Save parameter values every 100 iterations of optimization. 
    if(counter %% 100){
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
      save(params, file = sprintf("%s%s/params.RData", path_DPE, DPE_list[sim]))
      rm(params)
    }
    
  }
  
  # Evaluate time needed until convergence was reached.
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
  save(params, file = sprintf("%s%s/params.RData", path_DPE, DPE_list[sim]))
  
  return(params)
}

# Function to calculate marginal posterior of gamma. 
calculate_marginal = function(result, params0, X, y){

  # Inputs:
  # - result: list containing the output of DPE
  # - params0: list containing intial parameter values
  # - X: input data (scalar covariates)
  # - y: output data (image -> lesion masks)
  
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
  # Logistic function which avoids numerical instabilities.
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
  
  # Function to sum neighboring sparsity parameters theta together ( = sum_si_sj E[theta(s_i)]).
  sum_si_sj_function = function(j){
    x = matrix(apply(matrix(theta[,list_ind[[j]]], nrow = P), 1, sum), nrow = P)
    return(x)
  }
  
  # Set values to posterior quantities evaluated with BLESS for evaluated spike variance v0.
  beta = structure(params$Beta, dim = c(P, M))
  beta0 = params$beta0
  theta = structure(params$theta, dim = c(P, M))
  xi = structure(params$xi, dim=c(P, M))
  Sigma_Inv = matrix(params$Sigma_Inv, P, P)
  exp_gamma = structure(params$expected_gamma, dim=c(P, M))
  # Acquire binary values after thresholding posterior gamma values.
  gamma = exp_gamma
  gamma[exp_gamma > 0.5] = 1
  gamma[exp_gamma <= 0.5] = 0
  
  # Sum neighboring sparsity parameters theta together ( = sum_si_sj E[theta(s_i)]).
  sum_si_sj = apply(matrix(1:M, nrow = M), 1, sum_si_sj_function)
  # Acquire E[sum_si_sj [[theta(s_i) - theta(s_j)][theta(s_i) - theta(s_j)]^T].
  term = matrix((theta[,A$x] - theta[,A$y]), P)%*%t(matrix((theta[,A$x] - theta[,A$y]), P))
  
  # Evaluate marginal posterior of gamma under prior v0 = 0.
  marginal = 0
  for(j in 1:M){
    # X_gamma
    X_gamma = X*matrix(rep(gamma[,j], N), nrow = N, ncol = P, byrow = T)
    X_gamma = X_gamma[,apply(X_gamma, 2, function(x) !all(x == 0))]
    Q = sum(gamma[,j])
    
    # Q=0 (all gamma = 0 for voxel s_j) vs Q>0 (at least one covariate has gamma = 1 f or voxel s_j)
    if(Q == 0){
      eta = rep(beta0[j], N)
      Phi = structure(pnorm(-eta, mean = 0, sd = 1), dim = dim(eta))
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
        0.5*t(eta)%*%eta + sum(y[,j]*log(1-Phi +(10^(-10))) + (1-y[,j])*log(Phi +(10^(-10)))) +
        0.5*log(det(var_beta)+(10^(-10))) + (Q/2)
    }
  }
  
  marginal = marginal - (N/2)*sum(((1/(N + 1/sigma_beta0^2)) + beta0^2))  -
    sum((1/(2*sigma_beta0^2))*((1/(N + 1/sigma_beta0^2)) + beta0^2)) +
    sum(gamma*log(logistic(theta) + 10^(-10)) + (1-gamma)*log(1 - logistic(theta) + 10^(-10))) - ((M*P)/2)*log(2*pi) +
    (M/2)*log(det(Sigma_Inv)) - 0.5*sum(diag(Sigma_Inv%*%term)) +
    ((v - P - 1)/2)*log(det(Sigma_Inv)) - ((v*P)/2)*log(2) - ((P*(P-1))/4)*log(pi) -
    sum(apply(matrix(1:P, nrow = P), 1, function(p) gammaln((v + 1 - p)/2))) - (v/2)*log(P) -
    0.5*sum(diag(Sigma_Inv%*%diag(P))) + (M/2) + 0.5*sum(log((1/(N + 1/sigma_beta0^2)) + beta0^2)) -
    (M/2)*log(sigma_beta0^2)
  return(marginal)
}

################
### Analysis ###
################

# Run backwards DPE procedure.
for (sim in 1:n_sim_v0){
  
  # Create directory for evaluated spike variance in backwards DPE procedure to store parameter estimates in.
  dir.create(paste0(path_DPE, DPE_list[sim]))
   
  # Intialization of DPE procedure (first run of backwards DPE)
  if(DPE_list[sim] == n_sim){
  
    if(init == 'random'){
      v0_list = exp(seq(-10, -1, length.out = n_sim))
      params0 = list()
      params0$beta0 = rnorm(M, 0, 1)
      params0$Beta = matrix(rnorm(P*M, 0, 1), nrow = P, ncol = M)
      params0$Sigma_Inv = 1*diag(P)
      params0$theta = matrix(rnorm(P*M, 0, 1), nrow = P, ncol = M)
      params0$v0 = max(v0_list)
      params0$v1 = v1 
      params0$v = P
      params0$xi = matrix(rnorm(P*M, 0, 1), nrow = P, ncol = M)
      params0$sigma_beta0 = sigma_beta0
      params0$expected_gamma = matrix(0.5, P, M)
      params0$Q = -3000000000000000000000000
    }
    
    if(init == 'Firth'){
      params0 = list()
      x = as.numeric(unlist(data.table::fread(paste0(path_Firth, 'beta0.csv'), header = T)[,2]))
      x[is.na(x)] = 0
      params0$beta0 = x
      x = data.matrix(data.table::fread(paste0(path_Firth, 'Beta.csv'), header = T)[,2:(M+1)])
      x[is.na(x)] = 0
      params0$Beta = x
      rm(x)
      x = data.matrix(data.table::fread(paste0(path_Firth, 'test_statistic.csv'), header = T)[,2:(M+1)])
      x[abs(x)>1.96] = 1
      x[is.na(x)] = 0
      x[x!=1] = 0
      params0$theta = matrix(rep(apply(x, 1, mean), M), P, M)
      rm(x)
      params0$Sigma_Inv = 1*diag(P)
      params0$expected_gamma = matrix(0.5, P, M)
      params0$v0 = max(v0_list)
      params0$v1 = v1
      params0$v = P
      params0$xi = matrix(sqrt(params0$theta^2), nrow = P, ncol = M)
      params0$sigma_beta0 = sigma_beta0
      params0$Q = -3000000000000000000000000
    }
    
    params = estimate_BLESS(X, Y, params0, eps)
    marginal = calculate_marginal(params, params0, X, Y)
    params$v0_list = v0_list
  
  # Runs after intialization of DPE procedure utilizing previous parameter estimates as warm-start.
  } else{
    load(sprintf("%s%s/params.RData", path_DPE, (DPE_list[sim]+1)))
    params0 = params
    rm(params)
    params0$v0 = v0_list[(DPE_list[sim])]
    params0$v1 = v1
    params0$v = P
    params0$sigma_beta0 = sigma_beta0
    params = estimate_BLESS(X, Y, params0, eps)
    marginal = calculate_marginal(params, params0, X, Y)
    params$v0_list = v0_list
  }
  
  # Save output.
  write.csv(marginal, sprintf("%s%s/marginal.csv", path_DPE, DPE_list[sim]))
  save(params, file = sprintf("%s%s/params.RData", path_DPE, DPE_list[sim]))

}
