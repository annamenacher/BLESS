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
library(oro.nifti)
library(LaplacesDemon)

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
# Convergence threshold of difference in parameter estimates.
eps = 0.001
# Initialization of variational infrence algorithm (options: 'DPE', 'Firth', 'random').
init = 'DPE'

# Hyperparameter: Spike variance 
v0 = exp(seq(-10,log(0.05),length.out=5))[1]
# Hyperparameter: Slab variance
v1 = 10
# Hyperparameter: degrees of freedom for Wishart prior over precision matrix in MCAR prior 
v = P
# Hyperparameter: standard deviation in prior over spatially-varying intercept
sigma_beta0 = 10
# Concentration parameter of Dirichlet prior. 
dir_alpha = 1
# Number of bootstrap samples to draw.
B_bootstrap = 1000
# Concentration parameter of Dirichlet prior to draw weights for every subject in the sample. 
alpha = rep(dir_alpha,N)

# Set paths.
# path: folder with data 
path_data = ""
# path: folder with mask
path_mask = ""
# path: folder with DPE results to use as intialization for BB-BLESS
path_DPE = ""
# path: folder with Firth results (if used for intialization)
path_Firth = ""
# path: output folder for BB-BLESS
path_BB = ""
# path: folder to save parameter estimates as nifti files to later use with neuroimaging software FSL 
# to acquire cluster size images.
path_BB_results = ""


# Read in indices of mask, adjacency indices matrix & vector containing number of adjacent neighbors.
ind_mask = as.numeric(unlist(read.csv(paste0(path_mask, "ind_mask.csv"))[,2]))
A = data.frame(read.csv(paste0(path_mask, "A_mask.csv"))[,2:3])
colnames(A) = c('x','y')
n_sj = as.numeric(unlist(read.csv(paste0(path_mask, "n_sj_mask.csv"))[,2]))
load(paste0(path_mask, "list_ind_mask.RData"))

# Total number of voxels after masking lesion masks.
M = length(ind_mask)

# Data
# X (standardise data)
X = data.matrix(data.table::fread(paste0(path_data, "X_standardised_N", N, ".csv"), header = T)[,2:(P+1)])

# Y (binary lesion masks: format NxM matrix where every 3D lesion mask (.nii.gz file was vectorised and masked)
Y = data.matrix(data.table::fread(paste0(path_data, "Y_mask_N", N, ".csv"), header = T)[,2:(M+1)])


if(init == 'random'){
  v0_list = exp(seq(-10,-1, length.out = n_sim))
  params0 = list()
  params0$beta0 = rnorm(M,0,1)
  params0$Beta = matrix(rnorm(P*M,0,1), nrow = P, ncol = M)
  params0$Sigma_Inv = 1*diag(P)
  params0$theta = matrix(rnorm(P*M,0,1), nrow = P, ncol = M)
  params0$v0 = max(v0_list)
  params0$v1 = v1 
  params0$v = P
  params0$xi = matrix(rnorm(P*M,0,1), nrow = P, ncol = M)
  params0$sigma_beta0 = sigma_beta0
  params0$expected_gamma = matrix(0.5,P,M)
  params0$Q = -3000000000000000000000000
}

if(init == 'Firth'){
  
  params0 = list()
  x = as.numeric(unlist(read.csv(paste0(path_Firth, 'beta0.csv'))[,2]))
  x[is.na(x)] = 0
  params0$beta0 = x
  x = data.matrix(read.csv(paste0(path_Firth, 'Beta.csv'))[,2:(M+1)])
  x[is.na(x)] = 0
  params0$Beta = x
  rm(x)
  params0$Sigma_Inv = 1*diag(P)
  params0$theta = matrix(rnorm(P*M,0,1), nrow = P, ncol = M)
  params0$v0 = v0
  params0$v1 = v1
  params0$v = P
  params0$xi = matrix(rnorm(P*M,0,1), nrow = P, ncol = M)
  params0$sigma_beta0 = sigma_beta0
  params0$expected_gamma = matrix(0.5,P,M)
  params0$Q = -3000000000000000000000000
}

if(init == 'DPE'){
  load(paste0(path_DPE, 'params.RData'))
  
  params0 = list()
  params0$beta0 = params$beta0
  params0$Beta = params$Beta
  params0$Sigma_Inv = params$Sigma_Inv
  params0$theta =   params$theta
  params0$v0 = v0
  params0$v1 = v1
  params0$v = P
  params0$xi = params$xi 
  params0$sigma_beta0 = sigma_beta0
  params0$expected_gamma = params$expected_gamma
  params0$Q = -3000000000000000000000000
}


#################
### Functions ###
#################

# Logit function
logit = function(x){
  x = log(x/(1-x))
  return(x)
}

# Function to estimate parameter estimates with BB-BLESS.
estimate_BB_BLESS = function(X, Y, params0, eps){
  
  # Time function call.
  time.start = Sys.time()
  # Hyperparameters
  # spike & slab variance parameter
  v0 = params0$v0
  v1 = params0$v1
  # degrees of freedom of Wishart prior on covariance matrix of MCAR prior
  v = params0$v
  # Term for evaluation of Wishart covariance matrix.
  term_ind = matrix(1:(P*P),P,P)
  
  # Draw N dirichlet weights to re-weight likelihood.
  w = N*c(rdirichlet(1,alpha))
  # Draw perturbation of mean in spike-and-slab prior.
  mu = rnorm(M*P, 0, sqrt(v0))
  mu = matrix(mu, P,M)
  
  # Functions & other quantities
  # Exp-normalize trick to avoid numerical issues with large exponential values.
  exp_normalize = function(x){
    b = max(x)
    y = exp(x-b)
    return(y / sum(y))
  }
  
  # Logistic function which avoids numerical instabilities.
  logistic = function(x){
    x = ifelse(x >= 0,1 / ( 1 + exp(-x) ), exp(x) / ( 1 + exp(x) ))
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
    x = solve(XWX + diag(expected_gamma2,P))
    return(x)
  }
  
  # Function to update spatially-varying parameters.
  Beta_function = function(j){
    x = matrix(var_Beta[,j],P,P) %*%(XW%*%(expected_Z[,j] - beta0[j])+ diag(expected_gamma2[,j],P)%*%matrix(mu[,j],P,1))
    return(x)
  }
  
  # Function to update gamma values (probability of inclusion). 
  gamma_function = function(j){
    active = (-0.5*log(v1) - (matrix(diag(matrix(var_Beta[,j],P,P)), nrow = P, ncol = 1) + Beta[,j]^2- 2*Beta[,j]*mu[,j] + mu[,j]^2)/(2*v1) + theta[,j])
    not_active = (-0.5*log(v0) - (matrix(diag(matrix(var_Beta[,j],P,P)), nrow = P, ncol = 1) + Beta[,j]^2- 2*Beta[,j]*mu[,j] + mu[,j]^2)/(2*v0))
    x = apply(matrix(1:P, nrow=P), 1, function(p) exp_normalize(c(active[p],not_active[p]))[1])
    return(x)
  }
  
  # Function to sum neighboring values together.
  sum_si_sj_function = function(j){
    x = matrix(apply(matrix(theta[,list_ind[[j]]],nrow=P),1,sum),nrow=P)
    return(x)
  }
  
  # Function to update variational parameter xi.
  xi_function = function(j){
    x = sqrt((diag(matrix(Sigma_theta[,j],P,P) + theta[,j]%*%t(theta[,j]))))
    return(x)
  }
  
  # Initialize counter for number of iterations of optimization.
  counter = 0
  # Initialize difference in ELBO values to large value.
  diff=100
  # Set parameters & hyperparameters to initial values.
  sigma_beta0 = params0$sigma_beta0
  beta0=params0$beta0
  Beta=params0$Beta
  Sigma_Inv = params0$Sigma_Inv
  theta = params0$theta
  xi = params0$xi
  expected_gamma = params0$expected_gamma
  expected_gamma2 = (expected_gamma/v1) + ((1-expected_gamma)/v0)
  ELBO = NA
  Q = -3000000000000000000000000000
  # Quantities that are not memory intensive to store and evaluated in every iteration of the optimization below.
  XW = t(sweep(X,1,w,FUN=`*`))
  XWX = t(sweep(X,1,sqrt(w),FUN=`*`))%*%sweep(X,1,sqrt(w),FUN=`*`)
  
  # Run optimization until convergence criteria is reached.
  while(diff>eps){ 
    
    # Avoid memory issues by clearing unnecessairly loaded quantities every iteration of the optimization.
    gc()    
    
    # Increase number of iterations by 1.
    counter = counter + 1
    # Save old quantities. 
    Beta_old = Beta
    Q_old = Q
    Beta_old = Beta
 
    # Update Z.
    expected_Z = structure(hutils::if_else(Y == 1,
                                           (cbind(rep(1,N), X) %*% rbind(beta0, Beta))  +
                                             (structure(dnorm(-(cbind(rep(1,N), X) %*% rbind(beta0, Beta)) , mean = 0, sd=1),dim=c(N,M)) /
                                                (structure(pnorm((cbind(rep(1,N), X) %*% rbind(beta0, Beta)), mean = 0, sd=1)+ 10^(-10),dim=c(N,M)))) ,
                                           (cbind(rep(1,N), X) %*% rbind(beta0, Beta)) - (structure(dnorm(-(cbind(rep(1,N), X) %*% rbind(beta0, Beta)), mean = 0, sd=1),dim=c(N,M)) /
                                                                                            (1 - structure(pnorm((cbind(rep(1,N), X) %*% rbind(beta0, Beta)), mean = 0, sd=1),dim=c(N,M))+ 10^(-10))) ),
                           dim=c(N,M))
    
    # If this quantity diverges to infinity then set to really large or small value.
    expected_Z[expected_Z == -Inf] = -15
    expected_Z[expected_Z == Inf] = 15
    
    # Update beta0.
    beta0 = (1/(sum(w) + 1/sigma_beta0^2))*(as.vector(matrix(w,1,N)%*%(expected_Z - X%*%Beta)))
    
    # Update variance of Beta.
    var_Beta  = apply(expected_gamma2, 2, FUN=var_beta_function)
    
    # Update beta.
    Beta=apply(matrix(1:M,nrow=M),1,Beta_function)
    
    # Update gamma.
    expected_gamma = apply(matrix(1:M,nrow=M),1,gamma_function)
    
    # Update expectation: E[gamma*(1/v1) + (1-gamma)*(1/v0)]
    expected_gamma2 = (expected_gamma/v1) + ((1-expected_gamma)/v0)
    
    # Update theta
    sum_si_sj = apply(matrix(1:M,nrow=M),1,sum_si_sj_function)
    Sigma_theta = apply(matrix(1:M,nrow=M),1,function(j) as.vector(solve(n_sj[j]*Sigma_Inv + 2*diag(c(lambda_xi_func(xi[,j])),P))))
    theta = apply(matrix(1:M,nrow=M),1,function(j) matrix(Sigma_theta[,j],P,P)%*%(Sigma_Inv%*%sum_si_sj[,j] + expected_gamma[,j] - 1/2))
    
    # Update Sigma Inverse.
    term = matrix((theta[,A$x] - theta[,A$y]),P)%*%t(matrix((theta[,A$x] - theta[,A$y]),P)) + matrix(rowSums(Sigma_theta[,A$y]),P)
    Sigma_Inv = (M + v - P - 1)*solve(diag(P) + term)
    
    # Update xi.
    xi = apply(matrix(1:M, nrow=M), 1, xi_function)
    
    # too costly to evaluate ! (slows down optimization)
    # Calculate ELBO. 
    # expected_Z = 0
    # log_p_y_z = 0
    # log_p_z_beta_beta0 = (-0.5)*sum(apply(matrix(1:M,nrow=M),1,function(j) sum(diag(t(X)%*%diag(w,N)%*%X%*%(matrix(var_Beta[,j],P,P) + Beta[,j]%*%t(Beta[,j])))))) - 
    #   (sum(w)/2)*sum((1/(sum(w) + 1/sigma_beta0^2)) + beta0^2) -
    #   sum(apply(matrix(1:M,nrow=M),1,function(j) t(Beta[,j])%*%t(diag(w,N)%*%X)%*%matrix(rep(beta0[j],N),N,1)))
    # log_p_beta_gamma = (-0.5)*sum(expected_gamma*log(v1) + (1-expected_gamma)*log(v0)) - 
    #   0.5*sum(apply(matrix(1:M,nrow=M),1,function(j) sum(diag(diag(expected_gamma2[,j],P)%*%(matrix(var_Beta[,j],P,P) + Beta[,j]%*%t(Beta[,j])))))) -
    #   0.5*sum(apply(matrix(1:M,nrow=M),1,function(j) sum(diag(diag(expected_gamma2[,j],P)%*%matrix(mu[,j],P,1)%*%t(matrix(mu[,j],P,1)))))) +
    #   sum(apply(matrix(1:M,nrow=M),1,function(j) sum(diag(diag(expected_gamma2[,j],P)%*%matrix(mu[,j],P,1)%*%t(Beta[,j])))))
    # log_p_beta0 = -(M/2)*log(sigma_beta0^2) - sum((1/(2*sigma_beta0^2))*((1/(N + 1/sigma_beta0^2)) + beta0^2))
    # log_p_gamma_theta = sum(theta*expected_gamma) + sum(log(logistic(xi) + 10^(-10))) - 0.5*sum(theta + xi) - 
    #   sum(t(apply(matrix(1:P,ncol=P),2,function(i) sum(lambda_xi_func(xi[i,])*(Sigma_theta[term_ind[i,i],] + theta[i,]^2 - xi[i,]^2)))))
    # log_p_theta_Sigma_Inv = (-0.5)*sum(diag(Sigma_Inv%*%term))
    # log_p_Sigma_Inv = (-0.5)*sum(diag(diag(P)%*%Sigma_Inv)) - ((P^2)/2)*log(2) - sum(apply(matrix(1:P,nrow=P),1, function(p) gammaln((P+1-p)/2))) - (P/2)*log(P)
    # Eta = cbind(rep(1,N), X) %*% rbind(beta0, Beta) 
    # log_q_z = (-0.5)*sum(apply(matrix(1:M,nrow=M),1,function(j) t(Eta[,j]*w)%*%Eta[,j])) - 
    #   sum(diag(w,N)%*%Y*log(structure(pnorm(Eta, mean = 0, sd=1),dim=dim(Eta)) + 10^(-10)) + diag(w,N)%*%(1-Y)*log(1-structure(pnorm(Eta, mean = 0, sd=1),dim=dim(Eta)) + 10^(-10)))
    # log_q_beta = (-0.5)*sum(apply(matrix(1:M,nrow=M),1,function(j) log(det(matrix(var_Beta[,j],P,P)) + 10^(-10)))) - (P*M)/2
    # log_q_beta0 = - (M/2)*log((1/(N + 1/sigma_beta0^2))) -(M/2) 
    # log_q_gamma = sum(expected_gamma*log(expected_gamma + 10^(-10)) + (1-expected_gamma)*log(1-expected_gamma + 10^(-10)))
    # log_q_theta = (-0.5)*sum(apply(matrix(1:M,nrow=M),1,function(j) sum(log(det(matrix(Sigma_theta[,j],P,P) + 10^(-10)))))) - (P*M)/2 
    # log_q_Sigma_Inv = -(M + P)/2 - (((M+P)*P)/2)*log(2) - ((M+P)/2)*log(det(solve(diag(P) + term))) - sum(apply(matrix(1:P,nrow=P),1, function(p) gammaln((M+P+1-p)/2)))
    # 
    # Q = log_p_y_z + log_p_z_beta_beta0 + log_p_beta_gamma + log_p_beta0 + log_p_gamma_theta + log_p_theta_Sigma_Inv + log_p_Sigma_Inv -
    #   log_q_z - log_q_beta - log_q_beta0 - log_q_gamma - log_q_theta - log_q_Sigma_Inv
    
    #diff = (Q - Q_old)
    Eta = 0
    diff = max(abs(Beta - Beta_old))
    #ELBO= append(ELBO, Q)
    
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
      params$diff = diff
      params$final = F
      save(params, file = sprintf("%sparams%04d.RData", path_BB, sim))
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
  params$final = T
  save(params, file = sprintf("%sparams%04d.RData", path_BB, sim))
  
  return(params)
}

#####################
### Data Analysis ###
#####################

# Note: This for loop can be run in parallel! Every optimization is independent of each other.
for(sim in 1:B){

# Run BB-BLESS optimaization.
params = estimate_BB_BLESS(X,Y,params0,eps)
params$v0 = v0
params$v1 = v1

# Save parameter estimates as .RData object
save(params, file = sprintf("%sparams%04d.RData", path_BB, sim))

# Save parameter estimates for covariate age as 3D image in the nifti format.
mask = readNIfTI(paste0(path_mask, 'Emp_map_2mm_mask.nii.gz'))
img_NA = mask
mask = as.vector(mask)
ind1 = which(mask == 1, arr.ind = T)

img_filler = rep(NA, dim1*dim2*dim3)
img_filler[ind1] = as.numeric(unlist(params$Beta[1,]))
img_NA@.Data = array(img_filler, c(dim1,dim2,dim3))
writeNIfTI(img_NA, sprintf('%sBeta%04d', path_BB_results, sim))
rm(img_NA)

}
