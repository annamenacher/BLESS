##############
### Script ###
##############

# Run independent Firth regression models for every voxel j = 1,...,M (Code for simulation study).

#################
### Libraries ###
#################

library(brglm2)

#################
### Constants ###
#################
# Number of dataset to perform analysis on.
sim = 1
# Number of voxels.
M = 2500
# Dimensions of 2D slices.
dimension = sqrt(M)
# Number of covariates.
P = 2
# Base rate intensity (lambda = {1,2,3})
lambda = 3
# Sample size (N= {500, 1000, 5000})
N = 1000

# Path: data folder
path_data = ""
# Path: output folder
path = ""

########################
### Simulation Study ###
########################

# Data
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)
Y = data.frame(read.csv(paste0(path_data, "Y", sim, ".csv"), header = T)[,2:(M+1)])

time.start = Sys.time()
# Set up empty matrix for parameter estimates ((P+1)xM), parameter variances ((P+1)xM), and predictions (4xM).
# -> predictions are fixed for every covariate combination 
# -> (group 1 - female, group 1 - male, group 2 - female, group 2 - male)
# -> P + 1: number of covariates + spatially-varying intercept
params_Firth = matrix(NA, nrow = P+1, ncol = M)
var_Firth = matrix(NA, nrow = P+1, ncol = M)
y_pred_Firth = matrix(NA, nrow = 4, ncol = M)

# Run independent Firth regression for every voxel s_j.
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
time.end = Sys.time()
(time = time.end - time.start)

# Save output.
write.csv(as.vector(params_Firth[1,]), sprintf("%sbeta0%03d.csv", path, sim))
write.csv(as.vector(params_Firth[2,]), sprintf("%sBeta1%03d.csv", path, sim))
write.csv(as.vector(params_Firth[3,]), sprintf("%sBeta2%03d.csv", path, sim))
write.csv(as.vector(var_Firth[1,]), sprintf("%sstd_error_beta0%03d.csv", path, sim))
write.csv(as.vector(var_Firth[2,]), sprintf("%sstd_error_Beta1%03d.csv", path, sim))
write.csv(as.vector(var_Firth[3,]), sprintf("%sstd_error_Beta2%03d.csv", path, sim))
write.csv(as.vector(y_pred_Firth[1,]), sprintf("%sy_pred_1f%03d.csv", path, sim))
write.csv(as.vector(y_pred_Firth[3,]), sprintf("%sy_pred_1m%03d.csv", path, sim))
write.csv(as.vector(y_pred_Firth[2,]), sprintf("%sy_pred_2f%03d.csv", path, sim))
write.csv(as.vector(y_pred_Firth[4,]), sprintf("%sy_pred_2m%03d.csv", path, sim))
write.csv(time, sprintf("%stime%03d.csv", path, sim))



