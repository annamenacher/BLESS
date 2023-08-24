##############
### Script ###
##############

# Run Firth Regression at each voxel independently. Due to the large data size this script needs to be run
# in parallel without the for loop at each voxel independently.
# (Code for UK Biobank application)

#################
### Libraries ###
#################

library(brglm2)
library(data.table)

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
M_full = dim1 * dim2 * dim3 
# Number of covariates (age, sex, headsize scaling factor, age by sex).
P = 4

# Set paths.
# path: folder containing data
path = path_data = ""
# path: folder with mask 
path_mask = ""

# Read in indices of mask
ind_mask = as.numeric(unlist(read.csv(paste0(path_mask, "ind_mask.csv"))[,2]))
# Total number of voxels after masking lesion masks.
M = length(ind_mask)

# Data
# X (standardise data)
X = data.matrix(data.table::fread(paste0(path_data, "X_standardised_N", N, ".csv"), header = T)[,2:(P+1)])

# Y (binary lesion masks: format NxM matrix where every 3D lesion mask (.nii.gz file was vectorised and masked)
Y = data.matrix(data.table::fread(paste0(path_data, "Y_mask_N", N, ".csv"), header = T)[,2:(M+1)])

#####################
### Data Analysis ###
#####################

# Setup matrices to store parameter estimates.
params_Firth = matrix(0, P+1, M)
var_Firth = matrix(0, P+1, M)

# Run Firth regression for each voxel independently.
for(j in 1:M){
  tryCatch({
    data = data.frame(cbind(Y[,j], X))
    colnames(data) = c('y','x1','x2','x3','x4')
    
    model = glm(y ~ 1 + x1 + x2 + x3 + x4, family = binomial(probit), data = data, method = "brglmFit", type = 'AS_mean', maxit = 10000, epsilon = 1e-05, slowit = 0.5)
                  model = summary(model)
    params_Firth[,j] = model$coefficients[,1]
    var_Firth[,j] = model$coefficients[,2]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# Save parameter estimates.
write.csv(params_Firth[1,], sprintf("%s/beta0.csv", path, sim))
write.csv(params_Firth[2:(P+1),], sprintf("%s/Beta.csv", path, sim))
write.csv(var_Firth[1,], sprintf("%s/Sigma_beta0.csv", path, sim))
write.csv(var_Firth[2:(P+1),], sprintf("%s/Sigma_Beta.csv", path, sim))
write.csv((params_Firth[2:(P+1),]/var_Firth[2:(P+1),]), sprintf("%s/test_statistic.csv", path, sim))



