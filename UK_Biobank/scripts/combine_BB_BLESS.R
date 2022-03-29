##############
### Script ###
##############

# Run code to combine bootstrap samples into matrices to acquire sample from posterior distribution.
# (Code for UK Biobank application)

#################
### Libraries ###
#################

library(oro.nifti)

#################
### Constants ###
#################

# Set paths.
# path: folder with BB-BLESS samples
path = ''
# path: folder for concatenated beta bootstrap samples and  posterior mean, posterior standard deviation and 
# standardized effects after concatenating results.
path_out = ''
# path: folder for standardized effect maps -> needed to acquire cluster size maps 
path_t = ''
# path: folder with mask & spatial adjacency list & vector with number of adjacent neighbors
path_mask = ""
# path: folder with nifti parameter maps
path_nifti = ''

# Read in indices of mask, adjacency indices matrix & vector containing number of adjacent neighbors.
ind_mask = as.numeric(unlist(read.csv(paste0(path_mask, "ind_mask.csv"))[,2]))

# Number of bootstrap samples to draw.
B = 1000
# Total number of voxels after masking lesion masks.
M = length(ind_mask)

# Create BxM matrix (number of bootstrap samples x number of voxels) to store bootstrap samples for parameter map.
Beta = matrix(NA, B, M)

# Concatenate bootstrap samples.
for(b in 1:B){
  tryCatch({
    # In case of non-convergence of a bootstrap sample continue concatenation.
    img = readNIfTI(sprintf('%sBeta%04d.nii.gz', path, b))
    Beta[b,] = as.vector(img@.Data)[ind]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

# Only keep samples that have converged.
Beta = Beta[complete.cases(Beta),]
# Calculate posterior mean, posterior standard deviation and standardized effect for BB-BLESS. 
beta_mean = apply(Beta, 2, mean)
beta_sd = apply(Beta, 2, sd)
beta_t = beta_mean / beta_sd

# Save concatenated bootstrap samples.
write.csv(Beta, paste0(path_out, 'Beta.csv'))
# Save posterior quantities.
write.csv(beta_mean, paste0(path_out, 'beta_mean.csv'))
write.csv(beta_sd, paste0(path_out, 'beta_sd.csv'))
write.csv(beta_t, paste0(path_out, 'beta_t.csv'))

# Save test statistic map in nifti format.
mask = readNIfTI(paste0(path_mask, 'Emp_map_2mm_mask.nii.gz'))
img_NA = mask
img_filler = rep(NA, dim1*dim2*dim3)
img_filler[ind] = as.numeric(unlist(beta_t))
img_NA@.Data = array(img_filler, c(dim1,dim2,dim3))
writeNIfTI(img_NA, sprintf('%st_final', path_t))
rm(img_NA)

t = sweep(Beta, 2, beta_sd, FUN = '/')
rm(Beta)

# Acquire nifti file for every bootstrap standardized effect map.
for(b in 1:nrow(t)){
  mask = readNIfTI(paste0(path_mask, 'Emp_map_2mm_mask.nii.gz'))
  img_NA = mask
  img_filler = rep(NA, dim1*dim2*dim3)
  img_filler[ind] = as.numeric(unlist(t[b,]))
  img_NA@.Data = array(img_filler, c(dim1,dim2,dim3))
  writeNIfTI(img_NA, sprintf('%s/t_%s', path_t, b))
  rm(img_NA)
}

# Save concatenated boostrap test statistics.
write.csv(t, paste0(path_t, 'test_statistics.csv'))







