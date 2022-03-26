##############
### Script ###
##############

# Generate 2D lesion masks for N = 100,000 (very large sample size, close to population size) 
# for every base rate intensity (lambda = {1,2,3}) in order to invert lesion probabilities to 
# acquire true parameters.

#################
### Libraries ###
#################

library(deldir)
library(spatstat)
library(MASS)
library(hutils)
library(logistf)
library(scales)
library(dplyr)
library(reshape2)

#################
### Constants ###
#################

# Base rate intensity (lambda = {1,2,3})
lambda = 1
# Effect size.
size = 1
# Number of covariates.
P = 2

path = ""
M = 2500
dimension = 50

# female:effect on 50% of image (right)
# group 2: effect on 25% of lower left of image
setting_1f = c(1, 4, 4, 1)
setting_1m = c(1, 1, 1, 1)
setting_2f = c(1, 4, 4, 4)
setting_2m = c(1, 1, 1, 4)

#################
### Functions ###
#################

Lesion_Prob_Matrix_Setup = function(Prob_hat, N_scenario){
  Prob_hat = colSums(Prob_hat) / N_scenario
  # Convert to matrix with shape of dimension x dimension for plotting.
  Lesion_prob_matrix = matrix(Prob_hat, ncol = dimension, nrow = dimension)
  return(Lesion_prob_matrix)
}

lesion_number = function(lambda, setting){
  # Generates the number of lesions in each quadrant of the image matrix.
  # n_l ~ Poisson(lambda*setting)
  #     -> lambda: average number of lesions in each quadrant
  #     -> setting: increased lesion occurance for group or sex
  n_l = numeric(4)
  for(i in 1:4){
    n_l[i] = rpois(1, lambda * setting[i])
  }
  return(n_l)
}


quadrant_map = function(n_l, M){
  # Generates a lesion map for a single quadrant.
  #   -> n_l: number of lesion occurances in a quadrant
  #   -> M: number of voxels in the total image (4 * number of voxels in a quadrant)

  # Set of lesion sizes
  lesion_size = c(1,3,5)
  # Lesion location ~ Uniform anywhere within the quadrant.
  location = floor(runif(n_l, min=1, max = M/4))
  # For every lesion location -> determine the size of the lesion (uniform distributed across the set)
  size = sample(lesion_size, n_l, replace = TRUE)
  # Create empty image vector.
  image = numeric(M/4)
  for(i in 1:n_l){
    image[location[i]] = size[i]
  }
  # Convert image vector to matrix (columnwise, dim=(sqrt(M/4), sqrt(M/4)))
  image = matrix(image, nrow = sqrt(M/4))
  return(image)
}

combine_quadrants = function(image_q1, image_q2, image_q3, image_q4){
  # Combines image quadrants that capture 4 different scenarios (lesion intensities).
  #    -> Input: 4x image with sqrt(M/4) x sqrt(M/4) dimensions
  #    -> Output: 1x image with sqrt(M/4) x sqrt(M/4) dimensions
  upper_half = cbind2(image_q1, image_q2)
  lower_half = cbind2(image_q4, image_q3)
  image = rbind2(upper_half, lower_half)
  return(image)
}
neighbor_square <- function(image, i, j, size) {
  # Identifies the size of each lesion and acquires indices across the image matrix.
  # -> Note: Lesions of various sizes are able to overlap.
  if(size == 1){
    idx = matrix(c(i, j), ncol = 2)
  } else if(size == 3){
    idx = matrix(c(i-1, i-1, i-1, i, i, i, i+1, i+1, i+1, j-1, j,
                   j+1, j-1, j, j+1, j-1, j, j+1), ncol = 2)
    # set out of bound indices to 0
    idx[idx[, 1] > nrow(image), 1] = 0
    idx[idx[, 2] > ncol(image), 2] = 0
    idx = idx[idx[,1] > 0,]
    idx = idx[idx[,2] > 0,]
  } else{
    idx = matrix(c(i-2, i-2, i-2, i-2, i-2, i-1, i-1, i-1, i-1, i-1, i, i, i, i, i, i+1, i+1, i+1, i+1, i+1, i+2, i+2, i+2, i+2, i+2,
                   j-2, j-1, j, j+1, j+2, j-2, j-1, j, j+1, j+2, j-2, j-1, j, j+1, j+2, j-2, j-1, j, j+1, j+2, j-2, j-1, j, j+1, j+2), ncol = 2)
    # set out of bound indices to 0
    idx[idx[, 1] > nrow(image), 1] = 0
    idx[idx[, 2] > ncol(image), 2] = 0
    idx = idx[idx[,1] > 0,]
    idx = idx[idx[,2] > 0,]
  }
  return(idx)
}


config_image = function(image, n_l){
  # Configure image and fill in the lesions accrodingly as identified by the
  # function "neighbor_square()".
  idx = which(image !=0,arr.ind = T)
  n_l_total = nrow(idx)
  if(n_l_total > 0){
    for(i in 1:n_l_total){
      k = as.integer(idx[i,1])
      l = as.integer(idx[i,2])
      y = neighbor_square(image, k, l, image[k, l])
      # Might be redundant -> please check -> make lesion map binary
      image[y] = 1
    }
  }
  return(image)
}


generate_sample = function(lambda, setting, N, M){
  # Generate 4 lesion map quadrants varying according to a predefined setting
  # for all subjects N.
  # Output: N images -> NxM matrix
  Y = matrix(NA, nrow = N, ncol = M)
  for(s in 1:N){
    n_l = lesion_number(lambda = lambda, setting = setting)

    image_q1 = quadrant_map(n_l = n_l[1], M = M)
    image_q2 = quadrant_map(n_l = n_l[2], M = M)
    image_q3 = quadrant_map(n_l = n_l[3], M = M)
    image_q4 = quadrant_map(n_l = n_l[4], M = M)

    n_l = sum(n_l)

    lesion_map = combine_quadrants(image_q1, image_q2, image_q3, image_q4)
    lesion_map = config_image(lesion_map, n_l)
    Y[s,] = as.vector(lesion_map)
  }
  return(Y)
}

generate_truth_P2 = function(lambda, setting_1m, setting_2m, setting_1f, setting_2f, M){
  
dimension = sqrt(M)
result = list()

# y
y_1m = generate_sample(setting = setting_1m, lambda = lambda, N = 100000, M = M)
y_1f = generate_sample(setting = setting_1f, lambda = lambda, N = 100000, M = M)
y_2m = generate_sample(setting = setting_2m, lambda = lambda, N = 100000, M = M)
y_2f = generate_sample(setting = setting_2f, lambda = lambda, N = 100000, M = M)

# beta0 & Beta
result$beta0 = matrix((qnorm(colSums(y_1m)/100000,0,1)), nrow = dimension, ncol = dimension)
result$beta1 = matrix((qnorm(colSums(y_1f)/100000,0,1) - result$beta0), nrow = dimension, ncol = dimension)
result$beta2 = matrix((qnorm(colSums(y_2m)/100000,0,1) - result$beta0), nrow = dimension, ncol = dimension)

result$y_1m = matrix(colSums(y_1m)/100000, nrow = dimension, ncol = dimension)
result$y_1f = matrix(colSums(y_1f)/100000, nrow = dimension, ncol = dimension)
result$y_2m = matrix(colSums(y_2m)/100000, nrow = dimension, ncol = dimension)
result$y_2f = matrix(colSums(y_2f)/100000, nrow = dimension, ncol = dimension)
  
  return(result)
}

truth = generate_truth_P2(lambda,setting_1m, setting_2m, setting_1f, setting_2f, M)

save(truth, file = paste0(path, "truth.RData"))
