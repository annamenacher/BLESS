##############
### Script ###
##############

# Run code to acquire indices of spatial neighbors for masked map of voxel locations. 

# Note: Acquiring an adjacency matrix for 91x109x91 voxels would have dimensions of 902,629 x 902,629 
# which is too computationally costly to compute and most importantly infeasible to store in memory.
# Hence, the following functions are run to only acquire the indices of neighboring voxels for the masked white 
# matter locations where at least one lesion occurred at every voxel in the mask among the UK Biobank population
# for which lesion masks were avaialble. 

# Functions: 
# - adjacency_matrix 
# - n_neighbors

# (Code for UK Biobank application)

#################
### Libraries ###
#################

library(oro.nifti)

#################
### Constants ###
#################

# Dimensions of lesion masks (3D MRI scan with dimensions 91x109x91)
dim1 = 91
dim2 = 109
dim3 = 91
# Total number of voxels in a 3D MRI scan / lesion mask.
M = dim1 * dim2 * dim3 

# path: folder with empirical binary lesion mask indicating where lesions occur at least once in the population
# of the UK Biobank.
path_data = " "

#################
### Functions ###
#################

# Function for deriving indices of adjacent neighbors of every single voxel location for 
# spatial MCAR prior (a voxel location is considered adjacent if they share the same face).
adjacency_matrix = function(dim1, dim2, dim3){
    
    if(missing(dim3)){
      A = data.frame(x = integer(),y = integer())
      ind = 1:(dim1*dim2)
      conv = as.vector(matrix(1:(dim1*dim2), dim1,dim2, byrow = T))
      
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

# Function to swap indices in a vector (i.e. vec = c(1,2,3) should actually be vec = c(2,5,1)
# hence: 1 -> 2, 2 -> 5, 3 -> 1).
swap <- function(vec, from, to) {
  tmp <- to[ match(vec, from) ]
  tmp[is.na(tmp)] <- vec[is.na(tmp)]
  return(tmp)
}

################
### Analysis ###
################

# Indices of adjacency matrix of 3D lattice.
A = adjacency_matrix(dim1, dim2, dim3)

# Number of neighbors of 3D lattice.
n_sj = n_neighbors(dim1, dim2, dim3)

# Get indices of upper triangular of "adjacency matrix" for 3D lattice. 
A = upper_triangular(A, dim1*dim2*dim3)

# Make a list with every voxel location (vector: 1 to M) and its neighbor indices 
# (for example: voxel 1 in a 5x5 image will have neighbors at voxels 2 & 6, 
# and has therefore 2 neighbors)
list_ind = list()
for(j in 1:M){
  list_ind[[j]] = c(A$y[which(A$x == j, arr.ind = T)], A$x[which(A$y == j, arr.ind = T)])
}

# Save adjacency indices matrix "A", vector with number of adjacent neighbors "n_sj" and list of adjacent indices
# "list_ind".
write.csv(A, paste0(path_data, "A.csv"))
write.csv(n_sj, paste0(path_data, "n_sj.csv"))
save(list_ind, file = paste0(path_data, "list_ind.RData"))

# Read in mask containing white matter tract where at least one lesion occurs in the UK Biobank population. 
mask = oro.nifti::readNIfTI(paste0(path_data, "Emp_map_2mm_mask.nii"))
mask = mask@.Data

# Get indicies of a vectorised version of the mask where there is at least 1 lesion and where there is no lesion.
id_mask = which(as.vector(mask) == 1, arr.ind = T)
id_mask_zero = which(as.vector(mask) == 0, arr.ind = T)

# Get total number of voxels within the mask.
M = length(id_mask)

# Swap all values of indices in adjacency indices matrix "A" that are not in mask with 0. 
x = swap(A$x, id_mask_zero, rep(0,length(id_mask_zero)))
y = swap(A$y, id_mask_zero, rep(0,length(id_mask_zero)))
A = cbind(x,y)
# Delete all entries from adjacency indices matrix "A" that are not in mask.
A[A == 0] = NA
A = A[complete.cases(A),]
A = data.frame(A)
colnames(A) = c('x', 'y')

# Swap all values of indices in adjacency indices matrix "A" that are in mask (ordered from lowest to highest index) 
# with values 1:M where M is the number of total number of voxels within the mask. 
x = swap(A$x, id_mask, 1:M)
y = swap(A$y, id_mask, 1:M)
A = data.frame(cbind(x, y))
colnames(A) = c('x', 'y')

# Make a list with every voxel location (vector: 1 to M) and its neighbor indices 
# (for example: voxel 1 in a 5x5 image will have neighbors at voxels 2 & 6, 
# and has therefore 2 neighbors)
list_ind = list()
n_sj = numeric(M)
for(j in 1:M){
    list_ind[[j]] =  c(A$y[which(A$x == id_mask[j], arr.ind = T)], A$x[which(A$y == id_mask[j], arr.ind = T)])
  n_sj[j] = length(list_ind[[j]])
  }

# Save masked versions of adjacency indices matrix "A", vector with number of adjacent neighbors "n_sj" and 
# list of adjacent indices "list_ind".
save(list_ind, file = paste0(path_data, 'list_ind_mask.RData'))
write.csv(n_sj, paste0(path_data, "n_sj_mask.csv"))
write.csv(A, paste0(path_data, "A_mask.csv"))
