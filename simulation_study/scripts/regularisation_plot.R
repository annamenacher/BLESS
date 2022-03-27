##############
### Script ###
##############

# Script to get regularization plots comparing estimated coefficients across range of spike variances in backwards DPE procedure.
# (Code for Simulation Study)

#################
### Libraries ###
#################

library(ggplot2)
library(latex2exp)

############
### Plot ###
############

# Range of spike variances DPE was evaluated over.
v0_list = exp(seq(-20,-1,length.out = 15))
# Length of spike variance range.
n_sim_v0 = length(v0_list)

# path: file containing estimated BLESS posterior mean values from spatially-varying coefficients
path_beta = ''
# path: file containing estimated BLESS posterior mean values from latent variable gamma
path_gamma = ''

# Read in coefficients & latents.
gamma1_list = data.matrix(read.csv(path_gamma)[,2:2501])
Beta1_list = data.matrix(read.csv(path_beta)[,2:2501])

result = list()
result$v0_list = v0_list
result$gamma_list = gamma1_list
result$beta_list = Beta1_list

regularisation_plot = function(result, voxels = 'some', border = F, truth = F, effect = 1){
  if(voxels == 'some'){
    p = 10
    idx = c(342,532,673,987,1109,2341,1745,1867,2034,1981)
  }
  if(voxels == 'all'){
    if(border == T){
      p = 2500
      idx = 1:2500
    }
    if(border==F){
      idx_border = c(1:100,2401:2500,1151:1350,seq(1,2451,50),seq(2,2452,50),seq(49,2499,50),seq(50,2500,50),seq(24,2474,50),seq(25,2475,50),seq(26,2476,50),seq(27,2477,50))
      idx = c(1:2500)[-idx_border]
      p = length(idx)
    }
  }
  v0_list = log(result$v0_list)
  
  gamma_list = result$gamma_list[,idx]
  beta_list = result$beta_list[,idx]
  
  if(truth == F){
    t = matrix(as.matrix(gamma_list), nrow = length(v0_list), ncol = p)
    t[t>0.5] = 1
    t[t<=0.5] = 0
  }
  if(truth == T){
    if(effect == 1){
      t = matrix(0, nrow = length(v0_list), ncol = M)
      t[,1:1250] = 0
      t[,1250:2500] = 1
      t = t[,idx]
    }
    
    if(effect == 2){
      t_eff = matrix(0,50,50)
      t_eff[26:50,1:25] = 1
      t = matrix(rep(t_eff,length(v0_list)), length(v0_list), p, byrow = T)
    }
  }
  
  plot(0,0,xlim = c((min(v0_list)),max(v0_list)),ylim = c((min(beta_list, na.rm = T) - 0.1), (max(beta_list, na.rm=T) + 0.1)), xlab =  TeX('$log(\\nu_0)$'), ylab = TeX('$\\beta$'), type = "n", cex.axis=1.15, cex.lab=1.15)
  for (j in 1:p){
    lines(v0_list,beta_list[,j], col=alpha(rgb(0.5,0.5,0.5),0.1))
    for(i in 1:length(v0_list)){
      tryCatch({
        if(t[i,j] == 1){
          points(v0_list[i],beta_list[i,j], col= 'red', pch=20)
        }
        if(t[i,j] == 0){
          points(v0_list[i],beta_list[i,j], col='blue', pch=20)
        }
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})    
    }
  }

}

regularisation_plot(result, voxels = 'all', border = T, truth = F, effect = 1)


