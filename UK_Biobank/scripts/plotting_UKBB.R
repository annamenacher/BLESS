########################
### Acknowledgements ###
########################

# Code for plotting 2D slices of MRI scans has kindly been provided by Petya Kindalova! Many thanks! 

#################################
### Slice plots: MRI template ###
#################################

#-----------------
## set constants

# Set file paths.
# path: path to mask
path_mask = paste0("", "Emp_map_2mm_mask.nii.gz")
# path: path to brain template
path_mni152 = paste0("", "MNI152_T1_2mm_brain.nii.gz")
# path: folder with BLESS parameter estimates
path_results_BLESS = ''
# path: folder with Firth parameter estimates
path_results_firth = ''
# path: folder for plots
path_output = ''

# Set working directory.
setwd(path_output)

v0_start = log(0.05)
# Set range of spike variance values.
v0_list = exp(seq(-10, v0_start, length.out = n_sim))
n_sim_v0 = length(v0_list)
  
# Dimensions of 3D MRI scan / lesion mask.
dimension1 = dim1 = 91
dimension2 = dim2 = 109
dimension3 = dim3 = 91 

#-----------------
## load libraries
detach("package:ggplot2", unload=TRUE)
library(oro.nifti)
library(lattice)
library(raster)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)
library(grid)

#--------------------
## plotting function
# inputs to the function
# slice             (number) currently axial but you can easily change that
# cols              colours for the coefficients
# cols_mask         colours for the brain mask in the background
# mni152            (nifti) brain for the background
# img               (nifti) with coefficients
# mask              (nifti) brain mask
# empir_prob        (nifti) empirical lesion probability in case you want to threshold based on that
# empir_threshold   (number between 0 and 1) threshold to use for empirical lesion probability
# name of pdf file  (char) here you provide the name and I have hard coded the name to end with the slice
# legend            (bool) T or F whether tou want the legend included, this means colour key here
# min and max       usually +-5 used or +-2

zscores_plot = function(slice, cols, cols_mask, 
                        mni152, img, mask, 
                        empir_prob, empir_threshold, 
                        name, legend, min, max) {
  img = img@.Data
  img = as.numeric(unlist(img))
  img = array(img, dim = c(dimension1, dimension2, dimension3))
  mni152 = mni152@.Data
  mni152 = as.numeric(unlist(mni152))
  mni152 = array(mni152, dim = c(dimension1, dimension2, dimension3))
  mask = mask@.Data
  mask = as.numeric(unlist(mask))
  mask = array(mask, dim = c(dimension1, dimension2, dimension3))
  empir_prob = empir_prob@.Data
  empir_prob = as.numeric(unlist(empir_prob))
  empir_prob = array(empir_prob, dim = c(dimension1, dimension2, dimension3))
  # 
  z_score = img
  z_score[empir_prob<=empir_threshold] = NA
  z_score[which(mask!=1)] = NA
  
  
  z_score[which(z_score>max)] = max
  z_score[which(z_score<min)] = min
  
  my.at = c(seq(min, max, length=length(cols)-1))
  my.key = list(at=my.at, col=cols, labels=list(cex=2, font=2))
  
  #, scales=list(draw=F), colorkey = my.key
  pdf(paste0(name, slice, '.pdf'))
  if(legend) {
    
    print(levelplot(mni152[,,slice], xlab="", ylab="", col.regions=cols_mask, scales=list(draw=F), colorkey = my.key, cuts=(length(cols_mask)-1), at=c(seq(min(mni152), max(mni152), length=length(cols_mask)-1)), main="", cex.main=2, useRaster=T) + 
            levelplot(z_score[,,slice], xlab="", ylab="", region=T, col.regions=cols, scales=list(draw=F), cuts=(length(cols)-1),at=c(seq(min, max, length=length(cols)-1)), useRaster=T)+layer(grid.text(plottitle, x=unit(0.03, "npc"), y=unit(0.98, "npc"), just=c("left", "top"), gp=gpar(fontsize=35, col="white"))))  
  } else {
    
    print(levelplot(mni152[,,slice], xlab="", ylab="", col.regions=cols_mask, scales=list(draw=F), cuts=length(cols_mask)-1, at=c(seq(min(mni152), max(mni152), length=length(cols_mask)-1)), colorkey = legend, main="", cex.main=2, useRaster=T) + levelplot(z_score[,,slice], xlab="", ylab="", col.regions=cols, cuts=length(cols)-1, scales=list(draw=F),at=c(seq(min, max, length=length(cols)-1)), useRaster=T) + layer(grid.text(plottitle, x=unit(0.03, "npc"), y=unit(0.95, "npc"), just=c("left", "top"), gp=gpar(fontsize=35, col="white"))))  
  }
  dev.off()
}

#--------------------
## example

# set palettes (red - blue color scale)
cols = rev(brewer.pal(n = 11, name = "RdBu"))
cols_mask = colorRampPalette(c("black", "white"))(40)
# color scale: virdis (used for creating empirical lesion rate images)
#library(viridisLite)
#cols <- viridis(100)

# Plot slice 45 for BLESS (parameter map & binary significance map).
for(i in 1:n_sim_v0){
file_name_beta = paste0('BLESS_', toString(i), '_beta_age_slice')
file_name_t = paste0('BLESS_', toString(i), '_t_age_slice')
brain_mask = readNIfTI(path_mask)
empir_prob = brain_mask
mni152 = readNIfTI(path_mni152)

load(paste0(path_results_BLESS, toString(i), "/params.RData"))

mask = as.vector(brain_mask)
ind0 = which(mask == 0, arr.ind = T)
ind1 = which(mask == 1, arr.ind = T)

beta = rep(NA, dimension1*dimension2*dimension3)
beta[ind1] = as.numeric(unlist(params$Beta[1,]))

beta_age = array(t, c(dimension1, dimension2, dimension3))
beta_age = as.nifti(beta_age)
beta_age = as.vector(beta_age@.Data)
beta_age[beta_age ==0] = NA
beta_age = as.nifti(array(beta_age, c(dim1,dim2,dim3)))

plottitle="BLESS" 
zscores_plot(slice=45, cols=cols, cols_mask = cols_mask, mni152 = mni152, img = beta_age, mask=brain_mask, 
             empir_prob=empir_prob, empir_threshold=0, name=file_name_beta, legend=T, min=-1, max=1)

t = rep(NA, dimension1*dimension2*dimension3)
t[ind1] = as.numeric(unlist(params$expected_gamma[1,]))

t_age = array(t, c(dimension1, dimension2, dimension3))
t_age = as.nifti(t_age)
t_age = as.vector(t_age@.Data)
t_age[t_age ==0] =NA
t_age = as.nifti(array(t_age, c(dim1,dim2,dim3)))
t_age[t_age > 0.5] = 1
t_age[t_age <= 0.5] = 0

plottitle="BLESS" 
zscores_plot(slice=45, cols=cols, cols_mask = cols_mask, mni152 = mni152, img = t_age, mask=brain_mask, 
             empir_prob=empir_prob, empir_threshold=0, name=file_name_t, legend=T, min=0, max=1)
}

# Plot slice 45 for Firth regression (parameter map & binary significance map).
file_name_beta = 'Firth_beta_corrected_slice'
file_name_t = 'Firth_t_corrected_slice'

brain_mask = readNIfTI(path_mask)
empir_prob = brain_mask
mni152 = readNIfTI(path_mni152)

params = list()
params$Beta = data.matrix(data.table::fread(paste0(path_results_firth, "test_statistic.csv"), header = T)[1,2:(length(ind1)+1)])
params$t = data.matrix(data.table::fread(paste0(path_results_firth, "Beta.csv"), header = T)[1,2:(length(ind1)+1)])

mask = as.vector(brain_mask)
ind0 = which(mask == 0, arr.ind = T)
ind1 = which(mask == 1, arr.ind = T)

beta = rep(NA, dimension1*dimension2*dimension3)
beta[ind1] = as.numeric(unlist(params$Beta[1,]))

beta_age = array(t, c(dimension1, dimension2, dimension3))
beta_age = as.nifti(beta_age)
beta_age = as.vector(beta_age@.Data)
beta_age[beta_age ==0] = NA
beta_age = as.nifti(array(beta_age, c(dim1,dim2,dim3)))

plottitle="Firth" 
zscores_plot(slice=45, cols=cols, cols_mask = cols_mask, mni152 = mni152, img = beta_age, mask=brain_mask, 
             empir_prob=empir_prob, empir_threshold=0, name=file_name_beta, legend=T, min=-1, max=1)

# FDR control
t = params$t
p = 2*pnorm(q=-abs(t))
t = p.adjust(p, method='fdr')
t[t<=0.05] = 1
t[t!=1] = 0

t = rep(NA, dimension1*dimension2*dimension3)
t[ind1] = as.numeric(unlist(t))

t_age = array(t, c(dimension1, dimension2, dimension3))
t_age = as.nifti(t_age)
t_age = as.vector(t_age@.Data)
t_age[t_age ==0] =NA
t_age = as.nifti(array(t_age, c(dim1,dim2,dim3)))

plottitle="Firth" 
zscores_plot(slice=45, cols=cols, cols_mask = cols_mask, mni152 = mni152, img = t_age, mask=brain_mask, 
             empir_prob=empir_prob, empir_threshold=0, name=file_name_t, legend=T, min=0, max=1)

########################################
### Scatterplot: Parameter Estimates ###
########################################

library(ggplot2)
library(latex2exp)
library(ggpointdensity)

# Function to plot a scatterplot between parameter values estimated via Firth regression and BLESS.
scatter_plot_comparison = function(params1, params2, model_names, title){
  
  df = cbind2(params1, params2)
  df = df[complete.cases(df),]
  df = data.frame(df)
  colnames(df) = c('params1', 'params2')
  
  params1_name = paste0(model_names[1])
  params2_name = paste0(model_names[2])
  
  x_min = y_min = min(c(min(params1), min(params2)))
  x_max = y_max = max(c(max(params1), max(params2)))
  
  ggplot(df, aes(x= params1, y=params2)) +
    #geom_point(shape=20,) + 
    geom_pointdensity(adjust = .1,show.legend = F) +
    #theme(plot.background=element_rect(fill = "#f7fafb", color=NA),
    #      panel.background = element_rect(fill = '#f7fafb'),
    #      legend.background = element_rect(fill = "#f7fafb", color = NA),
    #      panel.grid.major = element_blank(), 
    #      panel.grid.minor = element_blank()) + 
    theme(plot.background=element_rect(fill = "white", color=NA),
          panel.background = element_rect(fill = 'white'),
          legend.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    geom_abline(intercept = 0) + 
    xlim(-0.55,1.55) + 
    ylim(-0.55,1.55) +
    xlab(params1_name) + 
    ylab(params2_name) +
    theme(text = element_text(size = 35)) +
    ggtitle(title)+ guides(fill=guide_legend(title=""))
  
}

M = length(ind1)
P=4

t_Firth = read.csv(paste0(path_results_firth, "Beta.csv"), header = T)[,2:(M+1)]
t_Firth = as.numeric(unlist(t_Firth))
t_Firth = matrix(t_Firth, P, M)
beta_age_Firth = t_Firth[1,]

load(paste0(path_results_BLESS, "1/params.RData"))
beta_age_BLESS = params$Beta[1,]

slice = 45
name = 'FirthvsBLESS'
pdf(paste0(path_output, name, '.pdf'))
scatter_plot_comparison(beta_age_Firth, beta_age_BLESS, c("Firth","BLESS"), "")
dev.off()

#####################
### Marginal Plot ###
#####################

marginal_plot = function(marginal,v0){
  df = cbind(marginal, v0)
  df = data.frame(df)
  colnames(df) = c('marginal', 'v0')
  
  ggplot(df, aes(x= v0, y=marginal)) +
    geom_line()+
    geom_point(aes(x=df$v0[which.max(df$marginal)], y=df$marginal[which.max(df$marginal)]), colour="red", shape=19, size=5) +
    geom_point() + 
    theme_classic() +
    xlab('v0') + 
    ylab('marginal') +
    theme(text = element_text(size = 35),
          axis.text.y = element_blank(),
          plot.background=element_rect(fill = "#f7fafb", color=NA),
          panel.background = element_rect(fill = '#f7fafb'),
          legend.background = element_rect(fill = "#f7fafb", color = NA)) 
}

marginal = numeric(length(v0_list))
for(i in 1:length(v0)){
  marginal[i] = read.csv(paste0(path_results_BLESS, i, "/marginal.csv"))[,2]
}

marginal_plot(marginal, log(v0_list))

###########################
### Regularisation Plot ###
###########################

# Only plot voxels from slice 45 that are within the mask. Otherwise, the plot becomes too full 

regularisation_plot = function(result, voxels = 'some', border = F, truth = F, effect = 1){
  if(voxels == 'some'){
    p = 10
    idx = c(342,532,673,987,1109,2341,1745,1867,2034,1981)
  }
  if(voxels == 'all'){
    if(border == T){
      p = 2146
      idx = 1:2146
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
  
  plot(0,0,xlim = c((min(v0_list)),max(v0_list)), ylim = c((min(beta_list, na.rm = T) - 0.1), (max(beta_list, na.rm=T) + 0.1)), xlab =  TeX('$log(\\nu_0)$'), ylab = TeX('$\\beta$'), type = "n", cex.axis=1.15, cex.lab=1.15, , bty = "l")
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

mask = as.vector(brain_mask)
mask45 = as.vector(brain_mask[,,45])
ind0 = which(mask == 0, arr.ind = T)
ind1 = which(mask == 1, arr.ind = T)
ind145 = which(mask45 == 1, arr.ind = T)

load(paste0(path_results_BLESS, '1/params.RData'))
x = rep(NA, dimension1*dimension2*dimension3)
x[ind1] = as.numeric(unlist(params$expected_gamma[1,]))
x = array(x, c(dimension1, dimension2, dimension3))
x = as.vector(x[,,45])
x[x>0.5] = 1
x[x!=1] = 0
new_ind0 = which(x==0,arr.ind = T)
new_ind1 = which(x==1,arr.ind = T)

result = list()
result$v0_list = params$v0_list
result$gamma_list = matrix(NA, 5, 2147)
result$beta_list = matrix(NA, 5, 2147)
for(i in 1:5){
  load(paste0(path_results_BLESS, toString(i), '/params.RData'))
  x = rep(NA, dimension1*dimension2*dimension3)
  x[ind1] = as.numeric(unlist(params$expected_gamma[1,]))
  x = array(x, c(dimension1, dimension2, dimension3))
  x = as.vector(x[,,45])
  result$gamma_list[i,] = c(x[new_ind1], x[new_ind0])
  
  load(paste0(path_results_BLESS, toString(i), '/params.RData'))
  x = rep(NA, dimension1*dimension2*dimension3)
  x[ind1] = as.numeric(unlist(params$Beta[1,]))
  x = array(x, c(dimension1, dimension2, dimension3))
  x = as.vector(x[,,45])
  result$beta_list[i,] = c(x[new_ind1], x[new_ind0])
}

result$gamma_list = result$gamma_list[,-c(1008)]
result$beta_list = result$beta_list[,-c(1008)]

regularisation_plot(result, voxels = 'all', border = T, truth = F, effect = 1)

