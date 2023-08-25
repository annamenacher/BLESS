##############
### Script ###
##############

# File for plotting figures that contain comparison of posterior marginal distributions via KL-divergence & Wasserstein distance, 
# as well as comparison of posterior means and standard deviations, between BB-BLESS, BLESS-VI and BLESS-Gibbs
# (Code for Simulation Study)

# - random inactive voxel locations: 425, 1025
# - random active voxel locations: 1525, 2025

#################
### Libraries ###
#################

library(transport)
require(mvtnorm)
require(truncnorm)
require(brglm2)
require(LaplacesDemon)
library(ggplot2)
library(data.table)
library(MASS)
library(ggpointdensity)
library(latex2exp)

#################
### Constants ###
#################

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
# Hyperparameter: spike variance (plot marginal posterior of gamma list under prior v0=0 & identify plateau in marginal distribtuion /
# stablization of parameters)
# Index of best spike variance
v0_idx = 12

# Paths with folders of results.
# path: folder with true parameter values
path_truth = ''
# path: folder with BB-BLESS results
path_BB_BLESS = ''
# path: folder with BLESS-Gibbs results
path_Gibbs = ''
# path: folder with BLESS-VI results
path_BLESS_VI = ''
# path: output folder for plots
path_plots = ''

# Voxel locations that are not on the border of an image (due to setup of simulation study, they are not homogeneously distributed and should hence 
# be excluded from the evaluation)
idx = c(1:100, 2401:2500, 1151:1350, seq(1,2451,50), seq(2,2452,50), seq(49,2499,50), seq(50,2500,50), seq(24,2474,50), seq(25,2475,50), seq(26,2476,50), seq(27,2477,50))
 
# Data
X = matrix(cbind(c(rep(1, N/2), rep(0,N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

# Truth
load(paste0(path_truth, "/truth.RData"))
beta1_truth = truth$beta1
beta2_truth = truth$beta2

# BLESS-Gibbs
# Number of Gibbs samples: n_iter = 15,000 - burn_in = 5,000
n_iter = 10000
result_Gibbs_Probit = list()
result_Gibbs_Probit$beta_chain = data.table::fread(paste0(path_Gibbs,"beta1_chain.csv"))[,2:(M+1)]
result_Gibbs_Probit$beta_chain = data.matrix(result_Gibbs_Probit$beta_chain)
result_Gibbs_Probit$beta_chain = result_Gibbs_Probit$beta_chain[complete.cases(result_Gibbs_Probit$beta_chain),]
x = data.table::fread(paste0(path_Gibbs,"/gamma1_chain.csv"))[,2:(M+1)]
x = data.matrix(x)
result_Gibbs_Probit$gamma_chain = x[complete.cases(x),]

# BLESS VI
result_BLESS = list()
result_BLESS$Beta = matrix(0, P, M)
result_BLESS$Beta[1,] = as.numeric(unlist(data.table::fread(sprintf("%sBeta1%03d.csv", path_BLESS_VI, dataset))[v0_idx,2:(M+1)]))
result_BLESS$Beta[2,] = as.numeric(unlist(data.table::fread(sprintf("%sBeta2%03d.csv", path_BLESS_VI, dataset))[v0_idx,2:(M+1)]))
result_BLESS$expected_gamma = matrix(0, P, M)
result_BLESS$expected_gamma[1,] = as.numeric(unlist(data.table::fread(sprintf("%sgamma1%03d.csv", path_BLESS_VI, dataset))[v0_idx,2:(M+1)]))
result_BLESS$expected_gamma[2,] = as.numeric(unlist(data.table::fread(sprintf("%sgamma2%03d.csv", path_BLESS_VI, dataset))[v0_idx,2:(M+1)]))
result_BLESS$expected_gamma2 = matrix(0, P, M)
result_BLESS$expected_gamma2[1,] = as.numeric(unlist(data.table::fread(sprintf("%sgamma21%03d.csv", path_BLESS_VI, dataset))[v0_idx,2:(M+1)]))
result_BLESS$expected_gamma2[2,] = as.numeric(unlist(data.table::fread(sprintf("%sgamma22%03d.csv", path_BLESS_VI, dataset))[v0_idx,2:(M+1)]))
result_BLESS$theta = matrix(0, P, M)
result_BLESS$theta[1,] = as.numeric(unlist(data.table::fread(sprintf("%stheta1%03d.csv", path_BLESS_VI, dataset))[v0_idx,2:(M+1)]))
result_BLESS$theta[2,] = as.numeric(unlist(data.table::fread(sprintf("%stheta2%03d.csv", path_BLESS_VI, dataset))[v0_idx,2:(M+1)]))

mean_beta_emvs_bless = result_BLESS$Beta 

sd_beta_emvs_bless = matrix(NA, P*P, M)
for(j in 1:M){
sd_beta_emvs_bless[,j] = as.vector(solve(t(X)%*%X + diag(result_BLESS$expected_gamma2[,j], P)))
}

result_BLESS_sample_beta1 = matrix(NA, 1000, M)
result_BLESS_sample_beta2 = matrix(NA, 1000, M)

for(i in 1:1000){
  for(j in 1:M){
  samp = mvrnorm(1, mean_beta_emvs_bless[,j], matrix(sd_beta_emvs_bless[,j], P, P))
  result_BLESS_sample_beta1[i,j] = samp[1]
  result_BLESS_sample_beta2[i,j] = samp[2]
  } }

result_BLESS_sample_beta1 = result_BLESS_sample_beta1[complete.cases(result_BLESS_sample_beta1),] 
result_BLESS_sample_beta2 = result_BLESS_sample_beta2[complete.cases(result_BLESS_sample_beta2),] 
result_BLESS_sample = result_BLESS_sample_beta1

# BB-BLESS
result_BB_random = list()
x = data.table::fread(paste0(path_BB_BLESS,"beta1.csv"))[,2:(M+1)]
x = data.matrix(x)
result_BB_random$beta_chain = x[complete.cases(x), ]
 result_BB_random$beta_chain_mean = numeric(2500)
result_BB_random$beta_chain_sd = numeric(2500)
for(j in 1:M){
  x = result_BB_random$beta_chain[,j]
  result_BB_random$beta_chain_mean[j]=mean(x, na.rm = T)
  result_BB_random$beta_chain_sd[j]=sd(x, na.rm = T)
  }
x = data.table::fread(paste0(path_BB_BLESS,"gamm1.csv"))[,2:(M+1)]
x = x[complete.cases(x), ]
result_BB_random$gamma_chain = data.matrix(x)

print('data read in complete')

#################
### Histogram ###
#################

voxel = 425

options(repr.plot.width = 20, repr.plot.height = 8)
feature = 'n'
label_column = 'Legend'
idx_outlier = 1
a = data.frame(n=result_Gibbs_Probit$beta_chain[,voxel], Legend = rep('Gibbs', length(result_Gibbs_Probit$beta_chain[,voxel])))
b = data.frame(n=result_BLESS_sample_beta1[,voxel], Legend = rep('BLESS-VI', length(result_BLESS_sample_beta1[,voxel])))
c = data.frame(n=result_BB_random$beta_chain[,voxel], Legend = rep('BB-BLESS', length(result_BB_random$beta_chain[,voxel])))

df = do.call('rbind', list(a, b, c))
means = mean_df = c(mean(a$n), mean(b$n), mean(c$n))

if(voxel > 1251){
# effect
truth_voxel = (sum(beta1_truth[(dimension/2+3):(dimension-2),(dimension/2+3):(dimension-2)] + beta1_truth[(3):(dimension/2-2),(dimension/2+3):(dimension-2)])) / (((dimension/2 - 4)^2)*2)
}else{
# no effect
truth_voxel = (sum(beta1_truth[(dimension/2+3):(dimension-2),(3):(dimension/2-2)] + beta1_truth[(3):(dimension/2-2),(3):(dimension/2-2)])) / (((dimension/2 - 4)^2)*2)
}

pdf(paste0(path_plots, "BB_BLESS_hist_voxel", voxel, '.pdf'), width = 10, height = 10)

plt = ggplot(df, aes(x = eval(parse(text = feature)), fill = eval(parse(text = label_column)))) +
  
  geom_histogram(alpha = 0.5, position = "identity", aes(y = ..density..), color = "black") +
    geom_density(alpha = 0)+
    geom_vline(xintercept = means, color = c('#F8766D', '#00BA38', '#619CFF'), linetype = "solid", size = 1) +
    labs(x = 'Beta', y = "Density") +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = 'white'),
          legend.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = c(0.85, 0.87),
          text=element_text(size = 27.5)) + 
    ggtitle(paste0("Voxel ", voxel)) +
    xlab('Beta')
  plt = plt + guides(fill = guide_legend(title = label_column))
  plt = plt + geom_vline(data = data.frame(rbind(c(mean_beta_emvs_bless), c(sd_beta_emvs_bless))),aes(xintercept = truth_voxel),
               linetype = "solid", size = 1)
  plt
dev.off()

voxel = 1025

options(repr.plot.width = 20, repr.plot.height = 8)
feature = 'n'
label_column = 'Legend'
idx_outlier = 1
a = data.frame(n = result_Gibbs_Probit$beta_chain[,voxel], Legend = rep('Gibbs', length(result_Gibbs_Probit$beta_chain[,voxel])))
b = data.frame(n = result_BLESS_sample_beta1[,voxel], Legend = rep('BLESS-VI', length(result_BLESS_sample_beta1[,voxel])))
c = data.frame(n = result_BB_random$beta_chain[,voxel], Legend = rep('BB-BLESS', length(result_BB_random$beta_chain[,voxel])))

df = do.call('rbind', list(a, b, c))
means = mean_df = c(mean(a$n), mean(b$n), mean(c$n))

if(voxel > 1251){
  # effect
  truth_voxel = (sum(beta1_truth[(dimension/2+3):(dimension-2),(dimension/2+3):(dimension-2)] + beta1_truth[(3):(dimension/2-2),(dimension/2+3):(dimension-2)])) / (((dimension/2 - 4)^2)*2)
}else{
  # no effect
  truth_voxel = (sum(beta1_truth[(dimension/2+3):(dimension-2),(3):(dimension/2-2)] + beta1_truth[(3):(dimension/2-2),(3):(dimension/2-2)])) / (((dimension/2 - 4)^2)*2)
}

pdf(paste0(path_plots, "BB_BLESS_hist_voxel", voxel, '.pdf'), width = 10, height = 10)

plt = ggplot(df, aes(x = eval(parse(text = feature)), fill = eval(parse(text = label_column)))) +
  
  geom_histogram(alpha = 0.5, position = "identity", aes(y = ..density..), color = "black") +
  geom_density(alpha = 0)+
  geom_vline(xintercept = means, color = c('#F8766D', '#00BA38', '#619CFF'), linetype = "solid", size = 1) +
  labs(x = 'Beta', y = "Density") +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = 'white'),
        legend.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.87),
        text=element_text(size = 27.5)) + 
  ggtitle(paste0("Voxel ", voxel)) +
  xlab('Beta')
plt = plt + guides(fill = guide_legend(title = label_column))
plt = plt + geom_vline(data = data.frame(rbind(c(mean_beta_emvs_bless), c(sd_beta_emvs_bless))),aes(xintercept = truth_voxel),
                       linetype = "solid", size = 1)
plt
dev.off()

voxel = 1525

options(repr.plot.width = 20, repr.plot.height = 8)
feature = 'n'
label_column = 'Legend'
idx_outlier = 1
a = data.frame(n = result_Gibbs_Probit$beta_chain[,voxel], Legend = rep('Gibbs', length(result_Gibbs_Probit$beta_chain[,voxel])))
b = data.frame(n = result_BLESS_sample_beta1[,voxel], Legend = rep('BLESS-VI', length(result_BLESS_sample_beta1[,voxel])))
c = data.frame(n = result_BB_random$beta_chain[,voxel], Legend = rep('BB-BLESS', length(result_BB_random$beta_chain[,voxel])))

df = do.call('rbind', list(a, b, c))
means = mean_df = c(mean(a$n), mean(b$n), mean(c$n))

if(voxel > 1251){
  # effect
  truth_voxel = (sum(beta1_truth[(dimension/2+3):(dimension-2),(dimension/2+3):(dimension-2)] + beta1_truth[(3):(dimension/2-2),(dimension/2+3):(dimension-2)])) / (((dimension/2 - 4)^2)*2)
}else{
  # no effect
  truth_voxel = (sum(beta1_truth[(dimension/2+3):(dimension-2),(3):(dimension/2-2)] + beta1_truth[(3):(dimension/2-2),(3):(dimension/2-2)])) / (((dimension/2 - 4)^2)*2)
}

pdf(paste0(path_plots, "BB_BLESS_hist_voxel", voxel, '.pdf'), width = 10, height = 10)

plt = ggplot(df, aes(x = eval(parse(text = feature)), fill = eval(parse(text = label_column)))) +
  
  geom_histogram(alpha = 0.5, position = "identity", aes(y = ..density..), color = "black") +
  geom_density(alpha = 0)+
  geom_vline(xintercept = means, color = c('#F8766D', '#00BA38', '#619CFF'), linetype = "solid", size = 1) +
  labs(x = 'Beta', y = "Density") +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = 'white'),
        legend.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.87),
        text=element_text(size = 27.5)) + 
  ggtitle(paste0("Voxel ", voxel)) +
  xlab('Beta')
plt = plt + guides(fill = guide_legend(title = label_column))
plt = plt + geom_vline(data = data.frame(rbind(c(mean_beta_emvs_bless), c(sd_beta_emvs_bless))), aes(xintercept = truth_voxel),
                       linetype = "solid", size = 1)
plt
dev.off()

voxel = 2025

options(repr.plot.width = 20, repr.plot.height = 8)
feature = 'n'
label_column = 'Legend'
idx_outlier = 1
a = data.frame(n = result_Gibbs_Probit$beta_chain[,voxel], Legend = rep('Gibbs', length(result_Gibbs_Probit$beta_chain[,voxel])))
b = data.frame(n = result_BLESS_sample_beta1[,voxel], Legend = rep('BLESS-VI', length(result_BLESS_sample_beta1[,voxel])))
c = data.frame(n = result_BB_random$beta_chain[,voxel], Legend = rep('BB-BLESS', length(result_BB_random$beta_chain[,voxel])))

df = do.call('rbind', list(a, b, c))
means = mean_df = c(mean(a$n), mean(b$n), mean(c$n))

if(voxel > 1251){
  # effect
  truth_voxel = (sum(beta1_truth[(dimension/2+3):(dimension-2),(dimension/2+3):(dimension-2)] + beta1_truth[(3):(dimension/2-2),(dimension/2+3):(dimension-2)])) / (((dimension/2 - 4)^2)*2)
}else{
  # no effect
  truth_voxel = (sum(beta1_truth[(dimension/2+3):(dimension-2),(3):(dimension/2-2)] + beta1_truth[(3):(dimension/2-2),(3):(dimension/2-2)])) / (((dimension/2 - 4)^2)*2)
}

pdf(paste0(path_plots, "BB_BLESS_hist_voxel", voxel, '.pdf'), width = 10, height = 10)

plt = ggplot(df, aes(x = eval(parse(text = feature)), fill = eval(parse(text = label_column)))) +
  
  geom_histogram(alpha = 0.5, position = "identity", aes(y = ..density..), color = "black") +
  geom_density(alpha = 0)+
  geom_vline(xintercept = means, color = c('#F8766D', '#00BA38', '#619CFF'), linetype = "solid", size = 1) +
  labs(x = 'Beta', y = "Density") +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = 'white'),
        legend.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.87),
        text=element_text(size = 27.5)) + 
  ggtitle(paste0("Voxel ", voxel)) +
  xlab('Beta')
plt = plt + guides(fill = guide_legend(title = label_column))
plt = plt + geom_vline(data = data.frame(rbind(c(mean_beta_emvs_bless), c(sd_beta_emvs_bless))), aes(xintercept = truth_voxel),
                       linetype = "solid", size = 1)
plt
dev.off()

#####################
### KL-divergence ###
#####################

# Functions.
new_x = c()
for(i in 1:dimension){
  rep_x = rep(i, dimension)
  new_x = append(new_x, rep_x)
}
new_y = rep(rev(1:dimension), dimension)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# Plot 2D-image (jet-scale colors)
plot_image = function(Lesion_prob_matrix, limit_start = 0.0, limit_end = 0.1, legend.position = 'none'){

 # Input: 
 # - Lesion_prob_matrix: 2D image to plot
 # - limit_start: lower limit of values to plot
 # - limit_end: uppper limit of values to plot
 # - legend.position: legend position
 
  df <- reshape2::melt(Lesion_prob_matrix, varnames = c("x", "y"), value.name = "value")
  df$x = new_x
  df$y = new_y
  
  ggplot(df, aes_string(x = "x", y = "y", fill = "value")) + 
    geom_raster() +                        
    scale_fill_gradientn(colours = jet.colors(10),
                         limits = c(limit_start, limit_end),
                         guide = "colourbar",
                         oob = squish) +
    theme_void() + 
    theme(legend.position = legend.position)
}

# Plot 2D-image (black-white colors for binary values: black 0 and white 1)
plot_binary_image = function(empirical_matrix){
 
 # Input: 
 # - empirical_matrix: 2D binary image to plot
 
  df <- reshape2::melt(empirical_matrix, varnames = c("x", "y"), value.name = "value")
  df$x = new_x
  df$y = new_y
  df$value = as.factor(df$value)
  ggplot(df, aes(x = x, y = y)) + 
    geom_raster(aes(fill = factor(value))) +
    scale_fill_manual(values=c("0" = "black", "1" = "white"), guide = "none") +
    theme_void()
}

# Plot 2D-image (red - blue scale colors)
plot_image_rdbu = function(Lesion_prob_matrix, limit_start = 0.0, limit_end = 0.1, legend.position = 'none'){

 # Input: 
 # - Lesion_prob_matrix: 2D image to plot
 # - limit_start: lower limit of values to plot
 # - limit_end: uppper limit of values to plot
 # - legend.position: legend position
 
  df <- reshape2::melt(Lesion_prob_matrix, varnames = c("x", "y"), value.name = "value")
  df$x = new_x
  df$y = new_y
  
  cols = rev(brewer.pal(n = 11, name = "RdBu"))
  
  ggplot(df, aes_string(x = "x", y = "y", fill = "value")) + 
    geom_raster() +                        
    scale_fill_gradientn(colours = cols,
                         limits = c(limit_start, limit_end),
                         guide = "colourbar",
                         oob = squish) +
    theme_void() + 
    theme(legend.position = legend.position)
}

# Plot 2D-image (viridis color scale)
plot_image_viridis = function(Lesion_prob_matrix, limit_start = 0.0, limit_end = 0.1, legend.position = 'none'){
  df <- reshape2::melt(Lesion_prob_matrix, varnames = c("x", "y"), value.name = "value")
  df$x = new_x
  df$y = new_y
  
  ggplot(df, aes_string(x = "x", y = "y", fill = "value")) + 
    geom_raster() +                        
    scale_fill_viridis_c(limits = c(limit_start, limit_end), guide = 'colourbar') +
    theme_void() + 
    theme(legend.position = legend.position)
}
 
# Plot scatterplot comparing two sets of values params1 vs params2.
scatter_plot_comparison = function(params1, params2, model_names, title){
  
  df = cbind2(params1, params2)
  df = df[complete.cases(df),]
  df = data.frame(df)
  colnames(df) = c('params1', 'params2')

  params1_name = paste0(model_names[1])
  params2_name = paste0(model_names[2])
  
  x_min = y_min = min(c(min(params1), min(params2)))
  x_max = y_max = max(c(max(params1), max(params2)))
  
  ggplot(df, aes(x = params1, y = params2)) +
    geom_point(shape = 20) + 
    theme_classic() + 
    geom_abline(intercept = 0) + 
    xlim(x_min, x_max) + 
    ylim(y_min, y_max) +
    xlab(params1_name) + 
    ylab(params2_name) +
    theme(text = element_text(size = 18)) +
    ggtitle(title)
}

# Plot for creating boxplot comparing set of values params1 vs params2 vs params3 vs params4.
boxplot_comparison = function(params1, params2, params3, params4, model_names, title){
  
  model = c(rep(model_names[2], length(c(1:2500)[-idx])), rep(model_names[3], length(c(1:2500)[-idx])), rep(model_names[4], length(c(1:2500)[-idx])))
  params1 = rep(params1, 3)
  params2 = c(params2, params3, params4)
  df = cbind2(params1, params2)
  df = df[complete.cases(df),]
  df = data.frame(df)
  colnames(df) = c('Effect' , 'Beta')
  df$Effect = round(c(df$Effect), 3)
  df$Effect = as.factor(df$Effect)
  df$Model = model
  
  ggplot(df,aes(x = Model, y = Beta)) +
    geom_boxplot(aes(fill = Effect)) +
    theme_minimal() + 
    labs(x = 'Model', y = 'Beta') +
    ylim(-0.25, 1.25) +
    ggtitle(title)
}

# Logistic function
logistic = function(x){
  x = exp(x)/(1+exp(x))
  return(x)
}

# Truth
# Beta1
truth_image = matrix(0, dimension, dimension)
truth_image[1:(dimension), (dimension/2+1):(dimension)] = (sum(beta1_truth[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)]) + sum(beta1_truth[(dimension/2 + 3):(dimension-2),(dimension/2+3):(dimension-2)])) / (((dimension/2 - 4)^2)*2)
truth_image[1:(dimension), 1:(dimension/2)] = (sum(beta1_truth[3:(dimension/2-2),3:(dimension/2-2)]) + sum(beta1_truth[(dimension/2 + 3):(dimension-2),3:(dimension/2-2)])) / (((dimension/2 - 4)^2)*2)

# Beta2 (commented out here but in case the evaluation of beta2 was of interest -> comment out beta1 and use the truth_image defined below)
#truth_image = matrix(0, dimension, dimension)
#truth_image[1:(dimension/2), 1:(dimension/2)] = (sum(beta2_truth[3:(dimension/2-2),3:(dimension/2-2)]) + sum(beta2_truth[3:(dimension/2-2),(dimension/2+3):(dimension-2)])+ sum(beta2_truth[(dimension/2 + 3):(dimension/2-2),(dimension/2+3):(dimension-2)])) / (((dimension/2 - 4)^2)*3)
#truth_image[1:(dimension), (dimension/2+1):(dimension)] = (sum(beta2_truth[3:(dimension/2-2),3:(dimension/2-2)]) + sum(beta2_truth[3:(dimension/2-2),(dimension/2+3):(dimension-2)])+ sum(beta2_truth[(dimension/2 + 3):(dimension/2-2),(dimension/2+3):(dimension-2)])) / (((dimension/2 - 4)^2)*3)
#truth_image[(dimension/2+1):(dimension), 1:(dimension/2)] = (sum(beta2_truth[(dimension/2+3):(dimension-2),(3):(dimension/2-2)]) ) / (((dimension/2 - 4)^2)*1)

# BB-BLESS
mean_bb_bless = numeric(M)
sd_bb_bless = numeric(M)
for(j in 1:M){
  x=result_BB_random$beta_chain[,j]
    mean_bb_bless[j] = mean(x)
    sd_bb_bless[j] = sd(x)
}

### Mean ###
setwd(path_plots)
beta_idx = 1
# Comparison of Truth vs BB-BLESS / BLESS-VI / Gibbs (Boxplots)
pdf("boxplot_mean.pdf", width = 15, height = 10) 
boxplot_comparison(as.vector(truth_image)[-idx], apply(result_Gibbs_Probit$beta_chain[,],2,mean)[-idx], apply(result_BLESS_sample, 2, mean)[-idx], mean_bb_bless[-idx], c('Truth', 'Gibbs', 'BLESS-VI', 'BB-BLESS'), paste0("Mean (Beta", beta_idx, "): Truth vs BB-BLESS / BLESS-VI / Gibbs"))
dev.off()

# Comparison of Gibbs vs BB-BLESS
pdf("scatterplot_mean_Gibbs_BB_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(apply(result_Gibbs_Probit$beta_chain, 2, mean)[-idx], mean_bb_bless[-idx], c('Gibbs', 'BB-BLESS'), paste0('Mean: Beta', beta_idx))
dev.off()
# Comparison of Gibbs vs BLESS-VI
pdf("scatterplot_mean_Gibbs_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(apply(result_Gibbs_Probit$beta_chain, 2, mean)[-idx], apply(result_BLESS_sample, 2, mean)[-idx], c('Gibbs', 'BLESS-VI'), paste0('Mean: Beta', beta_idx))
dev.off()
# Comparison of BLESS-VI vs BB-BLESS
pdf("scatterplot_mean_BB_BLESS_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(apply(result_BLESS_sample, 2, mean)[-idx], mean_bb_bless[-idx], c('BLESS-VI', 'BB-BLESS'), paste0('Mean: Beta', beta_idx))
dev.off()

# Standard deviation

# Comparison of Gibbs vs BB-BLESS
pdf("scatterplot_sd_Gibbs_BB_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(apply(result_Gibbs_Probit$beta_chain, 2, sd)[-idx], sd_bb_bless[-idx], c('Gibbs', 'BB-BLESS'), paste0('Standard Deviation: Beta', beta_idx))
dev.off()
# Comparison of Gibbs vs BLESS-VI
pdf("scatterplot_sd_Gibbs_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(apply(result_Gibbs_Probit$beta_chain, 2, sd)[-idx], apply(result_BLESS_sample, 2, sd)[-idx], c('Gibbs', 'BLESS-VI'), paste0('Standard Deviation: Beta', beta_idx))
dev.off()
# Comparison of BLESS-VI vs BB-BLESS
pdf("scatterplot_sd_BB_BLESS_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(apply(result_BLESS_sample, 2, sd)[-idx], sd_bb_bless[-idx], c('BLESS-VI', 'BB-BLESS'), paste0('Standard Deviation: Beta', beta_idx))
dev.off()

### Bias ###
bias_gibbs = apply(result_Gibbs_Probit$beta_chain, 2, mean) - as.vector(truth_image)
bias_bless_vi = apply(result_BLESS_sample, 2, mean) - as.vector(truth_image)
bias_bb_bless =  mean_bb_bless - as.vector(truth_image)

pdf("scatterplot_bias_Gibbs_BB_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(bias_gibbs[-idx], bias_bb_bless[-idx], c('Gibbs', 'BB-BLESS'), paste0('Bias: Beta', beta_idx))
dev.off()
pdf("scatterplot_bias_Gibbs_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(bias_gibbs[-idx], bias_bless_vi[-idx], c('Gibbs', 'BLESS-VI'), paste0('Bias: Beta', beta_idx))
dev.off()
pdf("scatterplot_bias_BB_BLESS_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(bias_bless_vi[-idx], bias_bb_bless[-idx], c('BLESS-VI', 'BB-BLESS'), paste0('Bias: Beta', beta_idx))
dev.off()

### MSE ###
mse_gibbs = apply(result_Gibbs_Probit$beta_chain, 2, sd)^2 + bias_gibbs^2
mse_bless_vi = apply(result_BLESS_sample, 2, sd)^2 + bias_bless_vi^2
mse_bb_bless =  sd_bb_bless^2 + bias_bb_bless^2

pdf("scatterplot_mse_Gibbs_BB_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(mse_gibbs[-idx], mse_bb_bless[-idx], c('Gibbs', 'BB-BLESS'), paste0('MSE: Beta', beta_idx))
dev.off()
pdf("scatterplot_mse_Gibbs_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(mse_gibbs[-idx], mse_bless_vi[-idx], c('Gibbs', 'BLESS-VI'), paste0('MSE: Beta', beta_idx))
dev.off()
pdf("scatterplot_mse_BB_BLESS_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(mse_bless_vi[-idx], mse_bb_bless[-idx], c('BLESS-VI', 'BB-BLESS'), paste0('MSE: Beta', beta_idx))
dev.off()

### Probability of Inclusion ###
# Effect Truth Image (map of effect)

# Beta1
effect_truth_image = matrix(0, dimension, dimension)
effect_truth_image[1:(dimension), 1:(dimension/2)] = 0
effect_truth_image[1:dimension, (dimension/2+1):(dimension)] = 1

# Beta2 (commented out here but in case the evaluation of beta2 was of interest -> comment out beta1 and use the effect_truth_image defined below)
#effect_truth_image = matrix(0, dimension, dimension)
#effect_truth_image[1:(dimension/2), 1:(dimension/2)] = 0
#effect_truth_image[1:dimension, (dimension/2+1):(dimension)] = 1
#effect_truth_image[(dimension/2+1):dimension, 1:(dimension/2)] = 0


### Test statistic: t = beta / sd(beta) ###
t_gibbs = apply(result_Gibbs_Probit$beta_chain, 2, mean) / apply(result_Gibbs_Probit$beta_chain, 2, sd)
t_bless_vi = apply(result_BLESS_sample, 2, mean) / apply(result_BLESS_sample, 2, sd)
t_bb_bless = mean_bb_bless / sd_bb_bless

# Comparison of Gibbs vs BB-BLESS
pdf("scatterplot_test_stat_Gibbs_BB_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(t_gibbs[-idx], t_bb_bless[-idx], c('Gibbs', 'BB-BLESS'), paste0('Test Statistic: Beta', beta_idx))
dev.off()
# Comparison of Gibbs vs BLESS-VI
pdf("scatterplot_test_stat_Gibbs_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(t_gibbs[-idx], t_bless_vi[-idx], c('Gibbs', 'BLESS-VI'), paste0('Test Statistic: Beta', beta_idx))
dev.off()
# Comparison of BLESS-VI vs BB-BLESS
pdf("scatterplot_test_stat_BB_BLESS_BLESS.pdf", width = 10, height = 10) 
scatter_plot_comparison(t_bless_vi[-idx], t_bb_bless[-idx], c('BLESS-VI', 'BB-BLESS'), paste0('Test Statistic: Beta', beta_idx))
dev.off()

### Test statistic:Binary results ### 
 
t_gibbs_binary = t_gibbs
t_gibbs_binary[abs(t_gibbs) > 1.96] = 1
t_gibbs_binary[t_gibbs_binary != 1] = 0
t_bless_vi = t_bless_vi_binary = result_BLESS$expected_gamma[1,]
t_bless_vi_binary[t_bless_vi > 0.5] = 1
t_bless_vi_binary[t_bless_vi <= 0.5] = 0
t_bb_bless_binary = t_bb_bless
t_bb_bless_binary[abs(t_bb_bless) > 1.96] = 1
t_bb_bless_binary[t_bb_bless_binary != 1] = 0

# Comparison of Gibbs vs BB-BLESS
pdf("gibbs_binary_significance_map.pdf", width = 10, height = 10) 
plot_binary_image(matrix(t_gibbs_binary, 50, 50))
dev.off()
# Comparison of Gibbs vs BLESS-VI
pdf("bless_vi_binary_significance_map.pdf", width = 10, height = 10) 
plot_binary_image(matrix(t_bless_vi_binary, 50, 50))
dev.off()
# Comparison of BLESS-VI vs BB-BLESS
pdf("bless_bb_binary_significance_map.pdf", width = 10, height = 10) 
plot_binary_image(matrix(t_bb_bless_binary, 50, 50))
dev.off()

### Standard Deviation ###
params1 = apply(result_Gibbs_Probit$beta_chain, 2, sd)[-idx]
params2 = sd_bb_bless[-idx]

pdf(paste0(path_plots, 'sd_bb_vs_gibbs.pdf'), width = 10, height = 10)
plot(params1[abs(t_bb_bless[-idx]) > 1.96], params2[abs(t_bb_bless[-idx]) > 1.96], col = 'red', xlim = c(0.0, 0.5), ylim = c(0.0, 0.5), xlab = 'Gibbs', ylab = 'BB-BLESS', main = 'Standard Deviation: Beta1')
points(params1[abs(t_bb_bless[-idx]) < 1.96], params2[abs(t_bb_bless[-idx]) < 1.96], col = 'blue')
abline(0, 1)
dev.off()

# Density Estimation
px = matrix(NA, 512, M)
py = matrix(NA, 512, M)
pz = matrix(NA, 512, M)

for(j in 1:M){
voxel_string = paste0("X",j)
p = ggplot() +
  geom_density(data = data.frame(result_BLESS_sample_beta1), aes_string(voxel_string))
p = ggplot_build(p)
pz[,j] = p$data[[1]]$y

voxel_string = paste0("V",j)
p = ggplot() +
  geom_density(data=data.frame(result_Gibbs_Probit$beta_chain), aes_string(voxel_string))
p = ggplot_build(p)
px[,j] = p$data[[1]]$y

x=result_BB_random$beta_chain[,j]
idx_outlier = which((abs(x - median(x)) / mad(x)) > 5, arr.ind = T)
if(length(idx_outlier) > 0){
q = ggplot() +
   geom_density(data = data.frame(result_BB_random$beta_chain[-idx_outlier,]), aes_string(voxel_string))
q = ggplot_build(q)
py[,j] = q$data[[1]]$y
} else{
  q = ggplot() +
    geom_density(data = data.frame(result_BB_random$beta_chain), aes_string(voxel_string))
  q = ggplot_build(q)
  py[,j] = q$data[[1]]$y
}
}

# KL-divergence
KL_divergence = numeric(M)
for(j in 1:M){
  KL_divergence[j] = LaplacesDemon::KLD(px[,j], py[,j])$sum.KLD.px.py
}

pdf(paste0(path_plots, 'KLD_gibbs_vs_bb_bless.pdf'), width = 10, height = 10)
print(ggplot(data = data.frame(c(1:2500)[-idx], KL_divergence[-idx]), mapping = aes(x = c(1:2500)[-idx], y = KL_divergence[-idx])) +
        geom_pointdensity(adjust = .1) + xlab("Voxel") + ylab("KL-divergence") + labs(color="Count") + theme_light())
dev.off()

# Wasserstein distance
wasserstein_metric = numeric(M)
for(j in 1:M){
  wasserstein_metric[j] = transport::wasserstein1d(px[,j], py[,j])
}

pdf(paste0(path_plots, 'wasserstein_gibbs_vs_bb.pdf'), width = 10, height = 10)
print(ggplot(data = data.frame(c(1:2500)[-idx], wasserstein_metric[-idx]), mapping = aes(x = c(1:2500)[-idx], y = wasserstein_metric[-idx])) +
        geom_pointdensity(adjust = .1) + xlab("Voxel") + ylab("Wasserstein metric") + labs(color = "Count") + theme_light())
dev.off()

kl_divergence_bb = KL_divergence
wasserstein_metric_bb = wasserstein_metric

# Gibbs vs BLESS VI
# KL-divergence
KL_divergence = numeric(M)
for(j in 1:M){
  KL_divergence[j] = LaplacesDemon::KLD(px[,j], pz[,j])$sum.KLD.px.py
}

pdf(paste0(path_plots, 'KLD_gibbs_vs_bless_vi.pdf'), width = 10, height = 10)
print(ggplot(data = data.frame(c(1:2500)[-idx], KL_divergence[-idx]), mapping = aes(x = c(1:2500)[-idx], y = KL_divergence[-idx])) +
        geom_pointdensity(adjust = .1) + xlab("Voxel") + ylab("KL-divergence") + labs(color = "Count") + theme_light())
dev.off()

# Wasserstein distance
wasserstein_metric = numeric(M)
for(j in 1:M){
  wasserstein_metric[j] = transport::wasserstein1d(px[,j], pz[,j])
}

pdf(paste0(path_plots, 'wasserstein_gibbs_vs_bless_vi.pdf'), width = 10, height = 10) 
print(ggplot(data = data.frame(c(1:2500)[-idx], wasserstein_metric[-idx]), mapping = aes(x = c(1:2500)[-idx], y = wasserstein_metric[-idx])) +
        geom_pointdensity(adjust = .1) + xlab("Voxel") + ylab("Wasserstein metric") + labs(color = "Count") + theme_light())

kl_divergence_vi = KL_divergence
wasserstein_metric_vi = wasserstein_metric


options(repr.plot.width = 20, repr.plot.height = 8)
metric = 'Wasserstein'
  feature = 'n'
  label_column = 'Legend'
  if(metric == 'Wasserstein'){
  pdf(paste0(path_plots, 'wasserstein_histogram_bb_vs_vi.pdf'), width = 10, height = 10)
  a = data.frame(n = wasserstein_metric_vi, Legend = rep('BLESS-VI', length(wasserstein_metric_vi)))
b = data.frame(n = wasserstein_metric_bb, Legend = rep('BB-BLESS', length(wasserstein_metric_bb)))
  } else{
    pdf(paste0(path_plots, 'kld_bb_vs_vi.pdf'), width = 10, height = 10)
  a = data.frame(n = kl_divergence_vi, Legend = rep('BLESS-VI', length(kl_divergence_vi)))
b = data.frame(n = kl_divergence_bb, Legend = rep('BB-BLESS', length(kl_divergence_bb)))
  }
  
df = do.call('rbind', list(a,b))

  plt = ggplot(df, aes(x = eval(parse(text = feature)), fill = eval(parse(text = label_column)))) +
    geom_histogram(alpha = 0.5, position = "identity", aes(y = ..density..), color = "black") +
    labs(x = 'Wasserstein distance', y = "Count") +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = 'white'),
          legend.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = c(0.85, 0.87),
          text=element_text(size = 27.5)) +
    ggtitle('Comparison of Gibbs vs BLESS') 
  plt = plt + guides(fill=guide_legend(title = label_column))
  plt
  dev.off()

metric = 'KL-divergence'
 feature = 'n'
  label_column = 'Legend'
  if(metric == 'Wasserstein'){
  pdf(paste0(path_plots, 'wasserstein_histogram_bb_vs_vi.pdf'), width = 10, height = 10)
  a = data.frame(n = wasserstein_metric_vi, Legend = rep('BLESS-VI', length(wasserstein_metric_vi)))
b = data.frame(n = wasserstein_metric_bb, Legend = rep('BB-BLESS', length(wasserstein_metric_bb)))
  } else{
    pdf(paste0(path_plots, 'kld_bb_vs_vi.pdf'), width = 10, height = 10)
  a = data.frame(n = kl_divergence_vi, Legend = rep('BLESS-VI', length(kl_divergence_vi)))
b = data.frame(n = kl_divergence_bb, Legend = rep('BB-BLESS', length(kl_divergence_bb)))
  }

df = do.call('rbind', list(a, b))

  plt = ggplot(df, aes(x = eval(parse(text = feature)), fill = eval(parse(text = label_column)))) +
    geom_histogram(alpha = 0.5, position = "identity", aes(y = ..density..), color = "black") +
    labs(x = 'KL-divergence', y = "Count") +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = 'white'),
          legend.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.85, 0.87),
          text=element_text(size = 27.5)) +
    ggtitle('Comparison of Gibbs vs BLESS')
  plt = plt + guides(fill=guide_legend(title = label_column))
  plt
  dev.off()





