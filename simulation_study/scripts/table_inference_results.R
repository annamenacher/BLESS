##############
### Script ###
##############

# Script to generate tables containing inference results.

#################
### Libraries ###
#################

library(ggplot2)
library(viridisLite)
library(RColorBrewer)

library(scales)
library(dplyr)
library(reshape2)

library(abind)
library(xtable)

#################
### Constants ###
#################

n_sim = 100
M = 2500
dim1 = dim2 = dimension = sqrt(M)
P = 2

lambda_list = c(1, 2, 3)
N_list = c(500, 1000, 5000)

v0 = exp(seq(-20, -1, length.out = 15))

# Path containing the true values.
path_truth = ''
# Path contain BLESS-VI results.
path_results_BLESS = ''
# Path contain BSGLMM results.
path_results_BSGLMM = ''
# Path contain Firth results.
path_results_Firth = ''
# Path contain BB-BLESS results.
path_results_BB_BLESS = ''
# Path contain BLESS-Gibbs results.
path_results_Gibbs = ''
# Path to output table in.
path_logs = ''

#################
### Functions ###
#################

plot_binary_image = function(empirical_matrix){
  # For plotting images. 
  new_x = c()
  for(i in 1:dimension){
    rep_x = rep(i, dimension)
    new_x = append(new_x, rep_x)
  }
  new_y = rep(rev(1:dimension), dimension)
  
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  df <- reshape2::melt(empirical_matrix, varnames = c("x", "y"), value.name = "value")
  df$x = new_x
  df$y = new_y
  df$value = as.factor(df$value)
  ggplot(df, aes(x = x, y = y)) + 
    geom_raster(aes(fill=factor(value))) +
    scale_fill_manual(values=c("0" = "black", "1" = "white"), guide = "none") +
    theme_void()
}

plot_image_rdbu = function(Lesion_prob_matrix, limit_start = 0.0, limit_end = 0.1, legend.position = 'none'){
  new_x = c()
  for(i in 1:dimension){
    rep_x = rep(i, dimension)
    new_x = append(new_x, rep_x)
  }
  new_y = rep(rev(1:dimension), dimension)
  
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
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

calculate_bias = function(truth, prediction, n_sim = 1, effect_size = 0.5){
  
  if(effect_size == 0.5){
    bias = list()
    
    truth_avg = numeric(2)
    truth_avg[1] = (sum(truth[3:(dimension/2 - 2), 3:(dimension/2 - 2)])+sum(truth[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)])) / (((dimension/2 - 4)^2)*2)
    truth_avg[2] = (sum(truth[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])+sum(truth[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])) / (((dimension/2 - 4)^2)*2)
    
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        result[i,1] = mean(abind(pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],pred[(dimension/2+3):(dimension - 2), 3:(dimension/2 - 2)], along = 2), na.rm = T) - truth_avg[1]
        result[i,2] = mean(abind(pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)], along=2), na.rm = T) - truth_avg[2]
        result[i,3] = result[i,1] + result[i,2]
      }
    }
    
    bias$left = mean(result[,1], na.rm = T)
    bias$right = mean(result[,2], na.rm = T)
    bias$total = mean(result[,3], na.rm = T)
    
    return(bias)
  }
  
  if(effect_size == 0.25){
    bias = list()
    
    truth_avg = numeric(2)
    truth_avg[1] = (sum(truth[3:(dimension/2 - 2), 3:(dimension/2 - 2)]) + sum(truth[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)])+ sum(truth[(dimension/2 + 3):(dimension-2),(dimension/2 + 3):(dimension-2)])) / (((dimension/2 - 4)^2)*3)
    truth_avg[2] = sum(truth[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))]) / (((dimension/2 - 4)^2)*1)
    
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        result[i,1] = mean(abind(pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],pred[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)],pred[(dimension/2 + 3):(dimension-2),(dimension/2 + 3):(dimension-2)],along=2),na.rm = T)-truth_avg[1]
        result[i,2] = mean(abind(pred[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))],along=2), na.rm = T)  -truth_avg[2]
        result[i,3] = result[i,1] + result[i,2]
      }
    }
    
    bias$left = mean(result[,1], na.rm = T)
    bias$right = mean(result[,2], na.rm = T)
    bias$total = mean(result[,3], na.rm = T)
    
    return(bias)
    
  }
  
  if(effect_size == 0.75){
    bias = list()
    
    truth_avg = numeric(2)
    truth_avg[1] = (sum(truth[3:(dimension/2 - 2), 3:(dimension/2 - 2)])) / (((dimension/2 - 4)^2)*1)
    truth_avg[2] = (sum(truth[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])+sum(truth[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)]) +sum(truth[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)])) / (((dimension/2 - 4)^2)*3)
    
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        result[i,1] = mean(abind(pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],along=2), na.rm = T) - truth_avg[1]
        result[i,2] = mean(abind(pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2+3):(dimension - 2), 3:(dimension/2 - 2)],along=2),na.rm = T) - truth_avg[2]
        result[i,3] = result[i,1] + result[i,2]
      }
    }
    
    bias$left = mean(result[,1], na.rm = T)
    bias$right = mean(result[,2], na.rm = T)
    bias$total = mean(result[,3], na.rm = T)
    
    return(bias)
  }
  
  if(effect_size == 1.0){
    bias = list()
    
    truth_avg = (sum(truth[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])+sum(truth[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)]) +sum(truth[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)])+sum(truth[3:(dimension/2 - 2), 3:(dimension/2 - 2)])) / (((dimension/2 - 4)^2)*4)
    
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        result[i,1] = 0
        result[i,2] = mean(abind(pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2+3):(dimension - 2), 3:(dimension/2 - 2)],pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],along=2),na.rm = T) - truth_avg
        result[i,3] = result[i,1] + result[i,2]
      }
    }
    bias$left = mean(result[,1], na.rm = T)
    bias$right = mean(result[,2], na.rm = T)
    bias$total = mean(result[,3], na.rm = T)
    return(bias)
  }
  
}

calculate_variance = function(truth, prediction, n_sim = 1, effect_size = 0.5){
  
  if(effect_size == 0.5){
    variance = list()
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        
        est_avg = numeric(2)
        est_avg[1] = mean(abind(pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],pred[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)],along=2), na.rm = T)
        est_avg[2] = mean(abind(pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)],along=2), na.rm = T)
        
        result[i,1] = mean(abind((est_avg[1] - pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2 , (est_avg[1] - pred[(dimension/2 +3):(dimension - 2), 3:(dimension/2 - 2)])^2,along=2), na.rm=T)
        result[i,2] = mean(abind((est_avg[2] - pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])^2,(est_avg[2] - pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])^2,along=2), na.rm=T)
        result[i,3] = mean(abind((est_avg[1] - pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2,(est_avg[1] - pred[(dimension/2 +3):(dimension - 2), 3:(dimension/2 - 2)])^2,
                                 (est_avg[2] - pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])^2,(est_avg[2] - pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])^2,along=2), na.rm = T)
      }
    }
    
    variance$left = mean(result[,1], na.rm = T)
    variance$right = mean(result[,2], na.rm = T)
    variance$total = mean(result[,3], na.rm = T)
    
  }
  
  if(effect_size == 0.25){
    variance = list()
    
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        
        est_avg = numeric(2)
        est_avg[1] = mean(abind(pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)] , pred[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)], pred[(dimension/2 +3):(dimension-2),(dimension/2 + 3):(dimension-2)],along=2), na.rm = T)
        est_avg[2] = mean(abind(pred[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))],along=2), na.rm = T) 
        
        result[i,1] = mean(abind((est_avg[1] - pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2 , (est_avg[1] - pred[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)])^2, (est_avg[1] - pred[(dimension/2+3):(dimension-2),(dimension/2 + 3):(dimension-2)])^2,along=2), na.rm=T)
        result[i,2] = mean(abind((est_avg[2] - pred[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))])^2,along=2), na.rm = T) 
        result[i,3] = mean(abind((est_avg[1] - pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2 ,
                                 (est_avg[1] - pred[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)])^2 ,
                                 (est_avg[1] - pred[(dimension/2+3):(dimension-2),(dimension/2 + 3):(dimension-2)])^2 ,
                                 (est_avg[2] - pred[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))])^2,along=2), na.rm = T)
      }
    }
    
    variance$left = mean(result[,1], na.rm = T)
    variance$right = mean(result[,2], na.rm = T)
    variance$total = mean(result[,3], na.rm = T)
    
  }
  
  if(effect_size == 0.75){
    variance = list()
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        
        est_avg = numeric(2)
        est_avg[1] = mean(abind(pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],along=2), na.rm = T )
        est_avg[2] = mean(abind(pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)] ,pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)],along=2), na.rm = T)
        
        result[i,1] = mean(abind((est_avg[1] - pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2,along=2), na.rm = T)
        result[i,2] = mean(abind((est_avg[2] - pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])^2,(est_avg[2] - pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])^2,(est_avg[2] - pred[(dimension/2 +3):(dimension - 2), 3:(dimension/2 - 2)])^2,along=2), na.rm = T)
        result[i,3] = mean(abind((est_avg[1] - pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2,(est_avg[2] - pred[(dimension/2 +3):(dimension - 2), 3:(dimension/2 - 2)])^2, 
                                 (est_avg[2] - pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])^2,(est_avg[2] - pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])^2,along=2), na.rm = T)
      }
    }
    
    variance$left = mean(result[,1], na.rm = T)
    variance$right = mean(result[,2], na.rm = T)
    variance$total = mean(result[,3], na.rm = T)
    
  }
  
  if(effect_size == 1.0){
    variance = list()
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        
        est_avg = mean(abind(pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)] , pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)], pred[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)],pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],along=2), na.rm = T)
        
        result[i,1] = 0
        result[i,2] = mean(abind((est_avg - pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])^2,(est_avg - pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])^2,(est_avg - pred[(dimension/2 +3):(dimension - 2), 3:(dimension/2 - 2)])^2 ,(est_avg - pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2,along=2), na.rm = T)
        result[i,3] = result[i,2]
      }
    }
    variance$left = mean(result[,1], na.rm = T)
    variance$right = mean(result[,2], na.rm = T)
    variance$total = mean(result[,3], na.rm = T)
    
  }
  return(variance)
}

calculate_mse = function(truth, prediction, n_sim = 1, effect_size = 0.5){
  
  if(effect_size == 0.5){
    mse = list()
    
    truth_avg = numeric(2)
    truth_avg[1] = (sum(truth[3:(dimension/2 - 2), 3:(dimension/2 - 2)])+sum(truth[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)])) / (((dimension/2 - 4)^2)*2)
    truth_avg[2] = (sum(truth[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])+sum(truth[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])) / (((dimension/2 - 4)^2)*2)
    
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        result[i,1] = mean((abind(pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],pred[(dimension/2+3):(dimension - 2), 3:(dimension/2 - 2)], along = 2) - truth_avg[1])^2, na.rm = T)
        result[i,2] = mean((abind(pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)], along=2) - truth_avg[2])^2, na.rm = T)
        result[i,3] = result[i,1] + result[i,2]
      }
    }
    
    mse$left = mean(result[,1], na.rm = T)
    mse$right = mean(result[,2], na.rm = T)
    mse$total = mean(result[,3], na.rm = T)
    
    return(mse)
  }
  
  if(effect_size == 0.25){
    mse = list()
    
    truth_avg = numeric(2)
    truth_avg[1] = (sum(truth[3:(dimension/2 - 2), 3:(dimension/2 - 2)]) + sum(truth[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)])+ sum(truth[(dimension/2 + 3):(dimension-2),(dimension/2 + 3):(dimension-2)])) / (((dimension/2 - 4)^2)*3)
    truth_avg[2] = sum(truth[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))]) / (((dimension/2 - 4)^2)*1)
    
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        result[i,1] = mean((abind(pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],pred[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)],pred[(dimension/2 + 3):(dimension-2),(dimension/2 + 3):(dimension-2)],along=2)-truth_avg[1])^2,na.rm = T)
        result[i,2] = mean((abind(pred[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))],along=2)  -truth_avg[2])^2, na.rm = T)
        result[i,3] = result[i,1] + result[i,2]
      }
    }
    
    mse$left = mean(result[,1], na.rm = T)
    mse$right = mean(result[,2], na.rm = T)
    mse$total = mean(result[,3], na.rm = T)
    
    return(mse)
    
  }
  
  if(effect_size == 0.75){
    mse = list()
    
    truth_avg = numeric(2)
    truth_avg[1] = (sum(truth[3:(dimension/2 - 2), 3:(dimension/2 - 2)])) / (((dimension/2 - 4)^2)*1)
    truth_avg[2] = (sum(truth[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])+sum(truth[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)]) +sum(truth[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)])) / (((dimension/2 - 4)^2)*3)
    
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        result[i,1] = mean((abind(pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],along=2) - truth_avg[1])^2, na.rm = T)
        result[i,2] = mean((abind(pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2+3):(dimension - 2), 3:(dimension/2 - 2)],along=2) - truth_avg[2])^2,na.rm = T)
        result[i,3] = result[i,1] + result[i,2]
      }
    }
    
    mse$left = mean(result[,1], na.rm = T)
    mse$right = mean(result[,2], na.rm = T)
    mse$total = mean(result[,3], na.rm = T)
    
    return(mse)
  }
  
  if(effect_size == 1.0){
    mse = list()
    
    truth_avg = (sum(truth[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])+sum(truth[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)]) +sum(truth[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)])+sum(truth[3:(dimension/2 - 2), 3:(dimension/2 - 2)])) / (((dimension/2 - 4)^2)*4)
    
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        result[i,1] = 0
        result[i,2] = mean((abind(pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)],pred[(dimension/2+3):(dimension - 2), 3:(dimension/2 - 2)],pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)],along=2) - truth_avg)^2,na.rm = T)
        result[i,3] = result[i,1] + result[i,2]
      }
    }
    mse$left = mean(result[,1], na.rm = T)
    mse$right = mean(result[,2], na.rm = T)
    mse$total = mean(result[,3], na.rm = T)
    return(mse)
  }
  
}

calculate_variance_beta = function(truth, prediction, n_sim = 1, effect_size = 0.5){
  
  if(effect_size == 0.5){
    variance = list()
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        
        result[i,1] = mean(abind((pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2 , (pred[(dimension/2 +3):(dimension - 2), 3:(dimension/2 - 2)])^2,along=2), na.rm=T)
        result[i,2] = mean(abind((pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])^2,(pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])^2,along=2), na.rm=T)
        result[i,3] = mean(abind((pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2,(pred[(dimension/2 +3):(dimension - 2), 3:(dimension/2 - 2)])^2,
                                 (pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])^2,(pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])^2,along=2), na.rm = T)
      }
    }
    
    variance$left = mean(result[,1], na.rm = T)
    variance$right = mean(result[,2], na.rm = T)
    variance$total = mean(result[,3], na.rm = T)
    
  }
  
  if(effect_size == 0.25){
    variance = list()
    
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        
        result[i,1] = mean(abind((pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2 , (pred[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)])^2, (pred[(dimension/2+3):(dimension-2),(dimension/2 + 3):(dimension-2)])^2,along=2), na.rm=T)
        result[i,2] = mean(abind((pred[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))])^2,along=2), na.rm = T) 
        result[i,3] = mean(abind((pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2 ,
                                 (pred[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)])^2 ,
                                 (pred[(dimension/2+3):(dimension-2),(dimension/2 + 3):(dimension-2)])^2 ,
                                 (pred[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))])^2,along=2), na.rm = T)
      }
    }
    
    variance$left = mean(result[,1], na.rm = T)
    variance$right = mean(result[,2], na.rm = T)
    variance$total = mean(result[,3], na.rm = T)
    
  }
  
  if(effect_size == 0.75){
    variance = list()
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        
        result[i,1] = mean(abind((pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2,along=2), na.rm = T)
        result[i,2] = mean(abind((pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])^2,(pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])^2,(pred[(dimension/2 +3):(dimension - 2), 3:(dimension/2 - 2)])^2,along=2), na.rm = T)
        result[i,3] = mean(abind((pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2,(pred[(dimension/2 +3):(dimension - 2), 3:(dimension/2 - 2)])^2, 
                                 (pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])^2,(pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])^2,along=2), na.rm = T)
      }
    }
    
    variance$left = mean(result[,1], na.rm = T)
    variance$right = mean(result[,2], na.rm = T)
    variance$total = mean(result[,3], na.rm = T)
    
  }
  
  if(effect_size == 1.0){
    variance = list()
    result = matrix(NA, nrow = n_sim, ncol = 3)
    for(i in 1:n_sim){
      pred = prediction[i,]
      if(!is.na(sum(pred))){
        
        pred = matrix(pred, nrow = dimension, ncol = dimension)
        
        result[i,1] = 0
        result[i,2] = mean(abind((pred[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)])^2,(pred[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])^2,(pred[(dimension/2 +3):(dimension - 2), 3:(dimension/2 - 2)])^2 ,( pred[3:(dimension/2 - 2), 3:(dimension/2 - 2)])^2,along=2), na.rm = T)
        result[i,3] = result[i,2]
      }
    }
    variance$left = mean(result[,1], na.rm = T)
    variance$right = mean(result[,2], na.rm = T)
    variance$total = mean(result[,3], na.rm = T)
    
  }
  return(variance)
}

calculate_mse_bias_variance = function(truth, prediction, prediction_sd, n_sim = 1 , effect_size = 0.5, type = 'beta'){
  
  if(type == 'y'){
    variance = calculate_variance(truth, prediction, n_sim, effect_size)
    bias = calculate_bias(truth, prediction, n_sim, effect_size)
    
    mse = list()
    mse$left = variance$left + bias$left^2
    mse$right = variance$right + bias$right^2
    mse$total = variance$total + bias$total^2
  }
  
  if(type == 'beta'){
    variance = calculate_variance_beta(truth=truth, prediction=prediction_sd, n_sim=n_sim,effect_size =  effect_size)
    bias = calculate_bias(truth, prediction, n_sim, effect_size)
    
    mse = list()
    mse$left = variance$left + bias$left^2
    mse$right = variance$right + bias$right^2
    mse$total = variance$total + bias$total^2
  }
  
  return(mse)
  
}


read_in_results_P2 = function(model, X, v0_list = seq(0.005, 0.15, 0.005), idx = 1:100){
  
  if(model == "BLESS_DPE_MCAR"){
    
    result = list()
    
    result$beta0_list = read.csv("beta0.csv")[idx,2:(M+1)]
    result$beta0_list = as.numeric(unlist(result$beta0_list))
    result$beta0_list = matrix(result$beta0_list, nrow = n_sim, ncol = M)
    
    result$Beta1_list = read.csv("Beta1.csv")[idx,2:(M+1)]
    result$Beta1_list = as.numeric(unlist(result$Beta1_list))
    result$Beta1_list = matrix(result$Beta1_list, nrow = n_sim, ncol = M)
    
    result$Beta2_list = read.csv("Beta2.csv")[idx,2:(M+1)]
    result$Beta2_list = as.numeric(unlist(result$Beta2_list))
    result$Beta2_list = matrix(result$Beta2_list, nrow = n_sim, ncol = M)
    
    n_sim = dim(result$Beta1_list)[1]
    M = dim(result$Beta1_list)[2]
    N = dim(X)[1]
    y_pred_1f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_1m = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2m = matrix(NA, nrow = n_sim, ncol = M)
    for (i in 1:n_sim) {
      tryCatch({
        Eta = cbind(rep(1,N), X) %*% rbind(as.vector(result$beta0_list[i,]), as.vector(result$Beta1_list[i,]), as.vector(result$Beta2_list[i,]))
        y_pred_1m[i,] = pnorm(Eta[(N*0.5 + 1),])
        y_pred_1f[i,] = pnorm(Eta[1,])
        y_pred_2m[i,] = pnorm(Eta[(N*0.75 + 1),])
        y_pred_2f[i,] = pnorm(Eta[(N*0.25 + 1),])
      }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
    }
    result$y_1m = y_pred_1m
    result$y_1f = y_pred_1f 
    result$y_2m = y_pred_2m
    result$y_2f = y_pred_2f 
    
    result$marginal_list = read.csv("marginal.csv")[,2:(length(v0)+1)]
    n_dim1 = dim(result$marginal_list)[1]
    result$marginal_list = as.numeric(unlist(result$marginal_list))
    result$marginal_list = matrix(result$marginal_list, n_dim1, length(v0))
    result$v0 = v0
    result$optimal_v0 = read.csv("optimal_v0.csv")[,2]
    result$optimal_v0_idx = read.csv("optimal_v0_idx.csv")[,2]
    
    result$ELBO_list = read.csv("ELBO.csv")[idx,2]
    result$ELBO_list = as.numeric(unlist(result$ELBO_list))
    result$ELBO_list = matrix(result$ELBO_list, nrow = n_sim, ncol = 1)
    
    result$gamma_list = read.csv("gamma1.csv")[idx,2:(M+1)]
    result$gamma_list = as.numeric(unlist(result$gamma_list))
    result$gamma_list = matrix(result$gamma_list, nrow = n_sim, ncol = M)
    
    result$gamma2_list = read.csv("gamma2.csv")[idx,2:(M+1)]
    result$gamma2_list = as.numeric(unlist(result$gamma2_list))
    result$gamma2_list = matrix(result$gamma2_list, nrow = n_sim, ncol = M)
    
    result$gamma21_list = read.csv("gamma21.csv")[idx,2:(M+1)]
    result$gamma21_list = as.numeric(unlist(result$gamma21_list))
    result$gamma21_list = matrix(result$gamma21_list, nrow = n_sim, ncol = M)
    
    result$gamma22_list = read.csv("gamma22.csv")[idx,2:(M+1)]
    result$gamma22_list = as.numeric(unlist(result$gamma22_list))
    result$gamma22_list = matrix(result$gamma22_list, nrow = n_sim, ncol = M)
    
    result$theta1_list = read.csv("theta1.csv")[idx,2:(M+1)]
    result$theta1_list = as.numeric(unlist(result$theta1_list))
    result$theta1_list = matrix(result$theta1_list, nrow = n_sim, ncol = M)
    
    result$theta2_list = read.csv("theta2.csv")[idx,2:(M+1)]
    result$theta2_list = as.numeric(unlist(result$theta2_list))
    result$theta2_list = matrix(result$theta2_list, nrow = n_sim, ncol = M)
    
    result$xi1_list = read.csv("xi1.csv")[idx,2:(M+1)]
    result$xi1_list = as.numeric(unlist(result$xi1_list))
    result$xi1_list = matrix(result$xi1_list, nrow = n_sim, ncol = M)
    
    result$xi2_list = read.csv("xi2.csv")[idx,2:(M+1)]
    result$xi2_list = as.numeric(unlist(result$xi2_list))
    result$xi2_list = matrix(result$xi2_list, nrow = n_sim, ncol = M)
    
    result$Sigma_Inv_list = read.csv("Sigma_Inv.csv")[idx,2:(P*P+1)]
    result$Sigma_Inv_list = as.numeric(unlist(result$Sigma_Inv_list))
    result$Sigma_Inv_list = matrix(result$Sigma_Inv_list, nrow = n_sim, ncol = P*P)
    
    result$t1 = result$gamma_list
    result$t1[result$t1 > 0.5] = 1
    result$t1[result$t1 <= 0.5] = 0
    
    result$t2 = result$gamma2_list
    result$t2[result$t2 > 0.5] = 1
    result$t2[result$t2 <= 0.5] = 0
    
    result$beta0_sd = matrix(NA, n_sim, M)
    for(i in 1:n_sim){
      result$beta0_sd[i,] = rep(sqrt(1/(N + 1/(10^2))), M)
    }
    result$Beta1_sd = matrix(NA, n_sim, M)
    result$Beta2_sd = matrix(NA, n_sim, M)
    for(i in 1:n_sim){
      for(j in 1:M){
        exp_gamma2 = c(result$gamma21_list[i,j],result$gamma22_list[i,j])
        x = solve(t(X)%*%X + diag(exp_gamma2,P))
        result$Beta1_sd[i,j] = sqrt(x[1,1])
        result$Beta2_sd[i,j] = sqrt(x[2,2])
      }
    }    
  } 
  
  if(model == "BLESS_DPE_CAR"){
    
    result = list()
    
    result$beta0_list = read.csv("beta0.csv")[idx,2:(M+1)]
    result$beta0_list = as.numeric(unlist(result$beta0_list))
    result$beta0_list = matrix(result$beta0_list, nrow = n_sim, ncol = M)
    
    result$Beta1_list = read.csv("Beta1.csv")[idx,2:(M+1)]
    result$Beta1_list = as.numeric(unlist(result$Beta1_list))
    result$Beta1_list = matrix(result$Beta1_list, nrow = n_sim, ncol = M)
    
    result$Beta2_list = read.csv("Beta2.csv")[idx,2:(M+1)]
    result$Beta2_list = as.numeric(unlist(result$Beta2_list))
    result$Beta2_list = matrix(result$Beta2_list, nrow = n_sim, ncol = M)
    
    n_sim = dim(result$Beta1_list)[1]
    M = dim(result$Beta1_list)[2]
    N = dim(X)[1]
    y_pred_1f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_1m = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2m = matrix(NA, nrow = n_sim, ncol = M)
    for (i in 1:n_sim) {
      tryCatch({
        Eta = cbind(rep(1,N), X) %*% rbind(as.vector(result$beta0_list[i,]), as.vector(result$Beta1_list[i,]), as.vector(result$Beta2_list[i,]))
        y_pred_1m[i,] = pnorm(Eta[(N*0.5 + 1),])
        y_pred_1f[i,] = pnorm(Eta[1,])
        y_pred_2m[i,] = pnorm(Eta[(N*0.75 + 1),])
        y_pred_2f[i,] = pnorm(Eta[(N*0.25 + 1),])
      }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
    }
    result$y_1m = y_pred_1m
    result$y_1f = y_pred_1f 
    result$y_2m = y_pred_2m
    result$y_2f = y_pred_2f 
    
    result$marginal_list = read.csv("marginal.csv")[,2:(length(v0)+1)]
    n_dim1 = dim(result$marginal_list)[1]
    result$marginal_list = as.numeric(unlist(result$marginal_list))
    result$marginal_list = matrix(result$marginal_list, n_dim1, length(v0))
    result$v0 = v0
    result$optimal_v0 = read.csv("optimal_v0.csv")[,2]
    result$optimal_v0_idx = read.csv("optimal_v0_idx.csv")[,2]
    
    result$ELBO_list = read.csv("ELBO.csv")[idx,2]
    result$ELBO_list = as.numeric(unlist(result$ELBO_list))
    result$ELBO_list = matrix(result$ELBO_list, nrow = n_sim, ncol = 1)
    
    result$gamma_list = read.csv("gamma1.csv")[idx,2:(M+1)]
    result$gamma_list = as.numeric(unlist(result$gamma_list))
    result$gamma_list = matrix(result$gamma_list, nrow = n_sim, ncol = M)
    
    result$gamma2_list = read.csv("gamma2.csv")[idx,2:(M+1)]
    result$gamma2_list = as.numeric(unlist(result$gamma2_list))
    result$gamma2_list = matrix(result$gamma2_list, nrow = n_sim, ncol = M)
    
    result$gamma21_list = read.csv("gamma21.csv")[idx,2:(M+1)]
    result$gamma21_list = as.numeric(unlist(result$gamma21_list))
    result$gamma21_list = matrix(result$gamma21_list, nrow = n_sim, ncol = M)
    
    result$gamma22_list = read.csv("gamma22.csv")[idx,2:(M+1)]
    result$gamma22_list = as.numeric(unlist(result$gamma22_list))
    result$gamma22_list = matrix(result$gamma22_list, nrow = n_sim, ncol = M)
    
    result$theta1_list = read.csv("theta1.csv")[idx,2:(M+1)]
    result$theta1_list = as.numeric(unlist(result$theta1_list))
    result$theta1_list = matrix(result$theta1_list, nrow = n_sim, ncol = M)
    
    result$xi1_list = read.csv("xi1.csv")[idx,2:(M+1)]
    result$xi1_list = as.numeric(unlist(result$xi1_list))
    result$xi1_list = matrix(result$xi1_list, nrow = n_sim, ncol = M)
    
    result$Sigma_Inv_list = read.csv("Sigma_Inv.csv")[idx,2]
    result$Sigma_Inv_list = as.numeric(unlist(result$Sigma_Inv_list))
    result$Sigma_Inv_list = matrix(result$Sigma_Inv_list, nrow = n_sim, ncol = 1)
    
    result$t1 = result$gamma_list
    result$t1[result$t1 > 0.5] = 1
    result$t1[result$t1 <= 0.5] = 0
    
    result$t2 = result$gamma2_list
    result$t2[result$t2 > 0.5] = 1
    result$t2[result$t2 <= 0.5] = 0
    
    result$beta0_sd = matrix(NA, n_sim, M)
    for(i in 1:n_sim){
      result$beta0_sd[i,] = rep(sqrt(1/(N + 1/(10^2))), M)
    }
    result$Beta1_sd = matrix(NA, n_sim, M)
    result$Beta2_sd = matrix(NA, n_sim, M)
    for(i in 1:n_sim){
      for(j in 1:M){
        exp_gamma2 = c(result$gamma21_list[i,j],result$gamma22_list[i,j])
        x = solve(t(X)%*%X + diag(exp_gamma2,P))
        result$Beta1_sd[i,j] = sqrt(x[1,1])
        result$Beta2_sd[i,j] = sqrt(x[2,2])
      }
    }
  } 
  
  if(model == "BLESS_DPE_PxCAR"){
    
    result = list()
    
    result$beta0_list = read.csv("beta0.csv")[idx,2:(M+1)]
    result$beta0_list = as.numeric(unlist(result$beta0_list))
    result$beta0_list = matrix(result$beta0_list, nrow = n_sim, ncol = M)
    
    result$Beta1_list = read.csv("Beta1.csv")[idx,2:(M+1)]
    result$Beta1_list = as.numeric(unlist(result$Beta1_list))
    result$Beta1_list = matrix(result$Beta1_list, nrow = n_sim, ncol = M)
    
    result$Beta2_list = read.csv("Beta2.csv")[idx,2:(M+1)]
    result$Beta2_list = as.numeric(unlist(result$Beta2_list))
    result$Beta2_list = matrix(result$Beta2_list, nrow = n_sim, ncol = M)
    
    n_sim = dim(result$Beta1_list)[1]
    M = dim(result$Beta1_list)[2]
    N = dim(X)[1]
    y_pred_1f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_1m = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2m = matrix(NA, nrow = n_sim, ncol = M)
    for (i in 1:n_sim) {
      tryCatch({
        Eta = cbind(rep(1,N), X) %*% rbind(as.vector(result$beta0_list[i,]), as.vector(result$Beta1_list[i,]), as.vector(result$Beta2_list[i,]))
        y_pred_1m[i,] = pnorm(Eta[(N*0.5 + 1),])
        y_pred_1f[i,] = pnorm(Eta[1,])
        y_pred_2m[i,] = pnorm(Eta[(N*0.75 + 1),])
        y_pred_2f[i,] = pnorm(Eta[(N*0.25 + 1),])
      }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
    }
    result$y_1m = y_pred_1m
    result$y_1f = y_pred_1f 
    result$y_2m = y_pred_2m
    result$y_2f = y_pred_2f 
    
    result$marginal_list = read.csv("marginal.csv")[,2:(length(v0)+1)]
    n_dim1 = dim(result$marginal_list)[1]
    result$marginal_list = as.numeric(unlist(result$marginal_list))
    result$marginal_list = matrix(result$marginal_list, n_dim1, length(v0))
    result$v0 = v0
    result$optimal_v0 = read.csv("optimal_v0.csv")[,2]
    result$optimal_v0_idx = read.csv("optimal_v0_idx.csv")[,2]
    
    result$ELBO_list = read.csv("ELBO.csv")[idx,2]
    result$ELBO_list = as.numeric(unlist(result$ELBO_list))
    result$ELBO_list = matrix(result$ELBO_list, nrow = n_sim, ncol = 1)
    
    result$gamma_list = read.csv("gamma1.csv")[idx,2:(M+1)]
    result$gamma_list = as.numeric(unlist(result$gamma_list))
    result$gamma_list = matrix(result$gamma_list, nrow = n_sim, ncol = M)
    
    result$gamma2_list = read.csv("gamma2.csv")[idx,2:(M+1)]
    result$gamma2_list = as.numeric(unlist(result$gamma2_list))
    result$gamma2_list = matrix(result$gamma2_list, nrow = n_sim, ncol = M)
    
    result$gamma21_list = read.csv("gamma21.csv")[idx,2:(M+1)]
    result$gamma21_list = as.numeric(unlist(result$gamma21_list))
    result$gamma21_list = matrix(result$gamma21_list, nrow = n_sim, ncol = M)
    
    result$gamma22_list = read.csv("gamma22.csv")[idx,2:(M+1)]
    result$gamma22_list = as.numeric(unlist(result$gamma22_list))
    result$gamma22_list = matrix(result$gamma22_list, nrow = n_sim, ncol = M)
    
    result$theta1_list = read.csv("theta1.csv")[idx,2:(M+1)]
    result$theta1_list = as.numeric(unlist(result$theta1_list))
    result$theta1_list = matrix(result$theta1_list, nrow = n_sim, ncol = M)
    
    result$theta2_list = read.csv("theta2.csv")[idx,2:(M+1)]
    result$theta2_list = as.numeric(unlist(result$theta2_list))
    result$theta2_list = matrix(result$theta2_list, nrow = n_sim, ncol = M)
    
    result$xi1_list = read.csv("xi1.csv")[idx,2:(M+1)]
    result$xi1_list = as.numeric(unlist(result$xi1_list))
    result$xi1_list = matrix(result$xi1_list, nrow = n_sim, ncol = M)
    
    result$xi2_list = read.csv("xi2.csv")[idx,2:(M+1)]
    result$xi2_list = as.numeric(unlist(result$xi2_list))
    result$xi2_list = matrix(result$xi2_list, nrow = n_sim, ncol = M)
    
    result$Sigma_Inv_list = read.csv("Sigma_Inv.csv")[idx,2]
    result$Sigma_Inv_list = as.numeric(unlist(result$Sigma_Inv_list))
    result$Sigma_Inv_list = matrix(result$Sigma_Inv_list, nrow = n_sim, ncol = P)
    
    result$t1 = result$gamma_list
    result$t1[result$t1 > 0.5] = 1
    result$t1[result$t1 <= 0.5] = 0
    
    result$t2 = result$gamma2_list
    result$t2[result$t2 > 0.5] = 1
    result$t2[result$t2 <= 0.5] = 0
    
    result$beta0_sd = matrix(NA, n_sim, M)
    for(i in 1:n_sim){
      result$beta0_sd[i,] = rep(sqrt(1/(N + 1/(10^2))), M)
    }
    result$Beta1_sd = matrix(NA, n_sim, M)
    result$Beta2_sd = matrix(NA, n_sim, M)
    for(i in 1:n_sim){
      for(j in 1:M){
        exp_gamma2 = c(result$gamma21_list[i,j],result$gamma22_list[i,j])
        x = solve(t(X)%*%X + diag(exp_gamma2,P))
        result$Beta1_sd[i,j] = sqrt(x[1,1])
        result$Beta2_sd[i,j] = sqrt(x[2,2])
      }
    }    
  } 
  
  if(model == 'BSGLMM'){
    
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
    
    n_sj = n_neighbors(dim1, dim2)
    
    result = list()
    
    result$beta0_list = read.csv("beta0.csv")[idx,2:(M+1)]
    result$beta0_list = as.numeric(unlist(result$beta0_list))
    result$beta0_list = matrix(result$beta0_list, nrow = n_sim, ncol = M)
    
    result$Beta1_list = read.csv("beta1.csv")[idx,2:(M+1)]
    result$Beta1_list = as.numeric(unlist(result$Beta1_list))
    result$Beta1_list = matrix(result$Beta1_list, nrow = n_sim, ncol = M)
    
    result$Beta2_list = read.csv("beta2.csv")[idx,2:(M+1)]
    result$Beta2_list = as.numeric(unlist(result$Beta2_list))
    result$Beta2_list = matrix(result$Beta2_list, nrow = n_sim, ncol = M)
    
    M = dim(result$Beta1_list)[2]
    N = dim(X)[1]
    y_pred_1f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_1m = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2m = matrix(NA, nrow = n_sim, ncol = M)
    for (i in 1:n_sim) {
      Eta = cbind(rep(1,N), X) %*% rbind(as.vector(result$beta0_list[i,]), as.vector(result$Beta1_list[i,]), as.vector(result$Beta2_list[i,]))
      y_pred_1m[i,] = pnorm(Eta[(N*0.5 + 1),])
      y_pred_1f[i,] = pnorm(Eta[1,])
      y_pred_2m[i,] = pnorm(Eta[(N*0.75 + 1),])
      y_pred_2f[i,] = pnorm(Eta[(N*0.25 + 1),])
    }
    result$y_1m = y_pred_1m
    result$y_1f = y_pred_1f 
    result$y_2m = y_pred_2m
    result$y_2f = y_pred_2f 
    
    #result$ELBO_list = read.csv("ELBO.csv")[,2]
    #result$ELBO_list = as.numeric(unlist(result$ELBO_list))
    #result$ELBO_list = matrix(result$ELBO_list, nrow = n_sim, ncol = 1)
    
    result$Sigma_Inv_list = read.csv("Sigma_Inv.csv")[idx,2:(P*P + 1)]
    result$Sigma_Inv_list = as.numeric(unlist(result$Sigma_Inv_list))
    result$Sigma_Inv_list = matrix(result$Sigma_Inv_list, nrow = n_sim, ncol = P*P)
    
    M = dim(result$Beta1_list)[2]
    dimension = sqrt(M)
    n_sj = matrix(4, nrow = dimension, ncol = dimension)
    n_sj[1,1] = n_sj[1,dimension] = n_sj[dimension,1] = n_sj[dimension,dimension] = 2
    n_sj[2:(dimension-1),1] = n_sj[2:(dimension-1),dimension] = n_sj[1,2:(dimension-1)] = n_sj[dimension,2:(dimension-1)] = 3
    n_sj = as.vector(n_sj)
    n_sim = dim(result$Beta1_list)[1]
  
    result$beta0_sd = read.csv("beta0_sd.csv")[idx,2:(M+1)]
    result$beta0_sd = as.numeric(unlist(result$beta0_sd))
    result$beta0_sd = matrix(result$beta0_sd, nrow = n_sim, ncol = M)
    
    result$Beta1_sd = read.csv("beta1_sd.csv")[idx,2:(M+1)]
    result$Beta1_sd = as.numeric(unlist(result$Beta1_sd))
    result$Beta1_sd = matrix(result$Beta1_sd, nrow = n_sim, ncol = M)
    
    result$Beta2_sd = read.csv("beta2_sd.csv")[idx,2:(M+1)]
    result$Beta2_sd = as.numeric(unlist(result$Beta2_sd))
    result$Beta2_sd = matrix(result$Beta2_sd, nrow = n_sim, ncol = M)
    
    result$t1 = read.csv("t_beta1.csv")[idx,2:(M+1)]
    result$t1 = as.numeric(unlist(result$t1))
    t = result$t1 = matrix(result$t1, nrow = n_sim, ncol = M)
    result$t1[abs(result$t1) > 1.96] = 1
    result$t1[result$t1 != 1] = 0
    
    t = 2*pnorm(q=-abs(t))
    result$t1_corrected = matrix(0, n_sim, M)
    for(i in 1:n_sim){
      result$t1_corrected[i,] = p.adjust(t[i,], method='fdr')
    }
    result$t1_corrected[result$t1_corrected <= 0.05] = 1
    result$t1_corrected[result$t1_corrected != 1] = 0
    
    result$t2 = read.csv("t_beta2.csv")[idx,2:(M+1)]
    result$t2 = as.numeric(unlist(result$t2))
    t = result$t2 = matrix(result$t2, nrow = n_sim, ncol = M)
    result$t2[abs(result$t2) > 1.96] = 1
    result$t2[result$t2 != 1] = 0
    
    t = 2*pnorm(q=-abs(t))
    result$t2_corrected = matrix(0, n_sim, M)
    for(i in 1:n_sim){
      result$t2_corrected[i,] = p.adjust(t[i,], method = 'fdr')
    }
    result$t2_corrected[result$t2_corrected <= 0.05] = 1
    result$t2_corrected[result$t2_corrected != 1] = 0
    
    print('BSGLMM fine')
  }
  
  if(model == 'Firth'){
    result = list()
    
    result$beta0_list = read.csv(paste0("beta0.csv"))[idx,2:(M+1)]
    result$beta0_list = as.numeric(unlist(result$beta0_list))
    result$beta0_list = matrix(result$beta0_list, nrow = n_sim, ncol = M)
    
    result$Beta1_list = read.csv(paste0("Beta1.csv"))[idx,2:(M+1)]
    result$Beta1_list = as.numeric(unlist(result$Beta1_list))
    result$Beta1_list = matrix(result$Beta1_list, nrow = n_sim, ncol = M)
    
    result$Beta2_list = read.csv(paste0("Beta2.csv"))[idx,2:(M+1)]
    result$Beta2_list = as.numeric(unlist(result$Beta2_list))
    result$Beta2_list = matrix(result$Beta2_list, nrow = n_sim, ncol = M)
    
    result$y_1m = read.csv(paste0("y_pred_1m.csv"))[idx,2:(M+1)]
    result$y_1m = as.numeric(unlist(result$y_1m))
    result$y_1m = matrix(result$y_1m, nrow = n_sim, ncol = M)
    
    result$y_1f = read.csv(paste0("y_pred_1f.csv"))[idx,2:(M+1)]
    result$y_1f = as.numeric(unlist(result$y_1f))
    result$y_1f = matrix(result$y_1f, nrow = n_sim, ncol = M)
    
    result$y_2m = read.csv(paste0("y_pred_2m.csv"))[idx,2:(M+1)]
    result$y_2m = as.numeric(unlist(result$y_2m))
    result$y_2m = matrix(result$y_2m, nrow = n_sim, ncol = M)
    
    result$y_2f = read.csv(paste0("y_pred_2f.csv"))[idx,2:(M+1)]
    result$y_2f = as.numeric(unlist(result$y_2f))
    result$y_2f = matrix(result$y_2f, nrow = n_sim, ncol = M)
    
    result$beta0_sd = read.csv(paste0("std_error_beta0.csv"))[idx,2:(M+1)]
    result$beta0_sd = as.numeric(unlist(result$beta0_sd))
    result$beta0_sd = matrix(result$beta0_sd, nrow = n_sim, ncol = M)
    
    result$Beta1_sd = read.csv(paste0("std_error_Beta1.csv"))[idx,2:(M+1)]
    result$Beta1_sd = as.numeric(unlist(result$Beta1_sd))
    result$Beta1_sd = matrix(result$Beta1_sd, nrow = n_sim, ncol = M)
    
    result$Beta2_sd = read.csv(paste0("std_error_Beta2.csv"))[idx,2:(M+1)]
    result$Beta2_sd = as.numeric(unlist(result$Beta2_sd))
    result$Beta2_sd = matrix(result$Beta2_sd, nrow = n_sim, ncol = M)
    
    t = result$Beta1_list /result$Beta1_sd
    t = 2*pnorm(q=-abs(t))
    result$t1 = matrix(0, n_sim, M)
    for(i in 1:n_sim){
      result$t1[i,] = p.adjust(t[i,], method='fdr')
    }
    result$t1[result$t1 <= 0.05] = 1
    result$t1[result$t1 != 1] = 0
    
    t = result$Beta2_list / result$Beta2_sd
    t = 2*pnorm(q=-abs(t))
    result$t2 = matrix(0, n_sim, M)
    for(i in 1:n_sim){
      result$t2[i,] = p.adjust(t[i,], method = 'fdr')
    }
    result$t2[result$t2 <= 0.05] = 1
    result$t2[result$t2 != 1] = 0
    
    result$t1_uncorrected = result$Beta1_list /result$Beta1_sd
    result$t1_uncorrected[abs(result$t1_uncorrected) > 1.96] = 1
    result$t1_uncorrected[result$t1_uncorrected != 1] = 0
    
    result$t2_uncorrected = result$Beta2_list /result$Beta2_sd
    result$t2_uncorrected[abs(result$t2_uncorrected) > 1.96] = 1
    result$t2_uncorrected[result$t2_uncorrected != 1] = 0
    
    print('Firth fine')
  }
  
  if(model == "BB_BLESS"){
    
    path = path_results_BB_BLESS
    
    result = list()
    
    result$beta0_list = read.csv(paste0(path,"beta0_mean.csv"))[,2:(M+1)]
    result$beta0_list = as.numeric(unlist(result$beta0_list))
    result$beta0_list = matrix(result$beta0_list, nrow = n_sim, ncol = M)
    
    result$beta0_sd = read.csv(paste0(path,"beta0_sd.csv"))[,2:(M+1)]
    result$beta0_sd = as.numeric(unlist(result$beta0_sd))
    result$beta0_sd = matrix(result$beta0_sd, nrow = n_sim, ncol = M)
    
    result$Beta1_list = read.csv(paste0(path,"beta1_mean.csv"))[,2:(M+1)]
    result$Beta1_list = as.numeric(unlist(result$Beta1_list))
    result$Beta1_list = matrix(result$Beta1_list, nrow = n_sim, ncol = M)
    
    result$Beta1_sd = read.csv(paste0(path,"beta1_sd.csv"))[,2:(M+1)]
    result$Beta1_sd = as.numeric(unlist(result$Beta1_sd))
    result$Beta1_sd = matrix(result$Beta1_sd, nrow = n_sim, ncol = M)
    
    result$Beta2_list = read.csv(paste0(path,"beta2_mean.csv"))[,2:(M+1)]
    result$Beta2_list = as.numeric(unlist(result$Beta2_list))
    result$Beta2_list = matrix(result$Beta2_list, nrow = n_sim, ncol = M)
    
    result$Beta2_sd = read.csv(paste0(path,"beta2_sd.csv"))[,2:(M+1)]
    result$Beta2_sd = as.numeric(unlist(result$Beta2_sd))
    result$Beta2_sd = matrix(result$Beta2_sd, nrow = n_sim, ncol = M)
    
    n_sim = dim(result$Beta1_list)[1]
    M = dim(result$Beta1_list)[2]
    N = dim(X)[1]
    y_pred_1f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_1m = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2m = matrix(NA, nrow = n_sim, ncol = M)
    for (i in 1:n_sim) {
      tryCatch({
        Eta = cbind(rep(1,N), X) %*% rbind(as.vector(result$beta0_list[i,]), as.vector(result$Beta1_list[i,]), as.vector(result$Beta2_list[i,]))
        y_pred_1m[i,] = pnorm(Eta[(N*0.5 + 1),])
        y_pred_1f[i,] = pnorm(Eta[1,])
        y_pred_2m[i,] = pnorm(Eta[(N*0.75 + 1),])
        y_pred_2f[i,] = pnorm(Eta[(N*0.25 + 1),])
      }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
    }
    result$y_1m = y_pred_1m
    result$y_1f = y_pred_1f 
    result$y_2m = y_pred_2m
    result$y_2f = y_pred_2f 
    
    result$theta1_list = read.csv(paste0(path,"theta1.csv"))[,2:(M+1)]
    result$theta1_list = as.numeric(unlist(result$theta1_list))
    result$theta1_list = matrix(result$theta1_list, nrow = n_sim, ncol = M)
    
    result$theta2_list = read.csv(paste0(path,"theta2.csv"))[,2:(M+1)]
    result$theta2_list = as.numeric(unlist(result$theta2_list))
    result$theta2_list = matrix(result$theta2_list, nrow = n_sim, ncol = M)
    
    result$Sigma_Inv_list = read.csv(paste0(path,"Sigma_Inv.csv"))[,2:(P*P+1)]
    result$Sigma_Inv_list = as.numeric(unlist(result$Sigma_Inv_list))
    result$Sigma_Inv_list = matrix(result$Sigma_Inv_list, nrow = n_sim, ncol = P*P)
    
    result$t1_list =  result$Beta1_list /  result$Beta1_sd
    result$t1_list = as.numeric(unlist(result$t1_list))
    result$t1_list = matrix(result$t1_list, nrow = n_sim, ncol = M)
    result$t1_list[abs(result$t1_list) > 1.96] = 1
    result$t1_list[result$t1_list != 1] = 0
    
    result$t2_list =  result$Beta2_list /  result$Beta2_sd
    result$t2_list = as.numeric(unlist(result$t2_list))
    result$t2_list = matrix(result$t2_list, nrow = n_sim, ncol = M)
    result$t2_list[abs(result$t2_list) > 1.96] = 1
    result$t2_list[result$t2_list != 1] = 0
    
    print('BB-BLESS fine')
  } 
  
  if(model == "Gibbs"){
    
    result = list()
    
    result$beta0_list = read.csv("beta0_mean.csv")[idx,2:(M+1)]
    result$beta0_list = as.numeric(unlist(result$beta0_list))
    result$beta0_list = matrix(result$beta0_list, nrow = n_sim, ncol = M)
    
    result$beta0_sd = read.csv("beta0_sd.csv")[idx,2:(M+1)]
    result$beta0_sd = as.numeric(unlist(result$beta0_sd))
    result$beta0_sd = matrix(result$beta0_sd, nrow = n_sim, ncol = M)
    
    result$Beta1_list = read.csv("beta1_mean.csv")[idx,2:(M+1)]
    result$Beta1_list = as.numeric(unlist(result$Beta1_list))
    result$Beta1_list = matrix(result$Beta1_list, nrow = n_sim, ncol = M)
    
    result$Beta1_sd = read.csv("beta1_sd.csv")[idx,2:(M+1)]
    result$Beta1_sd = as.numeric(unlist(result$Beta1_sd))
    result$Beta1_sd = matrix(result$Beta1_sd, nrow = n_sim, ncol = M)
    
    result$Beta2_list = read.csv("beta2_mean.csv")[idx,2:(M+1)]
    result$Beta2_list = as.numeric(unlist(result$Beta2_list))
    result$Beta2_list = matrix(result$Beta2_list, nrow = n_sim, ncol = M)
    
    result$Beta2_sd = read.csv("beta2_sd.csv")[idx,2:(M+1)]
    result$Beta2_sd = as.numeric(unlist(result$Beta2_sd))
    result$Beta2_sd = matrix(result$Beta2_sd, nrow = n_sim, ncol = M)
    
    n_sim = dim(result$Beta1_list)[1]
    M = dim(result$Beta1_list)[2]
    N = dim(X)[1]
    y_pred_1f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_1m = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2f = matrix(NA, nrow = n_sim, ncol = M)
    y_pred_2m = matrix(NA, nrow = n_sim, ncol = M)
    for (i in 1:n_sim) {
      tryCatch({
        Eta = cbind(rep(1,N), X) %*% rbind(as.vector(result$beta0_list[i,]), as.vector(result$Beta1_list[i,]), as.vector(result$Beta2_list[i,]))
        y_pred_1m[i,] = pnorm(Eta[(N*0.5 + 1),])
        y_pred_1f[i,] = pnorm(Eta[1,])
        y_pred_2m[i,] = pnorm(Eta[(N*0.75 + 1),])
        y_pred_2f[i,] = pnorm(Eta[(N*0.25 + 1),])
      }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
    }
    result$y_1m = y_pred_1m
    result$y_1f = y_pred_1f 
    result$y_2m = y_pred_2m
    result$y_2f = y_pred_2f 
    
    result$gamma_list = read.csv("gamma1.csv")[idx,2:(M+1)]
    result$gamma_list = as.numeric(unlist(result$gamma_list))
    result$gamma_list = matrix(result$gamma_list, nrow = n_sim, ncol = M)
    
    result$gamma2_list = read.csv("gamma2.csv")[idx,2:(M+1)]
    result$gamma2_list = as.numeric(unlist(result$gamma2_list))
    result$gamma2_list = matrix(result$gamma2_list, nrow = n_sim, ncol = M)
    
    result$theta1_list = read.csv("theta1.csv")[idx,2:(M+1)]
    result$theta1_list = as.numeric(unlist(result$theta1_list))
    result$theta1_list = matrix(result$theta1_list, nrow = n_sim, ncol = M)
    
    result$theta2_list = read.csv("theta2.csv")[idx,2:(M+1)]
    result$theta2_list = as.numeric(unlist(result$theta2_list))
    result$theta2_list = matrix(result$theta2_list, nrow = n_sim, ncol = M)
    
    result$xi1_list = read.csv("xi1.csv")[idx,2:(M+1)]
    result$xi1_list = as.numeric(unlist(result$xi1_list))
    result$xi1_list = matrix(result$xi1_list, nrow = n_sim, ncol = M)
    
    result$xi2_list = read.csv("xi2.csv")[idx,2:(M+1)]
    result$xi2_list = as.numeric(unlist(result$xi2_list))
    result$xi2_list = matrix(result$xi2_list, nrow = n_sim, ncol = M)
    
    result$Sigma_Inv_list = read.csv("Sigma_Inv.csv")[idx,2:(P*P+1)]
    result$Sigma_Inv_list = as.numeric(unlist(result$Sigma_Inv_list))
    result$Sigma_Inv_list = matrix(result$Sigma_Inv_list, nrow = n_sim, ncol = P*P)
    
    result$t1_list = read.csv("t_beta1.csv")[idx,2:(M+1)]
    result$t1_list = as.numeric(unlist(result$t1_list))
    result$t1_list = matrix(result$t1_list, nrow = n_sim, ncol = M)
    result$t1_list[abs(result$t1_list)>1.96] = 1
    result$t1_list[result$t1_list!=1] = 0
    
    result$t2_list = read.csv("t_beta2.csv")[idx,2:(M+1)]
    result$t2_list = as.numeric(unlist(result$t2_list))
    result$t2_list = matrix(result$t2_list, nrow = n_sim, ncol = M)
    result$t2_list[abs(result$t2_list)>1.96] = 1
    result$t2_list[result$t2_list!=1] = 0
    
    print('Gibbs fine')
    
  } 
  
  return(result)
}

evaluate_bias_mae_mse = function(truth, result, n_sim, quantity = 'all'){
  
  eval = list()
  
  if(quantity == 'beta0'){
    eval$bias = calculate_bias(truth$beta0, result$beta0_list, n_sim, effect_size = 1.0)
    eval$variance = calculate_variance_beta(truth$beta0, result$beta0_sd, n_sim, effect_size = 1.0)
    eval$mse = calculate_mse_bias_variance(truth$beta0, result$beta0_list, result$beta0_sd, n_sim, effect_size = 1.0, type = 'beta')
  }
  
  if(quantity == 'beta1'){
    eval$bias = calculate_bias(truth$beta1, result$Beta1_list, n_sim, effect_size = 0.5)
    eval$variance = calculate_variance_beta(truth$beta1, result$Beta1_sd, n_sim, effect_size = 0.5)
    eval$mse = calculate_mse_bias_variance(truth$beta1, result$Beta1_list, result$Beta1_sd, n_sim, effect_size = 0.5, type = 'beta')
  }
  
  if(quantity == 'beta2'){
    eval$bias = calculate_bias(truth$beta2, result$Beta2_list, n_sim, effect_size = 0.25)
    eval$variance = calculate_variance_beta(truth$beta2, result$Beta2_sd, n_sim,effect_size =  0.25)
    eval$mse = calculate_mse_bias_variance(truth$beta2, result$Beta2_list, result$Beta2_sd, n_sim, effect_size = 0.25, type = 'beta')
  }
  
  if(quantity == 'y_1m'){
    eval$bias = calculate_bias(truth$y_1m, result$y_1m, n_sim, effect_size = 1.0 )
    eval$variance = calculate_variance(truth$y_1m, result$y_1m, n_sim, effect_size = 1.0)
    eval$mse = calculate_mse(truth$y_1m, result$y_1m, n_sim,effect_size =  1.0)
  }
  
  if(quantity == 'y_1f'){
    eval$bias = calculate_bias(truth$y_1f, result$y_1f, n_sim, effect_size = 0.5)
    eval$variance = calculate_variance(truth$y_1f, result$y_1f, n_sim, effect_size = 0.5)
    eval$mse = calculate_mse(truth$y_1f, result$y_1f, n_sim, effect_size = 0.5)
  }
  
  if(quantity == 'y_2m'){
    eval$bias = calculate_bias(truth$y_2m, result$y_2m, n_sim, effect_size = 0.25)
    eval$variance = calculate_variance(truth$y_2m, result$y_2m, n_sim, effect_size = 0.25)
    eval$mse = calculate_mse(truth$y_2m, result$y_2m, n_sim, effect_size = 0.25)
  }
  
  if(quantity == 'y_2f'){
    eval$bias = calculate_bias(truth$y_2f, result$y_2f, n_sim, effect_size = 0.75)
    eval$variance = calculate_variance(truth$y_2f, result$y_2f, n_sim, effect_size = 0.75)
    eval$mse = calculate_mse(truth$y_2f, result$y_2f, n_sim,effect_size =  0.75)
  }
  return(eval)
}

evaluate_TPR_TDR_FPR_FDR = function(result, n_sim, effect_size){
  
  # TPR = TP / (TP + FN)
  # TDR = TP / (TP + FP) = 1 - FDR
  # FPR = FP / (FP + TN)
  # FDR = FP / (FP + TP)
  
  eval = list()
  if(effect_size == 0.5){
    t_list = result
    n_sim = dim(t_list)[1]
    M = dim(t_list)[2]
    dimension = sqrt(M)
    
    TPR = numeric(n_sim)
    TDR = numeric(n_sim)
    FPR = numeric(n_sim)
    FDR = numeric(n_sim)
    for(i in 1:n_sim){
      
      t = t_list[i,]
      if(!is.na(sum(t))){
        t = matrix(t, nrow = dimension, ncol = dimension)
        
        positives = sum(t[3:(dimension/2 - 2), 3:(dimension/2 - 2)]) + sum(t[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)]) + 
          sum(t[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)]) + sum(t[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])
        
        TPR[i] = (sum(t[3:(dimension/2 - 2),((dimension/2) + 3):(dimension-2)]) + sum(t[(dimension/2 + 3):(dimension-2),((dimension/2) + 3):(dimension-2)])) / (((dimension/2 - 4)^2)*2)
        TDR[i] = (sum(t[3:(dimension/2 - 2),((dimension/2) + 3):(dimension-2)]) + sum(t[(dimension/2 + 3):(dimension-2),((dimension/2) + 3):(dimension-2)])) / positives
        FPR[i] = (sum(t[3:(dimension/2 - 2),3:(dimension/2 - 2)]) +  sum(t[(dimension/2 + 3):(dimension-2),3:(dimension/2 - 2)]))/ (((dimension/2 - 4)^2)*2)
        FDR[i] = (sum(t[3:(dimension/2 - 2),3:(dimension/2 - 2)]) +  sum(t[(dimension/2 + 3):(dimension-2),3:(dimension/2 - 2)])) / positives
      }
    }
    
    idx = which(TPR == 0, arr.ind = T)
    if(length(idx)>0){
      eval$TPR = mean(TPR[-idx])
      eval$TDR = mean(TDR[-idx])
      eval$FPR = mean(FPR[-idx])
      eval$FDR = mean(FDR[-idx])
    } else{
      eval$TPR = mean(TPR)
      eval$TDR = mean(TDR)
      eval$FPR = mean(FPR)
      eval$FDR = mean(FDR)
    }
  }
  if(effect_size == 0.25){
    t_list = result
    n_sim = dim(t_list)[1]
    M = dim(t_list)[2]
    dimension = sqrt(M)
    
    TPR = numeric(n_sim)
    TDR = numeric(n_sim)
    FPR = numeric(n_sim)
    FDR = numeric(n_sim)
    
    for(i in 1:n_sim){
      
      t = t_list[i,]
      if(!is.na(sum(t))){
        t = matrix(t, nrow = dimension, ncol = dimension)
        positives = sum(t[3:(dimension/2 - 2), 3:(dimension/2 - 2)]) + sum(t[3:(dimension/2 - 2), (dimension/2 + 3):(dimension - 2)]) + 
          sum(t[(dimension/2 + 3):(dimension - 2), 3:(dimension/2 - 2)]) + sum(t[(dimension/2 + 3):(dimension - 2), (dimension/2 + 3):(dimension - 2)])
        
        TPR[i] = sum(t[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))]) / (((dimension/2 - 4)^2)*1)
        TDR[i] = sum(t[c((dimension/2 + 3):(dimension - 2)), c(3:(dimension/2 - 2))]) / positives
        FPR[i] = (sum(t[3:(dimension/2 - 2), 3:(dimension/2 - 2)]) + sum(t[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)]) + sum(t[(dimension/2 + 3):(dimension-2),(dimension/2 + 3):(dimension-2)])) / (((dimension/2 - 4)^2)*3)
        FDR[i] = (sum(t[3:(dimension/2 - 2), 3:(dimension/2 - 2)]) + sum(t[3:(dimension/2-2),(dimension/2 + 3):(dimension-2)]) + sum(t[(dimension/2 + 3):(dimension-2),(dimension/2 + 3):(dimension-2)])) / positives
      }
    }
    idx = which(TPR == 0, arr.ind = T)
    if(length(idx)>0){
      eval$TPR = mean(TPR[-idx])
      eval$TDR = mean(TDR[-idx])
      eval$FPR = mean(FPR[-idx])
      eval$FDR = mean(FDR[-idx])
    } else{
      eval$TPR = mean(TPR)
      eval$TDR = mean(TDR)
      eval$FPR = mean(FPR)
      eval$FDR = mean(FDR)
    }
  }
  
  return(eval)
}

bias_mae_mse_y_total = function(result1, result2, result3, result4){
  result = list()
  
  result$bias = result1$bias$total + result2$bias$total + result3$bias$total + result4$bias$total
  result$variance = result1$variance$total + result2$variance$total + result3$variance$total + result4$variance$total
  result$mse = result1$mse$total + result2$mse$total + result3$mse$total + result4$mse$total
  
  return(result)
}

#############
### BLESS ###
#############

N = N_list[1]
lambda = lambda_list[1]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BLESS, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

result_BLESS_MCAR = read_in_results_P2('BLESS_DPE_MCAR', X, v0, idx = c(1:100))

marginal_MCAR = result_BLESS_MCAR$marginal_list

BLESS_confusion_matrix1_MCAR1 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t1, n_sim = n_sim, effect_size = 0.5)
BLESS_confusion_matrix2_MCAR1 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t2, n_sim = n_sim, effect_size = 0.25)
rm(result_BLESS_MCAR)

N = N_list[1]
lambda = lambda_list[2]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BLESS, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

result_BLESS_MCAR = read_in_results_P2('BLESS_DPE_MCAR', X, v0, idx = c(1:100))

marginal_MCAR = result_BLESS_MCAR$marginal_list

BLESS_confusion_matrix1_MCAR2 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t1, n_sim = n_sim, effect_size = 0.5)
BLESS_confusion_matrix2_MCAR2 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t2, n_sim = n_sim, effect_size = 0.25)
rm(result_BLESS_MCAR)

N = N_list[1]
lambda = lambda_list[3]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BLESS, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

result_BLESS_MCAR = read_in_results_P2('BLESS_DPE_MCAR', X, v0, idx = c(1:100))

marginal_MCAR = result_BLESS_MCAR$marginal_list

BLESS_confusion_matrix1_MCAR3 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t1, n_sim = n_sim, effect_size = 0.5)
BLESS_confusion_matrix2_MCAR3 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t2, n_sim = n_sim, effect_size = 0.25)
rm(result_BLESS_MCAR)

N = N_list[2]
lambda = lambda_list[1]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BLESS, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

result_BLESS_MCAR = read_in_results_P2('BLESS_DPE_MCAR', X, v0, idx = c(1:100))

marginal_MCAR = result_BLESS_MCAR$marginal_list

BLESS_confusion_matrix1_MCAR4 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t1, n_sim = n_sim, effect_size = 0.5)
BLESS_confusion_matrix2_MCAR4 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t2, n_sim = n_sim, effect_size = 0.25)
rm(result_BLESS_MCAR)

N = N_list[2]
lambda = lambda_list[2]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))),2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BLESS, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))
result_BLESS_MCAR = read_in_results_P2('BLESS_DPE_MCAR', X, v0, idx = c(1:100))

marginal_MCAR = result_BLESS_MCAR$marginal_list

BLESS_confusion_matrix1_MCAR5 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t1, n_sim = n_sim, effect_size = 0.5)
BLESS_confusion_matrix2_MCAR5 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t2, n_sim = n_sim, effect_size = 0.25)
rm(result_BLESS_MCAR)

N = N_list[2]
lambda = lambda_list[3]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BLESS, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))
result_BLESS_MCAR = read_in_results_P2('BLESS_DPE_MCAR', X, v0, idx = c(1:100))

marginal_MCAR = result_BLESS_MCAR$marginal_list

BLESS_confusion_matrix1_MCAR6 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t1, n_sim = n_sim, effect_size = 0.5)
BLESS_confusion_matrix2_MCAR6 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t2, n_sim = n_sim, effect_size = 0.25)
rm(result_BLESS_MCAR)

N = N_list[3]
lambda = lambda_list[1]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BLESS, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))
result_BLESS_MCAR = read_in_results_P2('BLESS_DPE_MCAR', X, v0, idx = c(1:100))

marginal_MCAR = result_BLESS_MCAR$marginal_list

BLESS_confusion_matrix1_MCAR7 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t1, n_sim = n_sim, effect_size = 0.5)
BLESS_confusion_matrix2_MCAR7 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t2, n_sim = n_sim, effect_size = 0.25)
rm(result_BLESS_MCAR)

N = N_list[3]
lambda = lambda_list[2]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BLESS, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))
result_BLESS_MCAR = read_in_results_P2('BLESS_DPE_MCAR', X, v0, idx = c(1:100))

marginal_MCAR = result_BLESS_MCAR$marginal_list

BLESS_confusion_matrix1_MCAR8 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t1, n_sim = n_sim, effect_size = 0.5)
BLESS_confusion_matrix2_MCAR8 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t2, n_sim = n_sim, effect_size = 0.25)
rm(result_BLESS_MCAR)

N = N_list[3]
lambda = lambda_list[3]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BLESS, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))
result_BLESS_MCAR = read_in_results_P2('BLESS_DPE_MCAR', X, v0,idx = c(1:100))

marginal_MCAR = result_BLESS_MCAR$marginal_list

BLESS_confusion_matrix1_MCAR9 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t1, n_sim = n_sim, effect_size = 0.5)
BLESS_confusion_matrix2_MCAR9 = evaluate_TPR_TDR_FPR_FDR(result_BLESS_MCAR$t2, n_sim = n_sim, effect_size = 0.25)
rm(result_BLESS_MCAR)

##############
### BSGLMM ###
##############

N = N_list[1]
lambda = lambda_list[1]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BSGLMM, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_BSGLMM = read_in_results_P2(model = 'BSGLMM', X)

BSGLMM_confusion_matrix11 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix21 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2, n_sim, effect_size = 0.25)
BSGLMM_confusion_matrix1_corrected1 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1_corrected, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix2_corrected1 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2_corrected, n_sim, effect_size = 0.25)
rm(result_BSGLMM)

N = N_list[1]
lambda = lambda_list[2]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BSGLMM, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))
n_sim = 100
result_BSGLMM = read_in_results_P2(model = 'BSGLMM', X)

BSGLMM_confusion_matrix12 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix22 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2, n_sim, effect_size = 0.25)
BSGLMM_confusion_matrix1_corrected2 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1_corrected, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix2_corrected2 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2_corrected, n_sim, effect_size = 0.25)
rm(result_BSGLMM)

N = N_list[1]
lambda = lambda_list[3]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BSGLMM, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))
n_sim = 100
result_BSGLMM = read_in_results_P2(model = 'BSGLMM', X)

BSGLMM_confusion_matrix13 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix23 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2, n_sim, effect_size = 0.25)
BSGLMM_confusion_matrix1_corrected3 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1_corrected, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix2_corrected3 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2_corrected, n_sim, effect_size = 0.25)
rm(result_BSGLMM)

N = N_list[2]
lambda = lambda_list[1]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BSGLMM, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_BSGLMM = read_in_results_P2(model = 'BSGLMM', X)

BSGLMM_confusion_matrix14 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix24 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2, n_sim, effect_size = 0.25)
BSGLMM_confusion_matrix1_corrected4 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1_corrected, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix2_corrected4 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2_corrected, n_sim, effect_size = 0.25)
rm(result_BSGLMM)

N = N_list[2]
lambda = lambda_list[2]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BSGLMM, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_BSGLMM = read_in_results_P2(model = 'BSGLMM', X)

BSGLMM_confusion_matrix15 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix25 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2, n_sim, effect_size = 0.25)
BSGLMM_confusion_matrix1_corrected5 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1_corrected, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix2_corrected5 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2_corrected, n_sim, effect_size = 0.25)
rm(result_BSGLMM)

N = N_list[2]
lambda = lambda_list[3]
X = matrix(cbind(c(rep(1, N/2),rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BSGLMM, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_BSGLMM = read_in_results_P2(model = 'BSGLMM', X)

BSGLMM_confusion_matrix16 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix26 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2, n_sim, effect_size = 0.25)
BSGLMM_confusion_matrix1_corrected6 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1_corrected, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix2_corrected6 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2_corrected, n_sim, effect_size = 0.25)
rm(result_BSGLMM)

N = N_list[3]
lambda = lambda_list[1]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BSGLMM, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_BSGLMM = read_in_results_P2(model = 'BSGLMM', X)

BSGLMM_confusion_matrix17 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix27 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2, n_sim, effect_size = 0.25)
BSGLMM_confusion_matrix1_corrected7 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1_corrected, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix2_corrected7 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2_corrected, n_sim, effect_size = 0.25)
rm(result_BSGLMM)

N = N_list[3]
lambda = lambda_list[2]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BSGLMM, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_BSGLMM = read_in_results_P2(model = 'BSGLMM', X)

BSGLMM_confusion_matrix18 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix28 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2, n_sim, effect_size = 0.25)
BSGLMM_confusion_matrix1_corrected8 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1_corrected, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix2_corrected8 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2_corrected, n_sim, effect_size = 0.25)
rm(result_BSGLMM)

N = N_list[3]
lambda = lambda_list[3]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_BSGLMM, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_BSGLMM = read_in_results_P2(model = 'BSGLMM', X)

BSGLMM_confusion_matrix19 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix29 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2, n_sim, effect_size = 0.25)
BSGLMM_confusion_matrix1_corrected9 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t1_corrected, n_sim, effect_size = 0.5)
BSGLMM_confusion_matrix2_corrected9 = evaluate_TPR_TDR_FPR_FDR(result_BSGLMM$t2_corrected, n_sim, effect_size = 0.25)
rm(result_BSGLMM)

#############
### Firth ###
#############

N = N_list[1]
lambda = lambda_list[1]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_Firth, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_Firth = read_in_results_P2(model = 'Firth', X)

Firth_confusion_matrix11 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1, n_sim, effect_size = 0.5)
Firth_confusion_matrix21 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2, n_sim, effect_size = 0.25)
Firth_confusion_matrix1_uncorrected1 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1_uncorrected, n_sim, effect_size = 0.5)
Firth_confusion_matrix2_uncorrected1 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2_uncorrected, n_sim, effect_size = 0.25)
rm(result_Firth)

N = N_list[1]
lambda = lambda_list[2]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_Firth, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_Firth = read_in_results_P2(model = 'Firth', X)

Firth_confusion_matrix12 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1, n_sim, effect_size = 0.5)
Firth_confusion_matrix22 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2, n_sim, effect_size = 0.25)
Firth_confusion_matrix1_uncorrected2 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1_uncorrected, n_sim, effect_size = 0.5)
Firth_confusion_matrix2_uncorrected2 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2_uncorrected, n_sim, effect_size = 0.25)
rm(result_Firth)

N = N_list[1]
lambda = lambda_list[3]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_Firth, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_Firth = read_in_results_P2(model = 'Firth', X)

Firth_confusion_matrix13 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1, n_sim, effect_size = 0.5)
Firth_confusion_matrix23 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2, n_sim, effect_size = 0.25)
Firth_confusion_matrix1_uncorrected3 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1_uncorrected, n_sim, effect_size = 0.5)
Firth_confusion_matrix2_uncorrected3 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2_uncorrected, n_sim, effect_size = 0.25)
rm(result_Firth)

N = N_list[2]
lambda = lambda_list[1]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_Firth, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim=100
result_Firth = read_in_results_P2(model = 'Firth', X)

Firth_confusion_matrix14 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1, n_sim, effect_size = 0.5)
Firth_confusion_matrix24 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2, n_sim, effect_size = 0.25)
Firth_confusion_matrix1_uncorrected4 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1_uncorrected, n_sim, effect_size = 0.5)
Firth_confusion_matrix2_uncorrected4 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2_uncorrected, n_sim, effect_size = 0.25)
rm(result_Firth)

N = N_list[2]
lambda = lambda_list[2]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_Firth, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_Firth = read_in_results_P2(model = 'Firth', X)

Firth_confusion_matrix15 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1, n_sim, effect_size = 0.5)
Firth_confusion_matrix25 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2, n_sim, effect_size = 0.25)
Firth_confusion_matrix1_uncorrected5 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1_uncorrected, n_sim, effect_size = 0.5)
Firth_confusion_matrix2_uncorrected5 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2_uncorrected, n_sim, effect_size = 0.25)
rm(result_Firth)

N = N_list[2]
lambda = lambda_list[3]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_Firth, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_Firth = read_in_results_P2(model = 'Firth', X)

Firth_confusion_matrix16 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1, n_sim, effect_size = 0.5)
Firth_confusion_matrix26 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2, n_sim, effect_size = 0.25)
Firth_confusion_matrix1_uncorrected6 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1_uncorrected, n_sim, effect_size = 0.5)
Firth_confusion_matrix2_uncorrected6 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2_uncorrected, n_sim, effect_size = 0.25)
rm(result_Firth)

N = N_list[3]
lambda = lambda_list[1]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_Firth, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_Firth = read_in_results_P2(model = 'Firth', X)

Firth_confusion_matrix17 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1, n_sim, effect_size = 0.5)
Firth_confusion_matrix27 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2, n_sim, effect_size = 0.25)
Firth_confusion_matrix1_uncorrected7 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1_uncorrected, n_sim, effect_size = 0.5)
Firth_confusion_matrix2_uncorrected7 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2_uncorrected, n_sim, effect_size = 0.25)
rm(result_Firth)

N = N_list[3]
lambda = lambda_list[2]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_Firth, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_Firth = read_in_results_P2(model = 'Firth', X)

Firth_confusion_matrix18 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1, n_sim, effect_size = 0.5)
Firth_confusion_matrix28 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2, n_sim, effect_size = 0.25)
Firth_confusion_matrix1_uncorrected8 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1_uncorrected, n_sim, effect_size = 0.5)
Firth_confusion_matrix2_uncorrected8 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2_uncorrected, n_sim, effect_size = 0.25)
rm(result_Firth)

N = N_list[3]
lambda = lambda_list[3]
X = matrix(cbind(c(rep(1, N/2), rep(0, N/2)), rep(c(rep(0, (N/4)), rep(1, (N/4))), 2)), nrow = N, ncol = P)

setwd(sprintf("%sN%slambda%sP%s/", path_results_Firth, N, lambda, P))
load(paste0(path_truth, "lambda", lambda, "/truth.RData"))

n_sim = 100
result_Firth = read_in_results_P2(model = 'Firth', X)

Firth_confusion_matrix19 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1, n_sim, effect_size = 0.5)
Firth_confusion_matrix29 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2, n_sim, effect_size = 0.25)
Firth_confusion_matrix1_uncorrected9 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t1_uncorrected, n_sim, effect_size = 0.5)
Firth_confusion_matrix2_uncorrected9 = evaluate_TPR_TDR_FPR_FDR(result_Firth$t2_uncorrected, n_sim, effect_size = 0.25)
rm(result_Firth)

##############
### Tables ###
##############

DPE_confusion_table_function = function(result_bless1, result_bless2, result_bless3 ,
                                        result_bless4, result_bless5, result_bless6 ,
                                        result_bless7, result_bless8, result_bless9 ,
                                        result_bsglmm1, result_bsglmm2, result_bsglmm3, 
                                        result_bsglmm4, result_bsglmm5, result_bsglmm6, 
                                        result_bsglmm7, result_bsglmm8, result_bsglmm9, 
                                        result_bsglmm_corrected1, result_bsglmm_corrected2, result_bsglmm_corrected3, 
                                        result_bsglmm_corrected4, result_bsglmm_corrected5, result_bsglmm_corrected6, 
                                        result_bsglmm_corrected7, result_bsglmm_corrected8, result_bsglmm_corrected9, 
                                        result_firth_uncorrected1, result_firth_uncorrected2, result_firth_uncorrected3, 
                                        result_firth_uncorrected4, result_firth_uncorrected5, result_firth_uncorrected6, 
                                        result_firth_uncorrected7, result_firth_uncorrected8, result_firth_uncorrected9, 
                                        result_firth1, result_firth2, result_firth3,
                                        result_firth4, result_firth5, result_firth6,
                                        result_firth7, result_firth8, result_firth9){
  eval_param = matrix(NA, nrow = 30, ncol = 6)
  colnames(eval_param) = c('lambda=1', 'lambda=2', 'lambda=3', 
                           'lambda=1', 'lambda=2', 'lambda=3')
  rownames(eval_param) = c('BLESS', 'BSGLMM (uncorrected)', 'BSGLMM (corrected)', 'Firth (uncorrected)', 'Firth (corrected)',
                           'BLESS', 'BSGLMM (uncorrected)', 'BSGLMM (corrected)', 'Firth (uncorrected)', 'Firth (corrected)', 
                           'BLESS', 'BSGLMM (uncorrected)', 'BSGLMM (corrected)', 'Firth (uncorrected)', 'Firth (corrected)',
                           'BLESS', 'BSGLMM (uncorrected)', 'BSGLMM (corrected)', 'Firth (uncorrected)', 'Firth (corrected)',
                           'BLESS', 'BSGLMM (uncorrected)', 'BSGLMM (corrected)', 'Firth (uncorrected)', 'Firth (corrected)', 
                           'BLESS', 'BSGLMM (uncorrected)', 'BSGLMM (corrected)', 'Firth (uncorrected)', 'Firth (corrected)')
  
  ### TPR 
  evaluation = 'TPR'
  # N = 500
  eval_param[1, 1] = result_bless1[[evaluation]]
  eval_param[2, 1] = result_bsglmm1[[evaluation]]
  eval_param[3, 1] = result_bsglmm_corrected1[[evaluation]]
  eval_param[4, 1] = result_firth_uncorrected1[[evaluation]]
  eval_param[5, 1] = result_firth1[[evaluation]]
  
  eval_param[1, 2] = result_bless2[[evaluation]]
  eval_param[2, 2] = result_bsglmm2[[evaluation]]
  eval_param[3, 2] = result_bsglmm_corrected2[[evaluation]]
  eval_param[4, 2] = result_firth_uncorrected2[[evaluation]]
  eval_param[5, 2] = result_firth2[[evaluation]]
  
  eval_param[1, 3] = result_bless3[[evaluation]]
  eval_param[2, 3] = result_bsglmm3[[evaluation]]
  eval_param[3, 3] = result_bsglmm_corrected3[[evaluation]]
  eval_param[4, 3] = result_firth_uncorrected3[[evaluation]]
  eval_param[5, 3] = result_firth3[[evaluation]]
  # N = 1000
  eval_param[6, 1] = result_bless4[[evaluation]]
  eval_param[7, 1] = result_bsglmm4[[evaluation]]
  eval_param[8, 1] = result_bsglmm_corrected4[[evaluation]]
  eval_param[9, 1] = result_firth_uncorrected4[[evaluation]]
  eval_param[10, 1] = result_firth4[[evaluation]]
  
  eval_param[6, 2] = result_bless5[[evaluation]]
  eval_param[7, 2] = result_bsglmm5[[evaluation]]
  eval_param[8, 2] = result_bsglmm_corrected5[[evaluation]]
  eval_param[9, 2] = result_firth_uncorrected5[[evaluation]]
  eval_param[10, 2] = result_firth5[[evaluation]]
  
  eval_param[6, 3] = result_bless6[[evaluation]]
  eval_param[7, 3] = result_bsglmm6[[evaluation]]
  eval_param[8, 3] = result_bsglmm_corrected6[[evaluation]]
  eval_param[9, 3] = result_firth_uncorrected6[[evaluation]]
  eval_param[10, 3] = result_firth6[[evaluation]]
  # N = 5000
  eval_param[11, 1] = result_bless7[[evaluation]]
  eval_param[12, 1] = result_bsglmm7[[evaluation]]
  eval_param[13, 1] = result_bsglmm_corrected7[[evaluation]]
  eval_param[14, 1] = result_firth_uncorrected7[[evaluation]]
  eval_param[15, 1] = result_firth7[[evaluation]]
  
  eval_param[11, 2] = result_bless8[[evaluation]]
  eval_param[12, 2] = result_bsglmm8[[evaluation]]
  eval_param[13, 2] = result_bsglmm_corrected8[[evaluation]]
  eval_param[14, 2] = result_firth_uncorrected8[[evaluation]]
  eval_param[15, 2] = result_firth8[[evaluation]]
  
  eval_param[11, 3] = result_bless9[[evaluation]]
  eval_param[12, 3] = result_bsglmm9[[evaluation]]
  eval_param[13, 3] = result_bsglmm_corrected9[[evaluation]]
  eval_param[14, 3] = result_firth_uncorrected9[[evaluation]]
  eval_param[15, 3] = result_firth9[[evaluation]]
  
  ### TDR 
  evaluation = 'TDR'
  # N = 500
  eval_param[1, 4] = result_bless1[[evaluation]]
  eval_param[2, 4] = result_bsglmm1[[evaluation]]
  eval_param[3, 4] = result_bsglmm_corrected1[[evaluation]]
  eval_param[4, 4] = result_firth_uncorrected1[[evaluation]]
  eval_param[5, 4] = result_firth1[[evaluation]]
  
  eval_param[1, 5] = result_bless2[[evaluation]]
  eval_param[2, 5] = result_bsglmm2[[evaluation]]
  eval_param[3, 5] = result_bsglmm_corrected2[[evaluation]]
  eval_param[4, 5] = result_firth_uncorrected2[[evaluation]]
  eval_param[5, 5] = result_firth2[[evaluation]]
  
  eval_param[1, 6] = result_bless3[[evaluation]]
  eval_param[2, 6] = result_bsglmm3[[evaluation]]
  eval_param[3, 6] = result_bsglmm_corrected3[[evaluation]]
  eval_param[4, 6] = result_firth_uncorrected3[[evaluation]]
  eval_param[5, 6] = result_firth3[[evaluation]]
  # N = 1000
  eval_param[6, 4] = result_bless4[[evaluation]]
  eval_param[7, 4] = result_bsglmm4[[evaluation]]
  eval_param[8, 4] = result_bsglmm_corrected4[[evaluation]]
  eval_param[9, 4] = result_firth_uncorrected4[[evaluation]]
  eval_param[10, 4] = result_firth4[[evaluation]]
  
  eval_param[6, 5] = result_bless5[[evaluation]]
  eval_param[7, 5] = result_bsglmm5[[evaluation]]
  eval_param[8, 5] = result_bsglmm_corrected5[[evaluation]]
  eval_param[9, 5] = result_firth_uncorrected5[[evaluation]]
  eval_param[10, 5] = result_firth5[[evaluation]]
  
  eval_param[6, 6] = result_bless6[[evaluation]]
  eval_param[7, 6] = result_bsglmm6[[evaluation]]
  eval_param[8, 6] = result_bsglmm_corrected6[[evaluation]]
  eval_param[9, 6] = result_firth_uncorrected6[[evaluation]]
  eval_param[10, 6] = result_firth6[[evaluation]]
  # N = 5000
  eval_param[11, 4] = result_bless7[[evaluation]]
  eval_param[12,4] = result_bsglmm7[[evaluation]]
  eval_param[13,4] = result_bsglmm_corrected7[[evaluation]]
  eval_param[14,4] = result_firth_uncorrected7[[evaluation]]
  eval_param[15,4] = result_firth7[[evaluation]]
  
  eval_param[11,5] = result_bless8[[evaluation]]
  eval_param[12,5] = result_bsglmm8[[evaluation]]
  eval_param[13,5] = result_bsglmm_corrected8[[evaluation]]
  eval_param[14,5] = result_firth_uncorrected8[[evaluation]]
  eval_param[15,5] = result_firth8[[evaluation]]
  
  eval_param[11,6] = result_bless9[[evaluation]]
  eval_param[12,6] = result_bsglmm9[[evaluation]]
  eval_param[13,6] = result_bsglmm_corrected9[[evaluation]]
  eval_param[14,6] = result_firth_uncorrected9[[evaluation]]
  eval_param[15,6] = result_firth9[[evaluation]]
  
  ### FPR 
  evaluation = 'FPR'
  # N = 500
  eval_param[16, 1] = result_bless1[[evaluation]]
  eval_param[17, 1] = result_bsglmm1[[evaluation]]
  eval_param[18, 1] = result_bsglmm_corrected1[[evaluation]]
  eval_param[19, 1] = result_firth_uncorrected1[[evaluation]]
  eval_param[20, 1] = result_firth1[[evaluation]]
  
  eval_param[16, 2] = result_bless2[[evaluation]]
  eval_param[17, 2] = result_bsglmm2[[evaluation]]
  eval_param[18, 2] = result_bsglmm_corrected2[[evaluation]]
  eval_param[19, 2] = result_firth_uncorrected2[[evaluation]]
  eval_param[20, 2] = result_firth2[[evaluation]]
  
  eval_param[16, 3] = result_bless3[[evaluation]]
  eval_param[17, 3] = result_bsglmm3[[evaluation]]
  eval_param[18, 3] = result_bsglmm_corrected3[[evaluation]]
  eval_param[19, 3] = result_firth_uncorrected3[[evaluation]]
  eval_param[20, 3] = result_firth3[[evaluation]]
  # N = 1000
  eval_param[21, 1] = result_bless4[[evaluation]]
  eval_param[22, 1] = result_bsglmm4[[evaluation]]
  eval_param[23, 1] = result_bsglmm_corrected4[[evaluation]]
  eval_param[24, 1] = result_firth_uncorrected4[[evaluation]]
  eval_param[25, 1] = result_firth4[[evaluation]]
  
  eval_param[21, 2] = result_bless5[[evaluation]]
  eval_param[22, 2] = result_bsglmm5[[evaluation]]
  eval_param[23, 2] = result_bsglmm_corrected5[[evaluation]]
  eval_param[24, 2] = result_firth_uncorrected5[[evaluation]]
  eval_param[25, 2] = result_firth5[[evaluation]]
  
  eval_param[21, 3] = result_bless6[[evaluation]]
  eval_param[22, 3] = result_bsglmm6[[evaluation]]
  eval_param[23, 3] = result_bsglmm_corrected6[[evaluation]]
  eval_param[24, 3] = result_firth_uncorrected6[[evaluation]]
  eval_param[25, 3] = result_firth6[[evaluation]]
  # N = 5000
  eval_param[26, 1] = result_bless7[[evaluation]]
  eval_param[27, 1] = result_bsglmm7[[evaluation]]
  eval_param[28, 1] = result_bsglmm_corrected7[[evaluation]]
  eval_param[29, 1] = result_firth_uncorrected7[[evaluation]]
  eval_param[30, 1] = result_firth7[[evaluation]]
  
  eval_param[26, 2] = result_bless8[[evaluation]]
  eval_param[27, 2] = result_bsglmm8[[evaluation]]
  eval_param[28, 2] = result_bsglmm_corrected8[[evaluation]]
  eval_param[29, 2] = result_firth_uncorrected8[[evaluation]]
  eval_param[30, 2] = result_firth8[[evaluation]]
  
  eval_param[26, 3] = result_bless9[[evaluation]]
  eval_param[27, 3] = result_bsglmm9[[evaluation]]
  eval_param[28, 3] = result_bsglmm_corrected9[[evaluation]]
  eval_param[29, 3] = result_firth_uncorrected9[[evaluation]]
  eval_param[30, 3] = result_firth9[[evaluation]]
  
  ### FDR 
  evaluation = 'FDR'
  # N = 500
  eval_param[16, 4] = result_bless1[[evaluation]]
  eval_param[17, 4] = result_bsglmm1[[evaluation]]
  eval_param[18, 4] = result_bsglmm_corrected1[[evaluation]]
  eval_param[19, 4] = result_firth_uncorrected1[[evaluation]]
  eval_param[20, 4] = result_firth1[[evaluation]]
  
  eval_param[16, 5] = result_bless2[[evaluation]]
  eval_param[17, 5] = result_bsglmm2[[evaluation]]
  eval_param[18, 5] = result_bsglmm_corrected2[[evaluation]]
  eval_param[19, 5] = result_firth_uncorrected2[[evaluation]]
  eval_param[20, 5] = result_firth2[[evaluation]]
  
  eval_param[16, 6] = result_bless3[[evaluation]]
  eval_param[17, 6] = result_bsglmm3[[evaluation]]
  eval_param[18, 6] = result_bsglmm_corrected3[[evaluation]]
  eval_param[19, 6] = result_firth_uncorrected3[[evaluation]]
  eval_param[20, 6] = result_firth3[[evaluation]]
  # N = 1000
  eval_param[21, 4] = result_bless4[[evaluation]]
  eval_param[22, 4] = result_bsglmm4[[evaluation]]
  eval_param[23, 4] = result_bsglmm_corrected4[[evaluation]]
  eval_param[24, 4] = result_firth_uncorrected4[[evaluation]]
  eval_param[25, 4] = result_firth4[[evaluation]]
  
  eval_param[21, 5] = result_bless5[[evaluation]]
  eval_param[22, 5] = result_bsglmm5[[evaluation]]
  eval_param[23, 5] = result_bsglmm_corrected5[[evaluation]]
  eval_param[24, 5] = result_firth_uncorrected5[[evaluation]]
  eval_param[25, 5] = result_firth5[[evaluation]]
  
  eval_param[21, 6] = result_bless6[[evaluation]]
  eval_param[22, 6] = result_bsglmm6[[evaluation]]
  eval_param[23, 6] = result_bsglmm_corrected6[[evaluation]]
  eval_param[24, 6] = result_firth_uncorrected6[[evaluation]]
  eval_param[25, 6] = result_firth6[[evaluation]]
  # N = 5000
  eval_param[26, 4] = result_bless7[[evaluation]]
  eval_param[27,4] = result_bsglmm7[[evaluation]]
  eval_param[28,4] = result_bsglmm_corrected7[[evaluation]]
  eval_param[29,4] = result_firth_uncorrected7[[evaluation]]
  eval_param[30,4] = result_firth7[[evaluation]]
  
  eval_param[26,5] = result_bless8[[evaluation]]
  eval_param[27,5] = result_bsglmm8[[evaluation]]
  eval_param[28,5] = result_bsglmm_corrected8[[evaluation]]
  eval_param[29,5] = result_firth_uncorrected8[[evaluation]]
  eval_param[30,5] = result_firth8[[evaluation]]
  
  eval_param[26,6] = result_bless9[[evaluation]]
  eval_param[27,6] = result_bsglmm9[[evaluation]]
  eval_param[28,6] = result_bsglmm_corrected9[[evaluation]]
  eval_param[29,6] = result_firth_uncorrected9[[evaluation]]
  eval_param[30,6] = result_firth9[[evaluation]]
  
  write.table(eval_param, paste0(path_logs, 'BLESS_BSGLMM_Firth_TPR_TDR_FPR_FDR.csv'))
  
  xtable(eval_param, digits = 4)
}


print('Inference results for Beta1')
DPE_confusion_table_function(BLESS_confusion_matrix1_MCAR1,BLESS_confusion_matrix1_MCAR2, BLESS_confusion_matrix1_MCAR3,
                             BLESS_confusion_matrix1_MCAR4,BLESS_confusion_matrix1_MCAR5, BLESS_confusion_matrix1_MCAR6, 
                             BLESS_confusion_matrix1_MCAR7,BLESS_confusion_matrix1_MCAR8, BLESS_confusion_matrix1_MCAR9, 
                             BSGLMM_confusion_matrix11, BSGLMM_confusion_matrix12, BSGLMM_confusion_matrix13, 
                             BSGLMM_confusion_matrix14, BSGLMM_confusion_matrix15, BSGLMM_confusion_matrix16, 
                             BSGLMM_confusion_matrix17, BSGLMM_confusion_matrix18, BSGLMM_confusion_matrix19, 
                             BSGLMM_confusion_matrix1_corrected1, BSGLMM_confusion_matrix1_corrected2, BSGLMM_confusion_matrix1_corrected3, 
                             BSGLMM_confusion_matrix1_corrected4, BSGLMM_confusion_matrix1_corrected5, BSGLMM_confusion_matrix1_corrected6, 
                             BSGLMM_confusion_matrix1_corrected7, BSGLMM_confusion_matrix1_corrected8, BSGLMM_confusion_matrix1_corrected9, 
                             Firth_confusion_matrix1_uncorrected1,Firth_confusion_matrix1_uncorrected2, Firth_confusion_matrix1_uncorrected3, 
                             Firth_confusion_matrix1_uncorrected4,Firth_confusion_matrix1_uncorrected5, Firth_confusion_matrix1_uncorrected6, 
                             Firth_confusion_matrix1_uncorrected7,Firth_confusion_matrix1_uncorrected8, Firth_confusion_matrix1_uncorrected9, 
                             Firth_confusion_matrix11,Firth_confusion_matrix12, Firth_confusion_matrix13,
                             Firth_confusion_matrix14,Firth_confusion_matrix15, Firth_confusion_matrix16,
                             Firth_confusion_matrix17,Firth_confusion_matrix18, Firth_confusion_matrix19)

print('Inference results for Beta2')
DPE_confusion_table_function(BLESS_confusion_matrix2_MCAR1,BLESS_confusion_matrix2_MCAR2, BLESS_confusion_matrix2_MCAR3,
                             BLESS_confusion_matrix2_MCAR4,BLESS_confusion_matrix2_MCAR5, BLESS_confusion_matrix2_MCAR6, 
                             BLESS_confusion_matrix2_MCAR7,BLESS_confusion_matrix2_MCAR8, BLESS_confusion_matrix2_MCAR9, 
                             BSGLMM_confusion_matrix21, BSGLMM_confusion_matrix22, BSGLMM_confusion_matrix23, 
                             BSGLMM_confusion_matrix24, BSGLMM_confusion_matrix25, BSGLMM_confusion_matrix26, 
                             BSGLMM_confusion_matrix27, BSGLMM_confusion_matrix28, BSGLMM_confusion_matrix29, 
                             BSGLMM_confusion_matrix2_corrected1, BSGLMM_confusion_matrix2_corrected2, BSGLMM_confusion_matrix2_corrected3, 
                             BSGLMM_confusion_matrix2_corrected4, BSGLMM_confusion_matrix2_corrected5, BSGLMM_confusion_matrix2_corrected6, 
                             BSGLMM_confusion_matrix2_corrected7, BSGLMM_confusion_matrix2_corrected8, BSGLMM_confusion_matrix2_corrected9, 
                             Firth_confusion_matrix2_uncorrected1,Firth_confusion_matrix2_uncorrected2, Firth_confusion_matrix2_uncorrected3, 
                             Firth_confusion_matrix2_uncorrected4,Firth_confusion_matrix2_uncorrected5, Firth_confusion_matrix2_uncorrected6, 
                             Firth_confusion_matrix2_uncorrected7,Firth_confusion_matrix2_uncorrected8, Firth_confusion_matrix2_uncorrected9, 
                             Firth_confusion_matrix21,Firth_confusion_matrix22, Firth_confusion_matrix23,
                             Firth_confusion_matrix24,Firth_confusion_matrix25, Firth_confusion_matrix26,
                             Firth_confusion_matrix27,Firth_confusion_matrix28, Firth_confusion_matrix29)


