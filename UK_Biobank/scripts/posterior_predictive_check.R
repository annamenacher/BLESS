# Libraries
library(ggplot2)
library(data.table)
library(oro.nifti)
library(ggpointdensity)

set.seed(1234)

# Sample size of held-out dataset generated (reduce to 1000 in paper)
N = 1500

# Define paths.
# Path for in sample subject data.
path_in_sample = ''
# Path for all sample subject data.
path_all = ''
# Path for lesion masks.
path = ""
# Path for white matter mask.
path_mask = ""
# Path for output for out of sample data.
path_out = ""
# Path with UKBB results for parameter estimates -> posterior draws to draw from posterior predictive distribution.
path_results = ""


# Read in data.
X = fread(path_all, 'X.csv')[,2:6]
id_out_sample = as.numeric(as.vector(unlist(X[,1])))
X_N2000 = fread(path_in_sample, 'X.csv')[,2:6]

# Get ID's for in-sample data vs out-of-sample held out data.
id_in_sample = as.numeric(as.vector(unlist(X_N2000[,1])))
id_out_sample = id_out_sample[-which(id_out_sample %in% id_in_sample)]
id = sample(1:length(id_out_sample), N, replace = F)

X = data.matrix(X)
X = X[id,]

X_N2000_mean_age = mean(as.numeric(as.vector(unlist(X_N2000[,3]))))
X_N2000_sd_age = sd(as.numeric(as.vector(unlist(X_N2000[,3]))))

X_N2000_mean_headsize = mean(as.numeric(as.vector(unlist(X_N2000[,5]))))
X_N2000_sd_headsize = sd(as.numeric(as.vector(unlist(X_N2000[,5]))))

age = (X[,3] - X_N2000_mean_age) / X_N2000_sd_age
sex = as.numeric(as.vector(unlist(X[,2])))
headsize = (X[,4] - X_N2000_mean_headsize) / X_N2000_sd_headsize
agebysex = age * sex

X_standardized = data.frame(matrix(0,N,4))
X_standardized[,1] = age
X_standardized[,2] = sex
X_standardized[,3] = headsize
X_standardized[,4] = agebysex
colnames(X_standardized) = c('age', 'sex', 'headsize', 'agebysex')

# Mask based on empirical lesion rate
mask = readNIfTI(paste0(path_mask, "Emp_map_2mm_mask.nii.gz"))
mask = as.vector(mask@.Data)
ind = which(mask == 1, arr.ind = T)
write.csv(mask, paste0(path_mask, "mask.csv"))
write.csv(ind, paste0(path_mask, "ind_mask.csv"))
M = length(ind)

# Lesion data
Y = matrix(NA, nrow = N, ncol = M)
for(i in 1:N){
  print(i)
  tryCatch({
    id = X[i,'eid']
    img = readNIfTI(paste0(path, toString(id), "_lesion_to_MNI_2mm.nii.gz"))
    Y[i,] = as.vector(img@.Data)[ind]
  }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

ind_not_avail = which(is.na(Y[,1]), arr.ind=T)

X = X[-ind_not_avail,]
X_standardized = X_standardized[-ind_not_avail,]
Y = Y[-ind_not_avail,]
X = X[1:1000,]
X_standardized = X_standardized[1:1000,]
Y = Y[1:1000,]

N = 1000

write.csv(Y, paste0(path_out, "Y_mask_N", toString(N), ".csv"))
write.csv(X_standardized, paste0(path_out, "X_standardised_N", toString(N), ".csv"))
write.csv(X, paste0(path_out, "X_N", toString(N), ".csv"))

beta0 = data.matrix(data.table::fread(paste(path_results, 'beta0.csv'))[,2:54729])
beta0 = apply(beta0,2,mean)
beta0 = as.numeric(as.vector(unlist(beta0)))
Beta1 = data.matrix(data.table::fread(paste(path_results, 'Beta1.csv'))[,2:54729])
Beta1 = apply(Beta1,2,mean)
Beta1 = as.numeric(as.vector(unlist(Beta1)))
Beta2 = data.matrix(data.table::fread(paste(path_results, 'Beta2.csv'))[,2:54729])
Beta2 = apply(Beta2,2,mean)
Beta2 = as.numeric(as.vector(unlist(Beta2)))
Beta3 = data.matrix(data.table::fread(paste(path_results, 'Beta3.csv'))[,2:54729])
Beta3 = apply(Beta3,2,mean)
Beta3 = as.numeric(as.vector(unlist(Beta3)))
Beta4 = data.matrix(data.table::fread(paste(path_results, 'Beta4.csv'))[,2:54729])
Beta4 = apply(Beta4,2,mean)
Beta4 = as.numeric(as.vector(unlist(Beta4)))
Beta = rbind(as.vector(beta0), as.vector(Beta1), as.vector(Beta2), as.vector(Beta3), as.vector(Beta4))
rm(beta0)
rm(Beta1)
rm(Beta2)
rm(Beta3)
rm(Beta4)

X = data.matrix(data.table::fread(paste0(path_out, 'X_standardised_N1000.csv')[,2:5])

y_pred = cbind(rep(1, 1000), X)%*%Beta

p = pnorm(as.vector(y_pred))
p = matrix(p, 1000, 54728)

write.csv(p, sprintf('%sout_of_sample_lesion_probability_map.csv', path_out))

y_pred = matrix(rbinom(1000*54728, 1, as.vector(p)), 1000, 54728)

write.csv(y_pred, sprintf('%sout_of_sample_binary_predicted_lesion_map.csv', path_out))

# Calibration plots
# specific voxel: 4907, 4963, 6287, 32484
j = 4907
  
x1 = data.matrix(fread(paste0(path_out, 'out_of_sample_lesion_probability_map.csv'))[,2:54729])
x2 = data.matrix(fread(paste0(path_out, 'Y_mask_N1000.csv'))[,2:54729])
x1 = x1[,j]
x2 = x2[,j]
x1 = as.numeric(as.vector(unlist(x1)))
x2 = as.numeric(as.vector(unlist(x2)))

bin.x = numeric(10)

bin.x[1] = mean(x1[which(x1 < 0.1, arr.ind = T)]) 
bin.x[2] = mean(x1[which((x1 > 0.1) & (x1 < 0.2), arr.ind = T)]) 
bin.x[3] = mean(x1[which((x1 > 0.2) & (x1 < 0.3), arr.ind = T)]) 
bin.x[4] = mean(x1[which((x1 > 0.3) & (x1 < 0.4), arr.ind = T)]) 
bin.x[5] = mean(x1[which((x1 > 0.4) & (x1 < 0.5), arr.ind = T)]) 
bin.x[6] = mean(x1[which((x1 > 0.5) & (x1 < 0.6), arr.ind = T)]) 
bin.x[7] = mean(x1[which((x1 > 0.6) & (x1 < 0.7), arr.ind = T)]) 
bin.x[8] = mean(x1[which((x1 > 0.7) & (x1 < 0.8), arr.ind = T)]) 
bin.x[9] = mean(x1[which((x1 > 0.8) & (x1 < 0.9), arr.ind = T)]) 
bin.x[10] = mean(x1[which(x1 > 0.9, arr.ind = T)]) 

bin.y = numeric(10)

bin.y[1] = sum(x2[which(x1 < 0.1, arr.ind = T)]) / length(which(x1 < 0.1, arr.ind = T))
bin.y[2] = sum(x2[which((x1 > 0.1) & (x1 < 0.2), arr.ind = T)]) / length(which((x1 > 0.1) & (x1 < 0.2), arr.ind = T))
bin.y[3] = sum(x2[which((x1 > 0.2) & (x1 < 0.3), arr.ind = T)]) / length(which((x1 > 0.2) & (x1 < 0.3), arr.ind = T))
bin.y[4] = sum(x2[which((x1 > 0.3) & (x1 < 0.4), arr.ind = T)]) / length(which((x1 > 0.3) & (x1 < 0.4), arr.ind = T))
bin.y[5] = sum(x2[which((x1 > 0.4) & (x1 < 0.5), arr.ind = T)]) / length(which((x1 > 0.4) & (x1 < 0.5), arr.ind = T))
bin.y[6] = sum(x2[which((x1 > 0.5) & (x1 < 0.6), arr.ind = T)]) / length(which((x1 > 0.5) & (x1 < 0.6), arr.ind = T))
bin.y[7] = sum(x2[which((x1 > 0.6) & (x1 < 0.7), arr.ind = T)]) / length(which((x1 > 0.6) & (x1 < 0.7), arr.ind = T))
bin.y[8] = sum(x2[which((x1 > 0.7) & (x1 < 0.8), arr.ind = T)]) / length(which((x1 > 0.7) & (x1 < 0.8), arr.ind = T))
bin.y[9] = sum(x2[which((x1 > 0.8) & (x1 < 0.9), arr.ind = T)]) / length(which((x1 > 0.8) & (x1 < 0.9), arr.ind = T))
bin.y[10] = sum(x2[which(x1 > 0.9, arr.ind = T)]) / length(which(x1 > 0.9, arr.ind = T))

pdat <- data.frame(y = x2, prob = x1, bin.x = bin.x, bin.y = bin.y)

ggplot(pdat, aes(prob, y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_line(aes(bin.x, bin.y), col = 'blue') +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  xlab("Mean Predicted Probability") +
  ylab("Fraction of Positives") +
  ggtitle(paste0("Calibration Plot: Voxel ", toString(j))) + 
  theme_classic() + 
  theme(text = element_text(size = 20))   + 
  xlim(0, 1) + 
  ylim(0, 1)
ggsave(paste0(path_out, "calibration_plot_wo_data_", toString(j), ".pdf"))

ggplot(pdat, aes(prob)) +
  geom_histogram(color = "black", fill = "white") +
  xlab("Mean Predicted Probability") +
  ylab("Count") +
  ggtitle(paste0("Histogram: Voxel ", toString(j)))+ 
  theme_classic() + 
  theme(text = element_text(size = 20))    + 
  xlim(0, 1) 

ggsave(paste0(path_out, "calibration_histogram_count_", toString(j), ".pdf"))

# Calibaration plot for all voxels in mask for all out of sample subjects

x1 = data.matrix(fread(path_out, 'out_of_sample_lesion_probability_map.csv')[,2:54729])
x2 = data.matrix(fread(path_out, 'Y_mask_N1000.csv')[,2:54729])

x1 = as.numeric(as.vector(unlist(x1)))
x2 = as.numeric(as.vector(unlist(x2)))

bin.x = numeric(10)

bin.x[1] = mean(x1[which(x1 < 0.1, arr.ind = T)]) 
bin.x[2] = mean(x1[which((x1 > 0.1) & (x1 < 0.2), arr.ind = T)]) 
bin.x[3] = mean(x1[which((x1 > 0.2) & (x1 < 0.3), arr.ind = T)]) 
bin.x[4] = mean(x1[which((x1 > 0.3) & (x1 < 0.4), arr.ind = T)]) 
bin.x[5] = mean(x1[which((x1 > 0.4) & (x1 < 0.5), arr.ind = T)]) 
bin.x[6] = mean(x1[which((x1 > 0.5) & (x1 < 0.6), arr.ind = T)]) 
bin.x[7] = mean(x1[which((x1 > 0.6) & (x1 < 0.7), arr.ind = T)]) 
bin.x[8] = mean(x1[which((x1 > 0.7) & (x1 < 0.8), arr.ind = T)]) 
bin.x[9] = mean(x1[which((x1 > 0.8) & (x1 < 0.9), arr.ind = T)]) 
bin.x[10] = mean(x1[which(x1 > 0.9, arr.ind = T)]) 

bin.y = numeric(10)

bin.y[1] = sum(x2[which(x1 < 0.1, arr.ind = T)]) / length(which(x1 < 0.1, arr.ind = T))
bin.y[2] = sum(x2[which((x1 > 0.1) & (x1 < 0.2), arr.ind = T)]) / length(which((x1 > 0.1) & (x1 < 0.2), arr.ind = T))
bin.y[3] = sum(x2[which((x1 > 0.2) & (x1 < 0.3), arr.ind = T)]) / length(which((x1 > 0.2) & (x1 < 0.3), arr.ind = T))
bin.y[4] = sum(x2[which((x1 > 0.3) & (x1 < 0.4), arr.ind = T)]) / length(which((x1 > 0.3) & (x1 < 0.4), arr.ind = T))
bin.y[5] = sum(x2[which((x1 > 0.4) & (x1 < 0.5), arr.ind = T)]) / length(which((x1 > 0.4) & (x1 < 0.5), arr.ind = T))
bin.y[6] = sum(x2[which((x1 > 0.5) & (x1 < 0.6), arr.ind = T)]) / length(which((x1 > 0.5) & (x1 < 0.6), arr.ind = T))
bin.y[7] = sum(x2[which((x1 > 0.6) & (x1 < 0.7), arr.ind = T)]) / length(which((x1 > 0.6) & (x1 < 0.7), arr.ind = T))
bin.y[8] = sum(x2[which((x1 > 0.7) & (x1 < 0.8), arr.ind = T)]) / length(which((x1 > 0.7) & (x1 < 0.8), arr.ind = T))
bin.y[9] = sum(x2[which((x1 > 0.8) & (x1 < 0.9), arr.ind = T)]) / length(which((x1 > 0.8) & (x1 < 0.9), arr.ind = T))
bin.y[10] = sum(x2[which(x1 > 0.9, arr.ind = T)]) / length(which(x1 > 0.9, arr.ind = T))

pdat <- data.frame(bin.x = bin.x, bin.y = bin.y)

ggplot(pdat) +
  geom_abline(slope = 1, intercept = 0) +
  geom_line(aes(bin.x, bin.y), col = 'blue') +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  xlab("Mean Predicted Probability") +
  ylab("Fraction of Positives") +
  ggtitle(paste0("Calibration Plot")) + 
  theme_classic() + 
  theme(text = element_text(size = 20)) + 
  xlim(0, 1) + 
  ylim(0, 1)
ggsave(paste0(path_out, "calibration_plot_all_voxels.pdf"))

rm(x2)
pdat <- data.frame(prob = x1)
rm(x1)

ggplot(pdat, aes(prob)) +
  geom_histogram(color = "black", fill = "white") +
  xlab("Mean Predicted Probability") +
  ylab("Count") +
  ggtitle(paste0("Histogram of Predicted Probability")) + 
  theme_classic() + 
  theme(text = element_text(size = 20))   + 
  xlim(0, 1) 
ggsave(paste0(path_out, "calibration_histogram_count_all_voxels.pdf"))

# Empirical lesion probabilites comparison of true lesion rates and predicted lesion rate estimated via posterior predictive distribution on held-out data.
x1 = data.matrix(fread(path_out, 'out_of_sample_lesion_probability_map.csv')[,2:54729])
x2 = data.matrix(fread(path_out, 'Y_mask_N1000.csv')[,2:54729])

x1 = apply(x1, 2, mean)
x2 = apply(x2, 2, mean)
pdat <- data.frame(y_pred = x1, y_true = x2)

ggplot(pdat, aes(x = y_pred, y = y_true)) +
  geom_pointdensity(adjust = .1, show.legend = F) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = 'white'),
        legend.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  geom_abline(intercept = 0) + 
  xlab('Predicted lesion rate') + 
  ylab('Empiricial lesion rate') +
  theme(text = element_text(size = 35)) +
  guides(fill = guide_legend(title = ""))
ggsave(paste0(path_out, "lesion_pred_lesion_empir_scatterplot.pdf"))


