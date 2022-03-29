##############
### Script ###
##############

# Run code to acquire cluster size distributions and credible intervals of cluster size.
# (Code for UK Biobank application)

#################
### Libraries ###
#################

library(oro.nifti)
library(ggplot2)

#################
### Constants ###
#################

# Number of bootstrap samples.
B = 1000
# Set paths.
# path: folder with nifti files containing cluster formations and their respective cluster size 
# (FSL output: option 1).
path_cluster_BB = ''
# path: folder with nifti file containing observed cluster formations and their respective cluster size from 
# cluster size determined via posterior mean / posterior std (FSL output: option 2).
path_cluster_observed = ''
# Set cluster you want to look at: largest cluster 1, second largest 2 , ...
c_identifier = 1

# Step 1: Get indices of clusters within the observed quantity.
cluster_size_observed = readNIfTI(paste0(path_cluster_observed, 'cluster_size.nii.gz'))
cluster_index_observed = readNIfTI(paste0(path_cluster_observed, 'cluster_index.nii.gz'))
cluster_index_observed = as.vector(cluster_index_observed@.Data)
cluster_unique = unique(cluster_index_observed)
c = c(max(cluster_unique):1)[c_identifier]
idx_cluster = which(cluster_index_observed == c, arr.ind = T)

# Step 2: Find overlap of indices in each concatenated row of cluster sizes.
cluster_BB = matrix(0, B, length(idx_cluster))
for(b in 1:B){
  cluster_BB[b,] = as.vector(readNIfTI(paste0(path_cluster_BB, 'cluster_size', toString(b), '.nii.gz'))@.Data)[idx_cluster]
}

cluster_BB_dist = apply(cluster_BB,1,max)

# Cluster Size: posterior mean
(cluster_BB_mean = mean(cluster_BB_dist))
# Cluster Size: posterior sd
(cluster_BB_sd = sd(cluster_BB_dist))
# Cluster Size: 95%-credible interval
(c(quantile(cluster_BB_dist, 0.025), quantile(cluster_BB_dist, 0.975)))

# Create directory for plots for each cluster.
dir.create(paste0(path_cluster_observed, 'cluster_', toString(c_identifier)))
write.csv(cluster_BB, paste0(path_cluster_observed, 'cluster_', toString(c_identifier), '/cluster_BB.csv'))
write.csv(cluster_BB_dist, paste0(path_cluster_observed, 'cluster_', toString(c_identifier), '/cluster_BB_dist.csv'))

# Plot cluster size distribution with mean +/- std.
pdf(paste0(path_cluster_observed, 'cluster_', toString(c_identifier), '/cluster_size_histogram_mean_sd.pdf'))
 p = ggplot(data.frame(cluster_BB_dist), aes(x=cluster_BB_dist)) + 
  geom_histogram(binwidth = 15, fill="skyblue", color="#e9ecef", alpha=0.9) +
  xlim(cluster_BB_mean - 3*cluster_BB_sd, cluster_BB_mean + 3*cluster_BB_sd) +
  ggtitle("Histogram of Posterior Cluster Size") +
  theme(panel.background = element_rect(fill='white', colour='black'), plot.title = element_text(size=15)) + 
  xlab("Cluster Size") + 
  ylab("Count") 
 p = p + geom_vline(xintercept=cluster_BB_mean, color = 'blue')
 p = p + geom_vline(xintercept=cluster_BB_mean - cluster_BB_sd, color = 'blue', linetype='dashed')
 p = p + geom_vline(xintercept=cluster_BB_mean + cluster_BB_sd, color = 'blue', linetype='dashed')
 p
dev.off()

# Plot cluster size distribution with 95% credible interval.
cluster_size_value = unique(as.vector(cluster_size_observed@.Data)[idx_cluster])
pdf(paste0(path_cluster_observed, 'cluster_', toString(c_identifier), '/cluster_size_histogram_credible_interval.pdf'))
p = ggplot(data.frame(cluster_BB_dist), aes(x=cluster_BB_dist)) + 
  geom_histogram(binwidth = 15, fill="skyblue", color="#e9ecef", alpha=0.9) +
  xlim(cluster_BB_mean - 3*cluster_BB_sd, cluster_BB_mean + 3*cluster_BB_sd) +
  ggtitle("Histogram of Posterior Cluster Size") +
  theme(panel.background = element_rect(fill='white', colour='black'), plot.title = element_text(size=15)) + 
  xlab("Cluster Size") + 
  ylab("Count") 
p = p + geom_vline(xintercept=cluster_size_value, color = 'blue')
p = p + geom_vline(xintercept=quantile(cluster_BB_dist, 0.025), color = 'blue', linetype='dashed')
p = p + geom_vline(xintercept=quantile(cluster_BB_dist, 0.975), color = 'blue', linetype='dashed')
p
dev.off()
