#!/bin/bash
#
#$ -N UKBB_cluster_size   #name for submission
#$ -q short.qc            # submit to short queues on cluster
#$ -t 1:1000              # number of bootstrap samples

### request maximum of 24 hours of compute time
#$ -l h_rt = 30:00:00
#$ -l s_rt = 30:00:00
#
### -cwd means work in the directory where submitted
#$ -cwd 
#
### -j option combines output and error messages
#$ -j y
#$ -o logs
#$ -e logs

echo "*********************************"
echo "SGE Job ID: $JOB_ID"
echo "SGE Job ID: $SGE_JOB_ID"
echo "Run on host: `hostname`"
echo "Operating system: `uname -s`"
echo "Username: `whoami`"
echo "Started at: `date`"
echo "*********************************"

# Setup
. /etc/profile
. ~/.bash_profile

module load fsl  # load FSLeyes software

# Switch < ... > for appropriate file path.
# path_t: path containing nifti files with standardized effects posterior mean / posterior std
# path_c_oindex: path for output (nifti files with clusters grouped by cluster index)
# path_c_osize: path for output (nifti files with clusters grouped by cluster size)
# path_c_txt: path for output (text file containing number of cluster & their cluster size)
# path_c_final: path for output of cluster size determined by posterior mean / posterior standard deviation from concatenated bootstrap replicates
# Cluster defining threshold of 2.3

# Option 1: Run to acquire cluster size output for every bootstrap sample (set -t 1:1000 because there are B = 1000 bootstrap samples)
cluster -i <path_t>t_${SGE_TASK_ID}.nii.gz -t 2.3 -o <path_c>oindex/cluster_index${SGE_TASK_ID} --osize = <path_c_osize>cluster_size${SGE_TASK_ID} > <path_c_txt>cluster_info${SGE_TASK_ID}.txt

# Option 2: Run to acquire cluster size determined by posterior mean / posterior standard deviation from concatenated bootstrap replicates (set -t 1).
cluster -i <path_t>t_final.nii.gz -t 2.3 -o <path_c_final>cluster_index.nii.gz --osize = <path_c_final>cluster_size.nii.gz > <path_c_final>cluster_info.txt

echo "*********************************"
echo "finished at: `date`"
echo "*********************************"
exit 0

