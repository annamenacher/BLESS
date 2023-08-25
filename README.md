# Bayesian Lesion Estimaton with a Structured Spike-and-Slab Prior

![alt text](https://github.com/annamenacher/BLESS/blob/main/bb_bless_maps.png?raw=true)

## Data

### Abstract
- Simulation study: 2D lesion masks generated via data generating process described in Supplementary Materials Sec. 5
- UK Biobank: 3D lesion masks

- X: NxP matrix: subject-specific covariates
- Y: NxM matrix: lesion masks

### Availability 
- Simulation study: Data can be generated via simulation study script 'data_generation.R'. We additionally provide an example simulation study containing 100 datasets (= lesion masks) generated with the setting of sample size N=1,000 and base rate intensity lambda = 3 to showcase the reproducibility of our method. 
- UK Biobank: Healthcare data requires access permission and cannot be attached here (https://www.ukbiobank.ac.uk/).

### Description
Simulated 2D lesion masks provided in this repositry under simulation_study/data/N1000lambda3P2 and 3D lesion masks provided from the UK Biobank (not uploaded due to access restrictions).

## Code

### Abstract 
Code to run the mass-univariate method Firth regression, which fits a model at each voxel independently, the Bayesian spatial model BSGLMM, and our method BLESS-VI & BB- BLESS. Additionally, we provide all code to generate figures & acquire cluster size distributions via FSL in a shell script.

### Description 

Simualtion study:
- **data_generation.R** : Generates simulation study data.
- **generate_truth.R** : Generates true parameter values and lesion rates.
- **firth_regression.R**: Performs parameter estimation and inference for indepdent mass-univariate Firth regression models at each voxel location. 
- **BSGLMM.R** : Performs parameter estimation and inference via BSGLMM (= Bayesian Spatial Generalized Linear Mixed Model).
- **BLESS.R** : Performs parameter estimation and inference via BLESS-VI (= Bayesian Lesion Estimation via Structured Spike-and-slab Prior), specifically deploys dynamic posterior exploration.
- **BB_BLESS.R** : Performs parameter estimation and inference via BB-BLESS (= Bayesian Bootstrap - BLESS) with output from DPE from BLESS.R as initialization.
- **BLESS_Gibbs.R** : Performs parameter estimation and inference via BLESS-Gibbs (BLESS estimated via Gibbs sampler, gold standard method).
- **BLESS_CAR.R** : Performs parameter estimation and inference via BLESS (with spatial prior: CAR).
- **BLESS_PxCAR.R** Performs parameter estimation and inference via BLESS (with spatial prior: PxCAR).
- **marginal_plot.R** : Plots marginal posterior of gamma under prior v0=0 over sequence of spike variances.
- **regularisation_plot.R** : Plots regularization plot of parameter estiamtes over sequence of spike variances.
- **plotting_results_bb_bless_simulation_study.R** : Plotting script for additional figures appearing in the manuscript (i.e. comparison of marginal posterior densities for a single active/inactive voxel, comparison of KL-divergence and Wasserstein metric across marginal densities for parameters of all voxels, ...). 
- **table_parameter_prediction_inference.R** : Ouputs all the tables of the main paper and supplementary material to summarize parameter estimation, prediction and inference results in latex table format.
- **plot_bar_plot_tpr_tdr_fpr_fdr.R** : Plot bar plots of simulation study to compare inference results (TPR, TDR, FPR, FDR) for all methods (BLESS, BSGLMM, Firth) for various sample sizes (N=500, N=1,000, N=5,000) and base rate intensities (lambda=1, lambda=2, lambda=3).
  
UK Biobank:
- **A.R** : Create list of neighboring indices for every voxel in a masked 3D lattice. 
- firth_regression.R: Performs parameter estimation and inference for indepdent mass-univariate Firth regression models at each voxel location. 
- **BLESS_3D.R** : Performs parameter estimation and inference via BLESS-VI (= Bayesian Lesion Estimation via Structured Spike-and-slab Prior), specifically deploys dynamic posterior exploration.
- **BB_BLESS_3D.R** : Performs parameter estimation and inference via BB-BLESS (= Bayesian Bootstrap - BLESS) with output from DPE from BLESS.R as initialization.
- **combine_BB_BLESS.R** : Combine the bootstrap samples from BB-BLESS (which were executed in parallel) for further analysis and calculation of cluster size distributions.
- **fsl_script.sh** : Bash script to execute FSL command -cluster (if the calculation of cluster size based imaging statistics are of interest).
- **cluster_size_distribtion.R** : Plotting of cluster size distribution with credible intervals of cluster size.
- **plotting_UKBB.R** : Plotting script for UK Biobank application of 2D axial slices (Code kindly provided by [Petya Kindalova](https://github.com/petyakindalova)).
- **posterior_predictive_check.R** : Perform posterior predictive checking and plot figures from Section 4.5 of the supplementary material. 

### Optional Information 

**R libraries:**
The following R packages are necessary to successfully use the codes:

- deldir
- spatstat
- MASS
- hutils
- scales
- logistf
- dplyr
- reshape2
- brglm2
- LaplacesDemon
- logOfGamma
- latex2exp
- oro.nifti

For cluster size inference: 
- Software package: FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes)

**Creation of figures and tables in the main paper and supplementary material:**

Paper:
- Fig. 1: FSL
- Fig. 2: a) regularisation_plot.R b) marginal_plot.R
- Fig. 3: plotting_results_bb_bless_simulation_study.R
- Fig. 4: plot_bar_plot_tpr_tdr_fpr_fdr.R
- Fig. 5: plotting_UKBB.R
- Fig. 6: plotting_UKBB.R
- Tab. 1,2: table_parameter_prediction_inference.R

Supplementary material:
- Fig. 1: Latex
- Fig. 2: FSL
- Fig. 3: a) regularisation_plot.R b) marginal_plot.R
- Fig. 4: plotting_UKBB.R
- Fig. 5: posterior_predictive_check.R
- Fig. 6: plotting_UKBB.R
- Fig. 7: posterior_predictive_check.R
- Fig. 8: plotting_UKBB.R
- Fig. 9: regularisation_plot.R
- Fig. 10: data_generation.R
- Fig. 11: marginal_plot.R
- Fig. 12, 13, 14, 15: plotting_results_bb_bless_simulation_study.R
- Fig. 16, 17, 18, 19: FSL
- Tab. 1-25: table_parameter_prediction_inference.R

## Instructions for Use

### Reproducibility
- Simulation study: example case provided for N=1,000 and lambda=3 (all other studies can also be reproduced by generating additional datasets with 'data_generation.R')
- Files are listed in order of execution for one single dataset. If one desires to replicate the entire simulation study, then execution of scripts for all datasets {1, ... , 100} is required. It is recommended to run Firth regression first as it is needed for initialization for BLESS. More so, BB-BLESS also requires the output from BLESS as initialization. 
- UK Biobank application: Code to run UK Biobank application with masked lesion data. We concatenate the masked nifti files, containing lesion masks, into a NxM matrix, where N represents sample size and M is the number of voxels in the mask, and save the output as a csv file. 

### Replication (Optional)
