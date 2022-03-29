# Bayesian Lesion Estimaton with a Structured Spike-and-Slab Prior

## Data

### Abstract
- Simulation study: 2D lesion masks generated via data generating process described in Supplementary Materials Sec. 5
- UK Biobank: 3D lesion masks

X: NxP matrix: subject-specific covariates
Y: NxM matrix: lesion masks

### Availability 
- Simulation study: Data can be generated via simulation study script 'data_generation.R'. We additionally provide an example simulation study containing 100 datasets (= lesion masks) generated with the setting of sample size N=1,000 and base rate intensity lambda = 3 to showcase the reproducibility of our method. 
- UK Biobank: Healthcare data requires access permission and cannot be attached here (https://www.ukbiobank.ac.uk/).

### Description
Simulated 2D lesion masks provided in this repositry under simulation_study/data/N1000lambda3P2 and 3D lesion masks provided from the UK Biobank (not uploaded due to access restrictions).

## Code

### Abstract 
Code to run the mass-univariate method Firth regression, which fits a model at each voxel independently, the Bayesian spatial model BSGLMM, and our method BLESS & BB- BLESS. Additionally, we provide all code to generate figures & acquire cluster size distributions via FSL in a shell script.

### Description 

Simualtion study:
- BB_BLESS.R
- BLESS.R
- BLESS_CAR.R
- BLESS_Gibbs.R
- BLESS_PxCAR.R
- BSGLMM.R
- data_generation.R
- firth_regression.R
- generate_truth.R
- marginal_plot.R
- plotting_results_bb_bless_simulation_study.R
- regularisation_plot.R

UK Biobank:
- A.R
- BB_BLESS_3D.R
- BLESS_3D.R
- cluster_size_distribtion.R
- combine_BB_BLESS.R
- firth_regression.R
- fsl_script.sh
- plotting_UKBB.R

### Optional Information 
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

## Instructions for Use

### Reproducibility
- Simulation study: example case provided for N=1,000 and lambda=3 (all other studies can also be reproduced by generating additional datasets with 'data_generation.R')

### Replication (Optional)
