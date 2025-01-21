# Multi-treatment RBF-net

This is an R package for the treatment effect estimation using radial basis function (RBF)-nets with shared neurons.

## How to install the R package
devtools::install_github("pshuwei/RBF")

## What the R package contains:

The 'R' folder contains the coded methodology of our proposed method in Section 3 of the paper.

The 'inst' folder contains all simulation code that we discuss in the paper.

The `rbf_4.1_mse.R` and `rbf_4.2_mse.R` code contains the simulation code discussed in Section 4.1 and 4.2 regarding the simulation settings and MSE results.

The `rbf_4.2_blp.R` code contains the simulation code discussed in Section 4.2 regarding the second simulation and obtaining the Best Linear Projections (BLPs) for each predictor for each $\tau(x)$.

The `rbf_4.2_threshold.R` contains the simulation code discussed in Section 4.2 regarding the thresholding inference step to identify important predictors based on the BLP.

Finally `results.R` is the code to run all the above code as well as obtain the MIMIC cohort ids.

### MIMIC Data

Our data comes from the Medical Information Mart for Intensive Care III (MIMIC-III) data. The data are not publicly available due to privacy or ethical restrictions, but can be accessed freely by submitting an application to PhysioNet. 

The 'data' folder contains the `rbf_mimic_cohort.csv` file which contains the subject IDs and the corresponding treatment groups we assign them based on the paper.