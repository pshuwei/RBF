#### MIMIC cohort data (for Section 4.2 and Section 5) ####

# We provide the Subject ID and Treatment Group based on the usage of only mechanical ventilation (group 1), only vasopressor (group 2), and both (group 3)

#if file is .rda format
load("rbf_mimic_cohort.rda")

#if file is .csv format
#data <- read.csv("rbf_mimic_cohort.csv")

#### Simulation Case 1 (Section 4.1) ####

source("rbf_4.1_mse.R")

#### Simulation Case 2 (Section 4.2) ####

# We borrow seven predictors from the MIMIC data based on Employing the marginal screening approach from \cite{xue2017robust} based on Length of Stay in the ICU and SOFA Score

# Predictors are age, weight, bicarbonate level, sodium level, systolic blood pressure (SBP), potassium level, and mean arterial pressure (MAP)

source("rbf_4.2_mse.R")
source("rbf_4.2_blp.R")
source("rbf_4.2_threshold.R")

#### MIMIC Analysis ####

#The set of predictors considered in our analysis are:
#age, weight, hemoglobin, white blood cells, blood urea nitrogen (BUN), potassium, sodium, bicarbonate, creatinine, platelet count, heart rate, systolic blood pressure (SBP), diastolic blood pressure (DBP), mean arterial pressure (MAP), and temperature.

#We perform the same BLP and thresholding approach for the MIMIC analysis