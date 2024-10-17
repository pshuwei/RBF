library(dplyr)
library(mvtnorm)
library(grf)
source("../R/rbf.R")

# getting mimic dataset for testing
#load(file='../data/dataset.rda')
data(dataset)

dataset2 <- dataset %>%
  group_by(subject_id) %>%
  #needs to have more than 12 observations
  filter(n() >= 12) %>%
  #discharged to HOME
  filter(any(discharge_location == "HOME")) %>%
  #selecting patients, time, and predictors of interest
  select(subject_id, charttime, sbp, bicar, na, age, weight, k, map, 
         los_icu, sofa, vent_ind, vaso_ind) %>%
  #creating treatment groups based on whether or they had a vasopressor and/or ventilation
  mutate(trt_group = case_when(any(vent_ind) & any(vaso_ind) ~ "1", 
                               any(vent_ind) & is.na(vaso_ind) ~ "2", 
                               is.na(vent_ind) & any(vaso_ind) ~ "3", 
                               is.na(vent_ind) & is.na(vaso_ind) ~ NA)) %>%
  #remove those that did not use either
  filter(!all(is.na(trt_group))) %>%
  #only collect the first 12 hours
  filter(charttime <= 12)

dataset2$trt_group <- as.numeric(dataset2$trt_group)

#apply the log transformation on the predictors
dataset2[, 3:9] <- apply(dataset2[, 3:9], 2, function(x) log(x))

#collect their average predictor values during their 12 hours
dataset3 <- dataset2 %>%
  group_by(subject_id) %>%
  summarize(across(sbp:map, list(mean = ~mean(., na.rm = T))))

dataset3[sapply(dataset3, is.nan)] <- NA

#obtaining their 12th hour LOS ICU and SOFA scores
dataset4 <- dataset2 %>%
  group_by(subject_id) %>%
  arrange(charttime) %>%
  slice_tail(n = 1) %>%
  select(subject_id, los_icu, sofa, trt_group)

mimic_data <- dataset3 %>%
  left_join(dataset4)

#perform imputation on the data
imp <- mice::mice(mimic_data)
imp_base <- mice::complete(imp)

#apply the min-max transformation on the predictors
imp_base[, 2:8] <- apply(imp_base[, 2:8], 2, function(x) (x - min(x))/(max(x) - min(x)))

n <- nrow(imp_base)

mimic_bs <- imp_base

n <- nrow(mimic_bs)

#obtain matrix for predictors
x <- as.matrix(mimic_bs[, 2:8])

#get a vector for treatment assignements
z <- unlist(mimic_train$trt_group)
gs <- length(unique(z)) #tells us how many groups there are

#individual regression functions
f10 <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + (10 * x[,4] + 5 * x[, 5])
f20 <- (5 * x[, 2])/(1 + x[, 1]^2) + 5 * sin(x[, 3] * x[, 4]) + x[, 5]  #2*x[,3]+(x[,4]) + x[,5]
f30 <- 0.1 * exp(4 * x[, 1]) + 4/(1 + exp(-20 * (x[, 2] - 0.5))) + 3 * x[,3]^2 + 2 * x[, 4] + x[, 5]

#Creating treatment-outcome functions
f1 <- f10 + f20  #/2
f2 <- (f10 + f20)/2 + f30/3
f3 <- (f20 + f30)/2 + f10/3

#Create list to store RBF results and predictor matrix for each repetition
rbf_ls <- list()

reps <- 10

for (rep in 1:reps) {
  set.seed(rep)
  
  #Generating outcomes
  y1 <- rnorm(length(which(mimic_train$trt_group == 1)), mean = f1, sd = 1)
  y2 <- rnorm(length(which(mimic_train$trt_group == 2)), mean = f2, sd = 1)
  y3 <- rnorm(length(which(mimic_train$trt_group == 3)), mean = f3, sd = 1)
  
  y <- c(y1, y2, y3)
  
  #Getting RBF results
  rbfest <- scheme_mult(y, x, x, z, gs = 3, add = 0, skip = 10, sigma = NULL, sigma_k = NULL,C = 16, Total_itr = 15000, burn = 5000)
  
  #Storing the BLP results
  rbf_ls[[rep]] <- rbfest
}

cutoffs <- seq(0.1, 2, by = 0.1)

#Array to get how many repetitions out of total repetitions is a predictor significant (dimensions are number of cutoffs, number of replicated datasets, and predictors including the intercept)
cutoff21 <- array(0, dim=c(length(cutoffs),length(rbf_ls), ncol(x) + 1))
cutoff31 <- array(0, dim=c(length(cutoffs),length(rbf_ls), ncol(x) + 1))

#Getting the threshold results
for (i in 1:length(cutoffs)) {#For each cutoff
  
  for (rep in 1:reps) {
    samps <- 1000
    
    blp21 <- data.frame(numeric(0))
    blp31 <- data.frame(numeric(0))
    
    # Initial BLP results without cutoffs
    for (samp in 1:samps) {
      #BLP estimates for tau_21
      blp21 <- rbind(blp21, t(summary(lm(rbf_ls[[rep]][["tau21"]][[samp]]  ~ x))$coefficients[,1]))
      #BLP estimates for tau_31
      blp31 <- rbind(blp31, t(summary(lm(rbf_ls[[rep]][["tau31"]][[samp]]  ~ x))$coefficients[,1]))
    }
    
    #Thresholding
    threshold <- cutoffs[i]
    
    #Apply threshold to BLP coefficients
    blp21 <- apply(blp21, 2, function(x) ifelse(abs(x) < threshold, 0 , x))
    blp31 <- apply(blp31, 2, function(x) ifelse(abs(x) < threshold, 0 , x))
    
    for (j in 1:(ncol(x) + 1)) {
      #Getting results for inference of \tau_21
      #Proportion of BLP projections falling outside the interval (s_{p,t} values)
      cutoff21[i, rep, j] <- sum(abs(blp21[,j]) > threshold)/nrow(blp21)
      
      #Repeat for \tau_31
      cutoff31[i, rep, j]  <- sum(abs(blp31[,j]) > threshold)/nrow(blp31)
    }
  }
}

#Get the average s_{p,t} values across all replicated datasets
res21 <- apply(cutoff21, c(1, 3), mean)
res31 <- apply(cutoff31, c(1, 3), mean)

res21
res31