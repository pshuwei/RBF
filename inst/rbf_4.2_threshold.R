library(dplyr)
library(mvtnorm)
library(grf)
source("../R/rbf.R")

mimic_bs <- data

n <- nrow(mimic_bs)

#obtain matrix for predictors
x <- as.matrix(mimic_bs[, 2:8])

#get a vector for treatment assignements
z <- unlist(mimic_bs$trt_group)
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
  y1 <- rnorm(length(which(mimic_bs$trt_group == 1)), mean = f1, sd = 1)
  y2 <- rnorm(length(which(mimic_bs$trt_group == 2)), mean = f2, sd = 1)
  y3 <- rnorm(length(which(mimic_bs$trt_group == 3)), mean = f3, sd = 1)
  
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

#Compile the results as presented in Simulation 4.2
row_names <- cutoffs

rownames(res21) <- rownames(res31) <- row_names
colnames(res21) <- colnames(res31) <- c("Intercept", "SBP", "Bicarbonate", "Sodium", "Age", "Weight", "Potassium", "MAP")

tab_4.2_threshold <- rbind(res31, res21)
rownames(tab_4.2_threshold) <- c(cutoffs, cutoffs)

View(tab_4.2_threshold)
