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

mimic_data <- dataset3 %>% left_join(dataset4)

#perform imputation on the data
imp <- mice::mice(mimic_data)
imp_base <- mice::complete(imp)

#apply the min-max transformation on the predictors
imp_base[, 2:8] <- apply(imp_base[, 2:8], 2, function(x) (x - min(x))/(max(x) - min(x)))

set.seed(1)

#Creating the training set from the 172 patients
#Getting the proportion of each treatment group
group1 <- length(which(imp_base$trt_group == 1))/nrow(imp_base)
group2 <- length(which(imp_base$trt_group == 2))/nrow(imp_base)
group3 <- length(which(imp_base$trt_group == 3))/nrow(imp_base)

#Set this n to be the training set size
n <- floor(2/3 * nrow(imp_base))

#Set the size for each treatment group for the training set.
n_group1 <- ceiling(group1 * n)
n_group2 <- ceiling(group2 * n)
n_group3 <- n - n_group1 - n_group2

# Create data frames for each training treatment group
group1_sampled <- subset(imp_base, trt_group == 1) %>%
  sample_n(size = n_group1)
group2_sampled <- subset(imp_base, trt_group == 2) %>%
  sample_n(size = n_group2)
group3_sampled <- subset(imp_base, trt_group == 3) %>%
  sample_n(size = n_group3)

#Combine all to make the training data set
mimic_train <- rbind(group1_sampled, group2_sampled, group3_sampled)

#Set the difference between the original and the training to be the test set
mimic_test <- anti_join(imp_base, mimic_train)

#Creating matrix for both training and testing predictors

x_train <- as.matrix(mimic_train[, 2:8])
x_test <- as.matrix(mimic_test[, 2:8])
x <- x_train

#treatment set up
n <- nrow(mimic_train)

#get a vector for treatment assignements
z <- unlist(mimic_train$trt_group)
gs <- length(unique(z)) #tells us how many groups there are

#individual regression functions
f10 <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + (10 * x[,4] + 5 * x[, 5])
f20 <- (5 * x[, 2])/(1 + x[, 1]^2) + 5 * sin(x[, 3] * x[, 4]) + x[, 5]  #2*x[,3]+(x[,4]) + x[,5]
f30 <- 0.1 * exp(4 * x[, 1]) + 4/(1 + exp(-20 * (x[, 2] - 0.5))) + 3 * x[,
                                                                         3]^2 + 2 * x[, 4] + x[, 5]

#Creating treatment-outcome functions
f1 <- f10 + f20  #/2
f2 <- (f10 + f20)/2 + f30/3
f3 <- (f20 + f30)/2 + f10/3

reps <- 10
results <- data.frame(
  rep = numeric(0),
  rbf21 = numeric(0),
  rbf31 = numeric(0),
  grf21 = numeric(0),
  grf31 = numeric(0)
)

for (rep in 1:reps) {
 
    print(rep)
    set.seed(rep)
    
    #Generate outcome for each group based on Normal(f_g, 1)
    y1 <- rnorm(n/gs, mean = f1, sd = 1)
    y2 <- rnorm(n/gs, mean = f2, sd = 1)
    y3 <- rnorm(n/gs, mean = f3, sd = 1)
    
    #Store the results
    y <- c(y1, y2, y3)

    x <- x_train
    
    #Run the function on the training set and get the estimated treatment outcomes \hat{f_g} for the testing set
    rbfest <- scheme_mult(y, x, x_test, z, gs = 3, add = 0, skip = 10, sigma_k = NULL, C = 16, Total_itr = 15000, burn = 5000)
    
    x <- x_test

    #individual regression functions
    f10 <- 10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5])
    f20 <- (5*x[,2])/(1+x[,1]^2)+ 5*sin(x[,3]*x[,4]) + x[,5] #2*x[,3]+(x[,4]) + x[,5]
    f30 <- 0.1*exp(4*x[,1]) + 4/(1 + exp(-20*(x[,2] - 0.5))) + 3*x[,3]^2 + 2*x[,4] + x[,5]
    
    #Get the "actual" treatment outcomes of the test set
    f1_test <- f10+f20
    f2_test <- (f10+f20)/2+f30/3
    f3_test <- (f20+f30)/2 + f10/3
    
    #Get MSE of the actual CATE estimator versus the estimated CATE estimator from our RBF method
    f1hat <- Reduce('+', rbfest$f1est)/length(rbfest$f1est)
    f2hat <- Reduce('+', rbfest$f2est)/length(rbfest$f2est)
    f3hat <- Reduce('+', rbfest$f3est)/length(rbfest$f3est)
    
    rbfest21 <- mean(((f2hat-f1hat) - (f2_test-f1_test))^2)
    rbfest31 <- mean(((f3hat-f1hat) - (f3_test-f1_test))^2)
    
    #TO do the same process for the Causal Random Forest
    ## Training the causal forest
    mc.forest <- multi_arm_causal_forest(x_train, y, as.factor(z))
    ## Testing the causal forest
    mc.pred <- predict(object = mc.forest, newdata = x_test, drop = TRUE)
    
    ##Obtaining the estimated CATE estimators 
    tau.hat <- matrix(mc.pred$predictions, ncol=2)
    
    #Get MSE of the actual CATE estimator versus the estimated CATE estimator from the casual forest method
    grfest21 <- mean((tau.hat[,1] - (f2_test-f1_test))^2)
    grfest31 <- mean((tau.hat[,2] - (f3_test-f1_test))^2)
    
    results <- rbind(results, 
                     data.frame(
                       rep = rep,
                       rbf21 = rbfest21,
                       rbf31 = rbfest31,
                       grf21 = grfest21,
                       grf31 = grfest31
                     ))
    
}

#Obtaining the median of MSE of the our RBF versus the Causal Forest
apply(results, 2, median)

#Store the results
res <- rbind(apply(results, 2, median)[1:2+1],apply(results, 2, median)[1:2+3])
rownames(res) <- c("Shared-neuron RBF", "Casual Forest")
res <- round(res, 2)

#Gives us the lowest MSE between the 2
apply(res[-c(1:4+2),], 2, which.min)

#To get the range of errors
min_max <- apply(results, 2, 
                 function(x) {
                   paste0("(", round(min(x), 2), ", ", round(max(x), 2), ")", sep = "")
                 }
)[-1]

# Compile the results as presented in Simulation 4.2
tab_4.2_mse <- data.frame(rbind(res[1,], min_max[1:2], res[2,], min_max[3:4]))
colnames(tab_4.2_mse) <- c("tau_21 MSE", "tau_31 MSE")
rownames(tab_4.2_mse)[c(1,3)] <- c("Shared-neuron RBF", "Casual Forest")

View(tab_4.2_mse)