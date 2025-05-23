#packages
library(mvtnorm)
library(dplyr)
library(grf)
source("../R/rbf.R")

set.seed(1)

n = 180 #for three groups

x <-matrix(runif(5*n),nrow=n)

#Creating test and train dataset
train_proportion <- 2/3
train_size <- round(train_proportion*n)
train_indices <- sample(1:n, train_size)

x_train <- x[train_indices,]
#x_test <- x[train_indices,] #for the estimated values of the training set
x_test <- x[-train_indices,] #for the predicted values of the testing set

n <- train_proportion*n

f10 <- 10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5])
f20 <- (5*x[,2])/(1+x[,1]^2)+ 5*sin(x[,3]*x[,4]) + x[,5] #2*x[,3]+(x[,4]) + x[,5]
f30 <- 0.1*exp(4*x[,1]) + 4/(1 + exp(-20*(x[,2] - 0.5))) + 3*x[,3]^2 + 2*x[,4] + x[,5]

f1 <- f10+f20#/2
f2 <- (f10+f20)/2+f30/3
f3 <- (f20+f30)/2 + f10/3
#f4 <- (5*x[,2])/(1+x[,1]^2)+3*x[,3]+2*x[,4]+x[,5]

#treatment set up
z <- rep(seq(1:3), each = n/3)
gs <- length(unique(z))

#Repetitions
reps <- 10

#Store results
results <- data.frame(
  rep = numeric(0),
  rbf21 = numeric(0),
  rbf31 = numeric(0),
  grf21 = numeric(0),
  grf31 = numeric(0)
)

#Start replications
for (rep in 1:reps) {
  print(rep)
  
  set.seed(rep)
  
  #Generating outcomes from Normal(f_g, 1)
  y1 <- rnorm(n/gs, mean = f1, sd = 1)
  y2 <- rnorm(n/gs, mean = f2, sd = 1)
  y3 <- rnorm(n/gs, mean = f3, sd = 1)
  
  #Store all outcomes as a vector
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
  grfest32 <- mean(((tau.hat[,2]-tau.hat[,1])- (f3_test-f2_test))^2)
  
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
## Range of errors for RBF v GRF
min_max <- apply(results, 2, 
                 function(x) {
                   paste0("(", round(min(x), 2), ", ", round(max(x), 2), ")", sep = "")
                 }
                 )[-1]

#Compile the results as presented in Simulation 4.1
tab_4.1_mse <- data.frame(rbind(res[1,], min_max[1:2], res[2,], min_max[3:4]))
colnames(tab_4.1_mse) <- c("tau_21 MSE", "tau_31 MSE")
rownames(tab_4.1_mse)[c(1,3)] <- c("Shared-neuron RBF", "Casual Forest")

View(tab_4.1_mse)