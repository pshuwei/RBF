#packages
library(mvtnorm)
library(dplyr)
library(grf)
source("../R/OthermethodsLow.R")
source("../R/rbf.R")

set.seed(1)

n = 180 #for three groups
#n = 200 #for four groups

x <-matrix(runif(5*n),nrow=n)

#Creating test and train dataset
train_proportion <- 2/3
train_size <- round(train_proportion*n)
train_indices <- sample(1:n, train_size)

x_train <- x[train_indices,]
#x_test <- x[train_indices,] #for the estimated values of the training set
x_test <- x[-train_indices,] #for the predicted values of the testing set

x <- x_train
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

reps <- 10
pos <- 2
results <- data.frame(#po = pos[poi],
  rep = numeric(0),
  #add = adds[j],
  rbf21 = numeric(0),
  rbf31 = numeric(0),
  rbf32 = numeric(0),
  grf21 = numeric(0),
  grf31 = numeric(0),
  grf32 = numeric(0)
)

#reps <- 1
mat <- array(0, dim=c(reps,7, 3))

#rownames(mat) <- c("rvmerror", "rferror", "svmerror", "OLSerror", "BART20", "BART50", "BART200")
#for (poi in 1:length(pos)) {
for (rep in 1:reps) {
  #for (j in 1:length(adds)) {
  #set.seed(123)
  #print(pos[poi])
  print(rep)
  #print(adds[j])
  
  set.seed(rep)
  
  #results[30*(poi-1) + rep,1] <- pos[poi]
  #results[30*(poi-1) + rep,2] <- rep
  
  y1 <- rnorm(n/gs, mean = f1, sd = 1)
  y2 <- rnorm(n/gs, mean = f2, sd = 1)
  y3 <- rnorm(n/gs, mean = f3, sd = 1)
  #y4 <- rnorm(n/gs, mean = f4, sd = 2)
  
  y <- c(y1, y2, y3)
  
  
  #results[30*(poi-1) + rep, 3] <- adds[j]
  
  x <- x_train
  rbfest <- scheme_mult(y, x, x_test, z, gs = 3, sy, add=0, skip = 10, sigma_k = NULL, po=2, C = 16, Total_itr = 15000, burn = 5000)
  x <- x_test
  
  #out <- BART::mbart(x_train, z, x_test, printevery=10000L)
  #emat <- matrix(colMeans(out$yhat.test), nrow = 3)
  
  #rbf mse est
  
  # f1_test <- 10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5])
  # f2_test <- f1_test/2 + (5*x[,2])+(1+x[,1]^2)+ 2*x[,3]+exp(x[,4]) + x[,5]
  # f3_test <- f2_test/3+0.1*exp(4*x[,1]) + 4/(1 + exp(-20*(x[,2] - 0.5))) + 3*x[,3]^2 + 2*x[,4] + x[,5]      #f4_test <- (5*x[,2])/(1+x[,1]^2)+3*x[,3]+2*x[,4]+x[,5]
  
  f10 <- 10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5])
  f20 <- (5*x[,2])/(1+x[,1]^2)+ 5*sin(x[,3]*x[,4]) + x[,5] #2*x[,3]+(x[,4]) + x[,5]
  f30 <- 0.1*exp(4*x[,1]) + 4/(1 + exp(-20*(x[,2] - 0.5))) + 3*x[,3]^2 + 2*x[,4] + x[,5]
  
  f1_test <- f10+f20
  f2_test <- (f10+f20)/2+f30/3
  f3_test <- (f20+f30)/2 + f10/3
  
  rbfest21 <- mean(((rbfest$f2est-rbfest$f1est) - (f2_test-f1_test))^2)
  rbfest31 <- mean(((rbfest$f3est-rbfest$f1est) - (f3_test-f1_test))^2)
  rbfest32 <- mean(((rbfest$f3est-rbfest$f2est) - (f3_test-f2_test))^2)
  
  #grf mse est
  mc.forest <- multi_arm_causal_forest(x_train, y, as.factor(z))
  mc.pred <- predict(object = mc.forest, newdata = x_test, drop = TRUE)
  
  tau.hat <- matrix(mc.pred$predictions, ncol=2)
  
  grfest21 <- mean((tau.hat[,1] - (f2_test-f1_test))^2)
  grfest31 <- mean((tau.hat[,2] - (f3_test-f1_test))^2)
  grfest32 <- mean(((tau.hat[,2]-tau.hat[,1])- (f3_test-f2_test))^2)
  
  results <- rbind(results, 
                   data.frame(#po = pos[poi],
                     rep = rep,
                     #add = adds[j],
                     rbf21 = rbfest21,
                     rbf31 = rbfest31,
                     rbf32 = rbfest32,
                     grf21 = grfest21,
                     grf31 = grfest31,
                     grf32 = grfest32
                   ))
  
  
  y1 <- y
  #Treatment 1
  ind <- which(z==1)
  out1 <- OthermethodLow(x_train[ind,], y1[ind], x_test)
  
  #Treatment 2
  ind <- which(z==2)
  out2 <- OthermethodLow(x_train[ind,], y1[ind], x_test)
  
  #Treatment 3
  ind <- which(z==3)
  out3 <- OthermethodLow(x_train[ind,], y1[ind], x_test)
  
  ###############Calculate the errors for different treatment contrasts on test set########################################
  
  mat[rep,, 1] <- colMeans((out1-out2 - (f1_test - f2_test))^2) #rvmerror, rferror, svmerror, OLSerror, BART20, BART50, BART200
  mat[rep,, 2] <- colMeans((out3-out1 - (f3_test - f1_test))^2) #rvmerror, rferror, svmerror, OLSerror, BART20, BART50, BART200
  mat[rep,, 3] <- colMeans((out2-out3 - (f2_test - f3_test))^2) #rvmerror, rferror, svmerror, OLSerror, BART20, BART50, BART200
  
  print(mat[rep,,])
  print(results[rep, ])
  #results[30*(poi-1) + rep,4] <- 
  #results[30*(poi-1) + rep,5] <- 
  #results[30*(poi-1) + rep,6] <- 
  #results[rep,5] <- mean((rbfest$f4est - f4_test)^2)
}
#}
#}

#colMeans(results)
apply(results, 2, median)
apply(mat, 2:3, median)

res <- rbind(apply(results, 2, median)[1:3+1],apply(results, 2, median)[1:3+4],
             apply(mat, 2:3, median))
rownames(res) <- c("Shared-neuron RBF", "grf", "rvmerror", "rferror", "svmerror", "OLSerror", "BART20", "BART50", "BART200")

apply(res[-c(1:4+2),], 2, which.min)

