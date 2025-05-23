library(dplyr)
library(mvtnorm)
library(grf)
source("../R/rbf.R")

mimic_train <- data

#obtain matrix for predictors
x <- as.matrix(mimic_train[, 2:8])

#get a vector for treatment assignements
z <- unlist(mimic_train$trt_group)
gs <- length(unique(z)) #tells us how many groups there are

#individual regression functions
f10 <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + (10 * x[,
    4] + 5 * x[, 5])
f20 <- (5 * x[, 2])/(1 + x[, 1]^2) + 5 * sin(x[, 3] * x[, 4]) + x[, 5]  #2*x[,3]+(x[,4]) + x[,5]
f30 <- 0.1 * exp(4 * x[, 1]) + 4/(1 + exp(-20 * (x[, 2] - 0.5))) + 3 * x[,
    3]^2 + 2 * x[, 4] + x[, 5]

#Creating treatment-outcome functions
f1 <- f10 + f20  #/2
f2 <- (f10 + f20)/2 + f30/3
f3 <- (f20 + f30)/2 + f10/3


#Store RBF results for multiple repetitions
rbf_ls <- list()

#Create dataframes for use of BLP estimation
results21ci <- results31ci <- results21coef <- results31coef <- 
  data.frame(intercept = numeric(0), sbp = numeric(0), bicar = numeric(0),
    na = numeric(0), age = numeric(0), weight = numeric(0), k = numeric(0),
    map = numeric(0))

reps <- 1

for (rep in 1:reps) {
  set.seed(rep)
  
  #Generating outcomes
  y1 <- rnorm(length(which(mimic_train$trt_group == 1)), mean = f1, sd = 1)
  y2 <- rnorm(length(which(mimic_train$trt_group == 2)), mean = f2, sd = 1)
  y3 <- rnorm(length(which(mimic_train$trt_group == 3)), mean = f3, sd = 1)

  y <- c(y1, y2, y3)

  #Getting RBF results
  rbfest <- scheme_mult(y, x, x, z, gs = 3, add = 0, skip = 10, sigma = NULL, sigma_k = NULL,C = 16, Total_itr = 15000, burn = 5000)

  rbf_ls[[rep]] <- rbfest

}

#BLP Process
for (rep in 1:reps) {
    samps <- 1000
    
    #dataframes to store BLP estimates
    blp21 <- data.frame(numeric(0))
    blp31 <- data.frame(numeric(0))
  
    for (samp in 1:samps) {
        #BLP estimates for tau_21
        blp21 <- rbind(blp21, t(summary(lm(rbf_ls[[rep]][["tau21"]][[samp]] ~ x))$coefficients[,
            1]))
        #BLP estimates for tau_31
        blp31 <- rbind(blp31, t(summary(lm(rbf_ls[[rep]][["tau31"]][[samp]] ~ x))$coefficients[,
            1]))
    }

    for (j in 1:8) {
        #Getting results for inference of \tau_21
        #Average BLP coefficient
        results21coef[rep, ] <- colMeans(blp21)
        #Get the posterior credible intervals
        quantiles <- quantile(blp21[, j], c(0.025, 0.975))
        #Obtain the actual credible intervals
        results21ci[rep, j] <- paste("(", round(quantiles[1], 3), ", ", round(quantiles[2],
            3), ")", sep = "")
        
        #Repeat for \tau_31
        results31coef[rep, ] <- colMeans(blp31)
        quantiles <- quantile(blp31[, j], c(0.025, 0.975))
        results31ci[rep, j] <- paste("(", round(quantiles[1], 3), ", ", round(quantiles[2],
            3), ")", sep = "")
    }
}

tab4.2_blp <- cbind(t(rbind(round(results21coef, 2), results21ci)),
                    t(rbind(round(results31coef,2), results31ci)))
colnames(tab4.2_blp) <- c("tau_21 BLP", "tau_21 95% CI", "tau_31 BLP", "tau_31 95% CI")

View(tab4.2_blp)
