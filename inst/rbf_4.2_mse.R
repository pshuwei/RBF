library(dplyr)
library(mvtnorm)
library(grf)
source("../R/rbf.R")

#getting mimic dataset for testing
#dataset <- readr::read_csv("./test_data_locf.csv")[, -1]
data(dataset)

dataset2 <- dataset %>% group_by(subject_id) %>%
  filter(n() >= 12) %>% group_by(subject_id) %>%
  filter(any(discharge_location == "HOME")) %>% 
  #ilter(charttime <= 4) 
  #select(subject_id, charttime, sbp, bicar, na, age, weight, los_icu, death, sofa, vent_ind, vaso_ind) %>%
  select(subject_id, charttime, sbp, bicar, na, age, weight, k, map, los_icu, death, sofa, vent_ind, vaso_ind) %>%
  mutate(trt_group = case_when(any(vent_ind) & any(vaso_ind) ~ "1",
                               any(vent_ind) & is.na(vaso_ind) ~ "2",
                               is.na(vent_ind) & any(vaso_ind) ~ "3",
                               is.na(vent_ind) & is.na(vaso_ind) ~ NA
  )) %>% filter(!all(is.na(trt_group))) %>% 
  filter(charttime <= 12) 

dataset2$trt_group <- as.numeric(dataset2$trt_group)

dataset2[, 3:9] <- apply(dataset2[, 3:9], 2, function(x) log(x))

dataset3 <- dataset2 %>% group_by(subject_id) %>%
  summarize(across(sbp:map, list(mean = ~mean(., na.rm = T))))

dataset3[sapply(dataset3, is.nan)] <- NA

dataset4 <- dataset2 %>% group_by(subject_id) %>% arrange(charttime) %>%
  slice_tail(n = 1) %>% select(subject_id, los_icu, sofa, trt_group)

mimic_data <- dataset3 %>% left_join(dataset4)

imp <- mice::mice(mimic_data)
imp_base <- mice::complete(imp)
imp_base[,2:8] <- apply(imp_base[,2:8], 2, function(x) (x - min(x))/(max(x)-min(x)))

set.seed(1)

group1 <- length(which(imp_base$trt_group == 1))/nrow(imp_base)
group2 <- length(which(imp_base$trt_group == 2))/nrow(imp_base)
group3 <- length(which(imp_base$trt_group == 3))/nrow(imp_base)

n <- floor(2/3 * nrow(imp_base))

n_group1 <- ceiling(group1 * n)
n_group2 <- ceiling(group2 * n)
n_group3 <- n - n_group1 - n_group2

# Create data frames for each group
group1_sampled <- subset(imp_base, trt_group == 1) %>%
  sample_n(size = n_group1)
group2_sampled <- subset(imp_base, trt_group == 2) %>%
  sample_n(size = n_group2)
group3_sampled <- subset(imp_base, trt_group == 3) %>%
  sample_n(size = n_group3)

mimic_train <- rbind(group1_sampled, group2_sampled, group3_sampled)

mimic_test <- anti_join(imp_base, mimic_train)

x_train <- as.matrix(mimic_train[, 2:8])
x_test <- as.matrix(mimic_test[, 2:8])
x <- x_train

#treatment set up
n <- nrow(mimic_train)
z <- unlist(mimic_train$trt_group)
gs <- length(unique(z))

f10 <- 10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5])
f20 <- (5*x[,2])/(1+x[,1]^2)+ 5*sin(x[,3]*x[,4]) + x[,5] #2*x[,3]+(x[,4]) + x[,5]
f30 <- 0.1*exp(4*x[,1]) + 4/(1 + exp(-20*(x[,2] - 0.5))) + 3*x[,3]^2 + 2*x[,4] + x[,5]

f1 <- f10+f20#/2
f2 <- (f10+f20)/2+f30/3
f3 <- (f20+f30)/2 + f10/3
#f4 <- (5*x[,2])/(1+x[,1]^2)+3*x[,3]+2*x[,4]+x[,5]

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
      rbfest <- scheme_mult(y, x, x_test, z, gs = 3, sy, add=0, skip = 10, sigma_k = NULL, po=2, C = 16, Total_itr = 5000, burn = 2500)
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
      

    }
  #}
#}

apply(results, 2, median)

res <- rbind(apply(results, 2, median)[1:3+1],apply(results, 2, median)[1:3+4])
rownames(res) <- c("Shared-neuron RBF", "grf")

apply(res[-c(1:4+2),], 2, which.min)

