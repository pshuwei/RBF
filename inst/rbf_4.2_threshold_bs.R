library(dplyr)
library(mvtnorm)
library(grf)
source("../R/OthermethodsLow.R")
source("../R/rbf_blp.R")

# getting mimic dataset for testing
#dataset <- readr::read_csv("./test_data_locf.csv")[, -1]
data(dataset)

dataset2 <- dataset %>%
    group_by(subject_id) %>%
    filter(n() >= 12) %>%
    group_by(subject_id) %>%
    filter(any(discharge_location == "HOME")) %>%
    # ilter(charttime <= 4) select(subject_id, charttime, sbp, bicar,
    # na, age, weight, los_icu, death, sofa, vent_ind, vaso_ind) %>%
select(subject_id, charttime, sbp, bicar, na, age, weight, k, map, los_icu,
    sofa, vent_ind, vaso_ind) %>%
    mutate(trt_group = case_when(any(vent_ind) & any(vaso_ind) ~ "1", any(vent_ind) &
        is.na(vaso_ind) ~ "2", is.na(vent_ind) & any(vaso_ind) ~ "3", is.na(vent_ind) &
        is.na(vaso_ind) ~ NA)) %>%
    filter(!all(is.na(trt_group))) %>%
    filter(charttime <= 12)

dataset2$trt_group <- as.numeric(dataset2$trt_group)

dataset2[, 3:9] <- apply(dataset2[, 3:9], 2, function(x) log(x))

dataset3 <- dataset2 %>%
    group_by(subject_id) %>%
    summarize(across(sbp:map, list(mean = ~mean(., na.rm = T))))

dataset3[sapply(dataset3, is.nan)] <- NA

dataset4 <- dataset2 %>%
    group_by(subject_id) %>%
    arrange(charttime) %>%
    slice_tail(n = 1) %>%
    select(subject_id, los_icu, sofa, trt_group)

mimic_data <- dataset3 %>%
    left_join(dataset4)

imp <- mice::mice(mimic_data)
imp_base <- mice::complete(imp)
imp_base[, 2:8] <- apply(imp_base[, 2:8], 2, function(x) (x - min(x))/(max(x) -
    min(x)))

set.seed(1)

group1 <- length(which(imp_base$trt_group == 1))/nrow(imp_base)
group2 <- length(which(imp_base$trt_group == 2))/nrow(imp_base)
group3 <- length(which(imp_base$trt_group == 3))/nrow(imp_base)

# imp_base2 <- imp_base %>% sample_n(size = 171)

n <- nrow(imp_base)

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

# mimic_test <- anti_join(imp_base, mimic_train) results21 <-
# data.frame(#intercept = numeric(0), sbp = numeric(0), bicar =
# numeric(0), na = numeric(0), age = numeric(0), weight = numeric(0), k
# = numeric(0), map = numeric(0)) results31 <- data.frame(#intercept =
# numeric(0), sbp = numeric(0), bicar = numeric(0), na = numeric(0), age
# = numeric(0), weight = numeric(0), k = numeric(0), map = numeric(0))

results21_2 <- data.frame(intercept = numeric(0), sbp = numeric(0), bicar = numeric(0),
    na = numeric(0), age = numeric(0), weight = numeric(0), k = numeric(0),
    map = numeric(0))

results31_2 <- data.frame(intercept = numeric(0), sbp = numeric(0), bicar = numeric(0),
    na = numeric(0), age = numeric(0), weight = numeric(0), k = numeric(0),
    map = numeric(0))

results21_2ci <- data.frame(intercept = numeric(0), sbp = numeric(0), bicar = numeric(0),
    na = numeric(0), age = numeric(0), weight = numeric(0), k = numeric(0),
    map = numeric(0))

results31_2ci <- data.frame(intercept = numeric(0), sbp = numeric(0), bicar = numeric(0),
    na = numeric(0), age = numeric(0), weight = numeric(0), k = numeric(0),
    map = numeric(0))

results21_2coef <- data.frame(intercept = numeric(0), sbp = numeric(0), bicar = numeric(0),
    na = numeric(0), age = numeric(0), weight = numeric(0), k = numeric(0),
    map = numeric(0))

results31_2coef <- data.frame(intercept = numeric(0), sbp = numeric(0), bicar = numeric(0),
    na = numeric(0), age = numeric(0), weight = numeric(0), k = numeric(0),
    map = numeric(0))

sig21 <- data.frame(intercept = numeric(0), sbp = numeric(0), bicar = numeric(0),
    na = numeric(0), age = numeric(0), weight = numeric(0), k = numeric(0),
    map = numeric(0))

sig31 <- data.frame(intercept = numeric(0), sbp = numeric(0), bicar = numeric(0),
    na = numeric(0), age = numeric(0), weight = numeric(0), k = numeric(0),
    map = numeric(0))

rbf_ls <- list()
x_ls <- list()

reps <- 10

for (rep in 1:reps) {
    set.seed(rep)

    mimic_bs <- mimic_train %>%
        slice_sample(n = 172, replace = TRUE)

    n <- nrow(mimic_bs)

    x <- as.matrix(mimic_bs[, 2:8])
    x_ls[[rep]] <- x

    z <- unlist(mimic_bs$trt_group)
    gs <- length(unique(z))

    f10 <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + (10 *
        x[, 4] + 5 * x[, 5])
    f20 <- (5 * x[, 2])/(1 + x[, 1]^2) + 5 * sin(x[, 3] * x[, 4]) + x[, 5]  #2*x[,3]+(x[,4]) + x[,5]
    f30 <- 0.1 * exp(4 * x[, 1]) + 4/(1 + exp(-20 * (x[, 2] - 0.5))) + 3 *
        x[, 3]^2 + 2 * x[, 4] + x[, 5]

    f1 <- f10 + f20  #/2
    f2 <- (f10 + f20)/2 + f30/3
    f3 <- (f20 + f30)/2 + f10/3

    y1 <- rnorm(length(which(mimic_bs$trt_group == 1)), mean = f1, sd = 1)
    y2 <- rnorm(length(which(mimic_bs$trt_group == 2)), mean = f2, sd = 1)
    y3 <- rnorm(length(which(mimic_bs$trt_group == 3)), mean = f3, sd = 1)

    y <- c(y1, y2, y3)

    rbfest <- scheme_mult(y, x, x, z, gs = 3, sy, add = 0, skip = 10, sigma_k = NULL,
        po = 2, C = 16, Total_itr = 15000, burn = 5000)

    rbf_ls[[rep]] <- rbfest
}

save(rbf_ls, file = "rbf_bs_samps.Rdata")

for (rep in 1:reps) {

    samps <- 1000

    blp21 <- data.frame(numeric(0))
    blp31 <- data.frame(numeric(0))

    for (samp in 1:samps) {
        blp21 <- rbind(blp21, t(summary(lm(rbf_ls[[rep]][[1]][[samp]] ~ x_ls[[rep]]))$coefficients[,
            1]))
        blp31 <- rbind(blp31, t(summary(lm(rbf_ls[[rep]][[2]][[samp]] ~ x_ls[[rep]]))$coefficients[,
            1]))
    }

    # thresholding
    threshold <- 0

    blp21 <- apply(blp21, 2, function(x) ifelse(abs(x) < threshold, 0, x))
    blp31 <- apply(blp31, 2, function(x) ifelse(abs(x) < threshold, 0, x))


    for (j in 1:8) {
        results21_2coef[rep, ] <- colMeans(blp21)

        quantiles <- quantile(blp21[, j], c(0.025, 0.975))
        results21_2[rep, j] <- ifelse(quantiles[1] <= 0 & 0 <= quantiles[2],
            0, 1)
        results21_2ci[rep, j] <- paste("(", round(quantiles[1], 3), ",", round(quantiles[2],
            3), ")", sep = "")

        results31_2coef[rep, ] <- colMeans(blp31)

        quantiles <- quantile(blp31[, j], c(0.025, 0.975))
        # print(quantiles)
        results31_2[rep, j] <- ifelse(quantiles[1] <= 0 & 0 <= quantiles[2],
            0, 1)
        results31_2ci[rep, j] <- paste("(", round(quantiles[1], 3), ",", round(quantiles[2],
            3), ")", sep = "")
    }
}

# sig21 <- rbind(sig21, colSums(results21_2)) colnames(sig21) <-
# c('Intercept', 'SBP', 'Bicarbonate', 'Sodium', 'Age', 'Weight',
# 'Potassium', 'MAP') rownames(sig21) <- c('Threshold = 0.5', 'Threshold
# = 0.8', 'Threshold = 1.0', 'Threshold = 1.5', 'Threshold = 2.0') sig31
# <- rbind(sig31, colSums(results31_2)) colnames(sig31) <-
# c('Intercept', 'SBP', 'Bicarbonate', 'Sodium', 'Age', 'Weight',
# 'Potassium', 'MAP') rownames(sig31) <- c('Threshold = 0.5', 'Threshold
# = 0.8', 'Threshold = 1.0', 'Threshold = 1.5', 'Threshold = 2.0')

sig21 <- cbind(sig21, colSums(results21_2))
sig31 <- cbind(sig31, colSums(results31_2))
