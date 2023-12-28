OthermethodLow <- function(Xtr, Ytr, Xte, Yte, Bayes = T) {
    library(randomForest)
    library(kernlab)
    library(e1071)
    library(glmnet)
    library(BART)

    Xtrs <- Xtr
    Xtes <- Xte
    # RVM
    rvmm <- NULL
    rvmerror <- NA
    out <- try(rvmm <- rvm(Xtr, Ytr, kernel = "rbfdot", type = "regression",
        list(sigma = 0.1)), silent = T)  #kpar='automatic', iterations  =100000)#
    multi <- 1
    while (class(out) == "try-error") {
        multi <- multi + 1
        out <- try(rvmm <- rvm(Xtr, Ytr, kernel = "rbfdot", type = "regression",
            list(sigma = 0.1 * multi)), silent = T)
        if (multi > 10) {
            break
        }
    }
    if (!is.null(rvmm)) {
        ytest <- predict(rvmm, Xte)
        rvmerror <- ytest  #-Yte)^2) 
    }

    # RF
    foo <- randomForest(x = data.frame(Xtr), y = Ytr, importance = TRUE)  #,
    ytest <- predict(foo, Xte)
    rferror <- ytest  #-Yte)^2)

    # SVM
    out <- svm(Xtr, Ytr, scale = TRUE, type = NULL, kernel = "radial")
    ytest <- predict(out, Xte)
    svmerror <- ytest  #-Yte)^2)

    # OLS
    fit <- lm(Ytr ~ Xtr - 1)
    beta <- array(fit$coefficients)
    Yp1 <- Xte %*% beta
    OLSerror <- Yp1  #-Yte)^2)

    if (Bayes == T) {
        post = BART::gbart(Xtr, Ytr, Xte, sparse = F, ntree = 20, nskip = 5000,
            ndpost = 5000, keepevery = 1, printevery = 10000L)
        ytest <- colMeans(post$yhat.test)
        resultB20 <- ytest

        out <- post  #list(ytrhat = post$yhat.train, ytehat = post$yhat.test, sigmahat = post$sigma)
        fitB20 <- out

        post = BART::gbart(Xtr, Ytr, Xte, sparse = F, ntree = 50, nskip = 5000,
            ndpost = 5000, keepevery = 1, printevery = 10000L)
        ytest <- colMeans(post$yhat.test)
        resultB50 <- ytest  #-Yte)^2)

        out <- post  #list(ytrhat = post$yhat.train, ytehat = post$yhat.test, sigmahat = post$sigma)
        fitB50 <- out

        post = BART::gbart(Xtr, Ytr, Xte, sparse = F, ntree = 200, nskip = 5000,
            ndpost = 5000, keepevery = 1, printevery = 10000L)
        ytest <- colMeans(post$yhat.test)
        resultB200 <- ytest  #-Yte)^2)

        out <- post  #list(ytrhat = post$yhat.train, ytehat = post$yhat.test, sigmahat = post$sigma)
        fitB200 <- out
    }

    finalout <- cbind(rvmerror, rferror, svmerror, OLSerror)
    if (Bayes == T) {
        finalout <- cbind(finalout, resultB20, resultB50, resultB200)
    }
    return(finalout)
}

# ###########data generation################## n = 180 #for three groups
# #n = 200 #for four groups x <-matrix(runif(5*n),nrow=n) #Creating test
# and train dataset train_proportion <- 2/3 train_size <-
# round(train_proportion*n) train_indices <- sample(1:n, train_size)
# x_train <- x[train_indices,] x_test <- x[-train_indices,] x <- x_train
# n <- train_proportion*n f1 <-
# 10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5]) f2 <-
# (5*x[,2])/(1+x[,1]^2)+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5]) f3 <-
# 10*sin(pi*x[,1]*x[,2])+3*x[,3]+2*x[,4]+x[,5] #f4 <-
# (5*x[,2])/(1+x[,1]^2)+3*x[,3]+2*x[,4]+x[,5] #treatment set up z <-
# rep(seq(1:3), each = n/3) gs <- length(unique(z)) y1 <- rnorm(n/gs,
# mean = f1, sd = 2) y2 <- rnorm(n/gs, mean = f2, sd = 2) y3 <-
# rnorm(n/gs, mean = f3, sd = 2) #y4 <- rnorm(n/gs, mean = f4, sd = 2) y
# <- c(y1, y2, y3) x <- x_test f1_test <-
# 10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5]) f2_test <-
# (5*x[,2])/(1+x[,1]^2)+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5]) f3_test <-
# 10*sin(pi*x[,1]*x[,2])+3*x[,3]+2*x[,4]+x[,5] ##########Fit the model
# and predict on test set################ #Treatment 1 ind <-
# which(z==1) out1 <- OthermethodLow(x_train[ind,], y[ind], x_test)
# #Treatment 2 ind <- which(z==2) out2 <- OthermethodLow(x_train[ind,],
# y[ind], x_test) #Treatment 3 ind <- which(z==3) out3 <-
# OthermethodLow(x_train[ind,], y[ind], x_test) ###############Calculate
# the errors for different treatment contrasts on test
# set######################################## mat <- matrix(0, 7, 3)
# rownames(mat) <- c('rvmerror', 'rferror', 'svmerror', 'OLSerror',
# 'BART20', 'BART50', 'BART200') mat[, 1] <- colMeans((out1-out2 -
# (f1_test - f2_test))^2) #rvmerror, rferror, svmerror, OLSerror,
# BART20, BART50, BART200 mat[, 2] <- colMeans((out3-out2 - (f3_test -
# f2_test))^2) #rvmerror, rferror, svmerror, OLSerror, BART20, BART50,
# BART200 mat[, 3] <- colMeans((out1-out3 - (f1_test - f3_test))^2)
# #rvmerror, rferror, svmerror, OLSerror, BART20, BART50, BART200
# ##########Fit the model and predict on train set################
# #Treatment 1 ind <- which(z==1) out1 <- OthermethodLow(x_train[ind,],
# y[ind], x_train) #Treatment 2 ind <- which(z==2) out2 <-
# OthermethodLow(x_train[ind,], y[ind], x_train) #Treatment 3 ind <-
# which(z==3) out3 <- OthermethodLow(x_train[ind,], y[ind], x_train) x
# <- x_train f1_test <-
# 10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5]) f2_test <-
# (5*x[,2])/(1+x[,1]^2)+20*(x[,3]-0.5)^2+(10*x[,4]+5*x[,5]) f3_test <-
# 10*sin(pi*x[,1]*x[,2])+3*x[,3]+2*x[,4]+x[,5] ###############Calculate
# the errors for different treatment contrasts on train
# set######################################## mat <- matrix(0, 7, 3)
# rownames(mat) <- c('rvmerror', 'rferror', 'svmerror', 'OLSerror',
# 'BART20', 'BART50', 'BART200') mat[, 1] <- colMeans((out1-out2 -
# (f1_test - f2_test))^2) #rvmerror, rferror, svmerror, OLSerror,
# BART20, BART50, BART200 mat[, 2] <- colMeans((out3-out2 - (f3_test -
# f2_test))^2) #rvmerror, rferror, svmerror, OLSerror, BART20, BART50,
# BART200 mat[, 3] <- colMeans((out1-out3 - (f1_test - f3_test))^2)
# #rvmerror, rferror, svmerror, OLSerror, BART20, BART50, BART200
