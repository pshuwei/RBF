scheme_mult <- function(y, x, x_test, z, gs, sy, add, skip=5, sigma_k=NULL, C = 1, Total_itr, burn) {
  ####Number of Covariates#### 
  p <- ncol(x)
  
  ####Apply min-max transformation on the response####
  my <- min(y)+0.5*(max(y)-min(y))
  sy  <- (max(y)-min(y))
  
  y <- (y-my)/sy
  
  ####Finding optimal K####
  gs <- length(unique(z))
  Kls <- rep(14, gs)
  
  index_list <- list()
  for (g in 1:gs) {
    #index list to group the response and predictors associated with the g-th treatment
    index_list[[g]] <- which(z == g)
    
    #utilize Gaussian RBF kernel-based relevance vector machine (RVM)
    #used to obtain the group-specific number of relevance vectors
    out <- try(foo <- kernlab::rvm(x[index_list[[g]],], y[index_list[[g]]], 
                                   kernel = "rbfdot"), silent = T)
    
    if(class(out)!="try-error"){Kls[g] <- foo@nRV}
  }
  
  #set $K = G \times \max_g v_g$
  K <- gs*max(Kls)+add
  
  ####Functions#####
  
  #Gaussian radial function for the RBF network for the model
  X_mat_fun2 <- function(X, mu) {
    X_mat <- matrix(0, nrow(X), K)
    for (k in 1:K) {
      X_mat[,k] <- exp(-(sqrt(rowSums(t(t(X) - mu[,k])^2))/sigma_k[k])^2)
    }
    
    return(X_mat)
  }
  
  #Gaussian radial function multiplied by the corresponding treatment group for the i-th patient, used for the posterior sampling of theta
  X_star_fun2 <- function(X_mat, gamma) {
    X_star_mat <- matrix(0, nrow(X_mat), K)
    for (g in 1:gs) {
      for (k in 1:K) {
        X_star_mat[index_list[[g]],k] <- gamma[g, k] * X_mat[index_list[[g]], k]
      }
    }
    return(X_star_mat)
  }
  
  #Similar function as above, but for a hypothetical different treatment group from the assigned treatment group, used to obtain CATE estimators
  X_star_est <- function(X_mat, gamma, trt) {
    X_star_e <- matrix(0, nrow(X_mat), K)
    for (k in 1:K) {
      X_star_e[,k] <- gamma[trt,k] * X_mat[, k]
    }
    return(X_star_e)
  }
  
  #Matrix multiplication of theta and gamma, used to simplify coding
  V_func <- function(gamma, theta) {
    V_mat <- matrix(0, gs, K)
    for (g in 1:gs) {
      for (k in 1:K) {
        V_mat[g, k] <-gamma[g, k] * theta[k]
      }
    }
    return(V_mat)
  }
  
  #Obtain errors, notice that from the model, \epsilon_i = y_i - f_g(x_i). We can use this for gradient sampling of \mu
  y_xv_func <- function(y, X_mat, V, g) {
    # y_xv <- matrix(0, gs)
    #for (g in 1:2) {
    y_xv <- sum((y[index_list[[g]]] - X_mat[index_list[[g]],] %*% V[g,])^2)
    #}
    return(y_xv)
  }
  
  #Conditional posterior log-likelihood function of model, used for gradient sampling of \mu
  llhood_fun <- function(y, mu) {
    X_mat <- X_mat_fun2(x, mu)
    y_xv <-sum(sapply(1:gs, function(g) y_xv_func(y, X_mat, V, g)))
    each <- -(1/(2*sigma^2))*y_xv #add prior information
    return(each - sum(mu^2/(2*(sigma_mu)^2)))
  }
  
  #Derivative of the above log-likelihood function
  derivmu2 <- function(y, x, V, mu, K, p, gs) {
    X_mat <- X_mat_fun2(x, mu)
    score_mat <- matrix(0, p, K)
    groups <- matrix(0, gs, p)
    for (k in 1:K) {
      for (g in 1:gs) {
        pr <- 2*(sqrt(rowSums(t(t(x[index_list[[g]],]) - mu[,k])^2))/sigma_k[k])
        groups[g,] <- (1/sigma^2/2)*V[g,k]*colSums((t(t(x[index_list[[g]],]) - mu[,k])/sigma_k[k]^2) * c((y[index_list[[g]]] - X_mat[index_list[[g]],] %*% array(V[g,])) * (X_mat[index_list[[g]],k])*(pr*sigma_k[k]/sqrt(rowSums(t(t(x[index_list[[g]],]) - mu[,k])^2)))))
      }
      score_mat[,k] <- colSums(groups)
    }
    return(score_mat- mu/sigma_mu^2)
  }
  
  ####Initializing####
  
  #initializing \mu with the entropy weighted k-means clustering algorithm
  mu_ini <- as.matrix(t(wskm::ewkm(x, K)$centers))+matrix(rnorm(p*K, mean=0, sd=0.01), p,K)
  
  sigma_mu <- 1#max(abs(x))
  
  #pre-specified bandwidth parameters \sigma_k's, discussed in "Other computational settings"
  if(is.null(sigma_k)){
    for(k in 1:K){ 
      mat <- as.matrix(dist(rbind(x, mu_ini[,k])))
      sigma_k[k] <- mean(dist(t(mu_ini)))*sqrt(2)#*mean(mat[n+1,1:n])#mean(dist(t(mu_ini)))
    }
  }
  
  p_j <- 0.5
  
  #initializing \Gamma based on v_g, K
  #creating G disjoints sets, and fixing small set of indices to be 1
  gamma_ini <- matrix(0, gs, K)
  
  for (g in 1:gs) {
    gamma_ini[g, c(1:Kls[g]+sum(Kls[1:g-1]))] <- 1
  }
  
  #initializing other values using the functions earlier created
  X_mat <- X_mat_fun2(x, mu_ini)
  X_test_mat <- X_mat_fun2(x_test, mu_ini)
  X_star_mat <- X_star_fun2(X_mat, gamma_ini)
  cx <- crossprod(X_star_mat)
  llhood <- llhood_fun(y, mu_ini) 

  #initializing theta
  theta_ini <- solve(cx+1e-2*mean(diag(cx))*diag(nrow(cx))) %*% colSums(X_star_mat * y)
  #\btheta_0 used for the full conditional of \btheta
  theta_0 <- rep(0, K)
  
  #initializing other functions that required initializion of theta
  V <- V_func(gamma_ini, theta_ini)
  y_xv <- sum(sapply(1:gs, function(g) y_xv_func(y, X_mat, V, g)))
  
  #error variance based on an estimate of \sigma from the data 
  #in our case is the error standard deviation from the linear regression model
  sigmahat <- sd(lm(y~x+z)$residuals)
  
  #set the first iteration values to the be the initial values
  theta <- theta_ini
  sigma <- sigmahat
  gammas <- gamma_ini
  mu <- mu_ini
  
  itr <- 0 #iterations
  epsilon1 <- 1e-1 #error
  flag1 <- 0 #flag for MCMC sampling
  
  p_g <- rep(0.5,gs)
  
  #Define lists to store postburn samples
  thetals <- list()
  gammals <- list()
  sigmals <- list()
  f1ls <- list()
  f2ls <- list()
  f3ls <- list()
  muls <- list()
  treatls <- list()
  
  #Define lists to store CATE estimators
  tau21ls <- list()
  tau31ls <- list()
  
  Total_itr <- Total_itr
  burn <- burn #Sample before itr==burn are in burn-in period.
  
  #Initialize K1 as the minimum number of bases from all groups
  K1 <- max(rowSums(gammas))
  
  #initialization of \sigma_d^{-2}
  D_0inv <- diag(rep(C*(K1+1), K))
  
  while(itr < Total_itr){
    
    start_time <- Sys.time()
    
    itr <- itr + 1
    
    if(itr==2000){
      K1 <- max(rowSums(gammas))
      
      D_0inv <- diag(rep(C*(K1+1), K))
    }
    
    
    #### Sampler ####
    #Sample \alpha (the intercept)
    alpha0 <- (y-X_star_mat %*% theta)/sigma^2
    
    sigma0 <- (n/sigma^2+D_0inv[1,1])^(-1)
    
    alpha <- rnorm(1, alpha0*sigma0, sqrt(sigma0))
    
    yred <- y-alpha #the residuals R
    
    #Sample \theta conditional given gamma, sigma, y#
    
    ##find which columns contain only 0, sample from prior
    zero_column_indices <- which(colSums(X_star_mat) == 0)
    
    is_zero_column <- seq_len(K) %in% zero_column_indices
    
    theta[is_zero_column] <- rnorm(sum(is_zero_column), mean = 0, sd = D_0inv[1, 1]^-1)
    
    ##sample from full conditional for non-zero columns
    D <- solve(D_0inv[!is_zero_column, !is_zero_column] + (1/sigma^2)*crossprod(X_star_mat[, !is_zero_column]))
    
    D  <- (D + t(D))/2
    
    theta_tilde <- D %*% ((D_0inv[!is_zero_column, !is_zero_column] %*% theta_0[ !is_zero_column]) + (1/sigma^2)*(t(X_star_mat[, !is_zero_column])%*%yred))
    
    theta[!is_zero_column] <- array(rmvnorm(1,mean=theta_tilde,sigma=D))
    
    #Sample gamma_j given gamma_-j, theta, sigma, y#
    ##for the first 1000 iterations, only sample where gamma_kg is fixed at 1
    if(itr<1000){
      if (itr %% 1==0) {
        for (g in 1:gs) {
          for (j in c(1:K)[-round(c(1:Kls[g]+sum(Kls[1:g-1])))]) { 
            gamma1 <- gamma0 <- gammas
            gamma1[g, j] <- 1
            gamma0[g, j] <- 0
            
            V_star <- V_func(gamma1, theta)
            V_2star <- V_func(gamma0, theta)
            
            y_xvstar <- y_xv_func(yred, X_mat, V_star, g)
            y_xv2star <- y_xv_func(yred, X_mat, V_2star, g)
            
            d_jDivc_j <- (1-p_g[g])*exp(-1/(2*sigma^2)*y_xv2star + 1/(2*sigma^2)*y_xvstar)/p_g[g]
            p_j_tilde <- 1/(1+d_jDivc_j)
            gammas[g,j] <- rbinom(1, 1, p_j_tilde)
          }
          p_g[g]<- rbeta(1,sum(gammas[g,])-Kls[g]+1,K-sum(gammas[g,])+1)
        }
      }
    }
    ##afterwards, sample all gamma_kg
    if(itr>=1000){
      if (itr %% 1==0) {
        for (g in 1:gs) {
          for (j in c(1:K)) {
            gamma1 <- gamma0 <- gammas
            gamma1[g, j] <- 1
            gamma0[g, j] <- 0

            V_star <- V_func(gamma1, theta)
            V_2star <- V_func(gamma0, theta)
            
            y_xvstar <- y_xv_func(yred, X_mat, V_star, g)
            y_xv2star <- y_xv_func(yred, X_mat, V_2star, g)
            
            d_jDivc_j <- (1-p_g[g])*exp(-1/(2*sigma^2)*y_xv2star + 1/(2*sigma^2)*y_xvstar)/p_g[g]
            p_j_tilde <- 1/(1+d_jDivc_j)
            gammas[g,j] <- rbinom(1, 1, p_j_tilde)
          }
          p_g[g]<- rbeta(1,sum(gammas[g,])+1,K-sum(gammas[g,])+1)
        }
      }
    }
    
    #updating certain values with new samples
    V <- V_func(gammas, theta)
    X_star_mat <- X_star_fun2(X_mat, gammas)
    
    y_xv <- sum(sapply(1:gs, function(g) y_xv_func(yred, X_mat, V, g)))
    
    #Sampling sigma  
    a <- 1/rgamma(1, 1 + n/2, y_xv/2)
    gen <- sigma
    
    ##Log transformation of the acceptance probability ratio
    D_check <- dcauchy(sqrt(a),0, sigmahat, log=T)-dcauchy(sigma,0, sigmahat, log=T)
    ##Log transformation of the Jacobian adjustment
    D_check <- D_check + 3*(log(sqrt(a))-log(sqrt(gen)))
    
    u <- runif(1)
    
    ##Acceptance step
    if(log(u) < D_check){ 
      gen <- sqrt(a)
    }
    sigma <- gen
    
    #Sampling mu_kp with gradient MH
    if(itr > 000){
      if (itr %% 1==0) {
        #generate derivative
        dermu <- derivmu2(yred, x, V, mu, K, p, gs)#, 0)
        #The update is like ’gradient-based update from the current value + some noise’
        temp <- mu + dermu * epsilon1/2 + matrix(rnorm(1, sd = sqrt(epsilon1)), p, K)
        #Set the update as the candidate value
        muc <- temp
        #log-likelihood of the candidate excluding the part not involving the candidate
        llhoodc <- llhood_fun(yred, muc
        #log-likelihood of the current value excluding the part not involving the current value
        llhood <- llhood_fun(yred, mu)
        
        #log{q(current value|candidate)}
        q1 <- -sum((mu-derivmu2(yred, x, V, temp, K, p, gs)* epsilon1/2- temp)^2/epsilon1) /2#, 0)
        #log{q(candidate|current value)}
        q2 <- -sum((temp-derivmu2(yred, x, V, mu, K, p, gs) * epsilon1/2- mu)^2/epsilon1) /2#, 0)
        
        ##Log transformation of the acceptance probability ratio
        D_check <- llhoodc - llhood + q1 - q2
        u <- runif(1)
        
        if(log(u) < D_check){
          mu <- muc
          llhood <- llhoodc #change the candidate value to the current value
          # Then update the function values based on the new samples
          X_mat <- X_mat_fun2(x, mu)
          X_test_mat <- X_mat_fun2(x_test, mu)
          X_star_mat <- X_star_fun2(X_mat, gammas)
          y_xv <- sum(sapply(1:gs, function(g) y_xv_func(yred, X_mat, V, g)))
          flag1 <- flag1 + 1 #To count number of accepts
          
        }
      }
    }
    
    #Keep the acceptance rate at an optimal rate (~0.6)
    if(itr %% 200==0){ 
      if(itr > 000){
        #Auto-tuning step
        ar <- flag1*1 / (itr-000) #Acceptance rate
        if(ar < 0.45){ #Acceptance is too low
          epsilon1 <- epsilon1 / 2
        }
        if(ar > 0.7){ #Acceptance is too high
          epsilon1 <- epsilon1 * 2
        }
      }
      #print(range(mu))
      if(itr %% 400==0){
        #print(D_0inv[1,1])
        for(k in 1:K){ 
          mat <- as.matrix(dist(rbind(x, mu[,k])))
        }}
      
      K1 <- K
    }
    
    ####Update the treatment list####
    f1est <- sy*(X_star_est(X_test_mat, gammas, 1) %*% theta+alpha) + my
    f2est <- sy*(X_star_est(X_test_mat, gammas, 2) %*% theta+alpha) + my
    f3est <- sy*(X_star_est(X_test_mat, gammas, 3) %*% theta+alpha) + my
    
    ####Store the posterior samples####
    if(itr > burn){
      if((itr-burn)%%skip==0){
        thetals[[(itr-burn)/skip]] <- theta
        gammals[[(itr-burn)/skip]] <- gammas
        sigmals[[(itr-burn)/skip]] <- sigma
        muls[[(itr-burn)/skip]] <- mu
        f1ls[[(itr-burn)/skip]] <- f1est
        f2ls[[(itr-burn)/skip]] <- f2est
        f3ls[[(itr-burn)/skip]] <- f3est
        
        #getting the CATE estimations samples
        tau21ls[[(itr - burn)/skip]] <- tau21est
        tau31ls[[(itr - burn)/skip]] <- tau31est
      }
    }
    end_time <- Sys.time()    
  }
  
  #Get the average sample based on all iterations
  thetas <- Reduce('+', thetals)/length(thetals)
  gammas <- Reduce('+', gammals)/length(gammals)
  sigmas <- Reduce('+', sigmals)/length(sigmals)
  mus <- Reduce('+', muls)/length(muls)
  f1hat <- Reduce('+', f1ls)/length(f1ls)
  f2hat <- Reduce('+', f2ls)/length(f2ls)
  f3hat <- Reduce('+', f3ls)/length(f3ls)
  
  #Store final samples in list
  f_list <- list(theta_est = thetas,
                 gamma_est = gammas,
                 sigma_est = sigmas,
                 mu_est = mus,
                 f1est = f1hat,
                 f2est = f2hat,
                 f3est = f3hat,
                 tau21 = tau21ls, 
                 tau31 = tau31ls
  )
  
  return(f_list)
}
