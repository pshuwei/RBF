scheme_mult <- function(y, x, x_test, z, gs, sy, add, skip = 5, sigma_k = NULL,
    po = 0.5, C = 1, Total_itr, burn) {
    #### Finding optimal K#### scaling y
    my <- min(y) + 0.5 * (max(y) - min(y))
    sy <- (max(y) - min(y))

    y <- (y - my)/sy

    gs <- length(unique(z))
    Kls <- rep(14, gs)

    index_list <- list()
    for (g in 1:gs) {
        index_list[[g]] <- which(z == g)

        out <- try(foo <- kernlab::rvm(x[index_list[[g]], ], y[index_list[[g]]],
            kernel = "rbfdot"), silent = T)

        if (class(out) != "try-error") {
            Kls[g] <- foo@nRV
        }
    }

    K <- gs * max(Kls) + add

    print(K)

    p <- ncol(x)
    p_j <- 0.5

    #### Functions#####

    # generating X_mat
    X_mat_fun2 <- function(X, mu) {
        X_mat <- matrix(0, nrow(X), K)
        for (k in 1:K) {
            X_mat[, k] <- exp(-(sqrt(rowSums(t(t(X) - mu[, k])^2))/sigma_k[k])^po)
        }

        return(X_mat)
    }

    # generating X_star_mat
    X_star_fun2 <- function(X_mat, gamma) {
        X_star_mat <- matrix(0, nrow(X_mat), K)
        for (g in 1:gs) {
            for (k in 1:K) {
                X_star_mat[index_list[[g]], k] <- gamma[g, k] * X_mat[index_list[[g]],
                  k]
            }
        }
        return(X_star_mat)
    }

    X_star_est <- function(X_mat, gamma, trt) {
        X_star_e <- matrix(0, nrow(X_mat), K)
        for (k in 1:K) {
            X_star_e[, k] <- gamma[trt, k] * X_mat[, k]
        }
        return(X_star_e)
    }

    # generating V

    V_func <- function(gamma, theta) {
        V_mat <- matrix(0, gs, K)
        for (g in 1:gs) {
            for (k in 1:K) {
                V_mat[g, k] <- gamma[g, k] * theta[k]
            }
        }
        return(V_mat)
    }

    # generating y-XV (for any situation, accounting for group)

    y_xv_func <- function(y, X_mat, V, g) {
        # y_xv <- matrix(0, gs) for (g in 1:2) {
        y_xv <- sum((y[index_list[[g]]] - X_mat[index_list[[g]], ] %*% V[g,
            ])^2)
        # }
        return(y_xv)
    }

    llhood_fun <- function(y, mu) {
        X_mat <- X_mat_fun2(x, mu)
        y_xv <- sum(sapply(1:gs, function(g) y_xv_func(y, X_mat, V, g)))
        each <- -(1/(2 * sigma^2)) * y_xv  #add prior information
        return(each - sum(mu^2/(2 * (sigma_mu)^2)))
    }

    derivmu2 <- function(y, x, V, mu, K, p, gs) {
        X_mat <- X_mat_fun2(x, mu)
        score_mat <- matrix(0, p, K)
        groups <- matrix(0, gs, p)
        for (k in 1:K) {
            for (g in 1:gs) {
                pr <- po * (sqrt(rowSums(t(t(x[index_list[[g]], ]) - mu[,
                  k])^2))/sigma_k[k])^(po - 1)
                groups[g, ] <- (1/sigma^2/2) * V[g, k] * colSums((t(t(x[index_list[[g]],
                  ]) - mu[, k])/sigma_k[k]^2) * c((y[index_list[[g]]] - X_mat[index_list[[g]],
                  ] %*% array(V[g, ])) * (X_mat[index_list[[g]], k]) * (pr *
                  sigma_k[k]/sqrt(rowSums(t(t(x[index_list[[g]], ]) - mu[,
                  k])^2)))))
            }
            score_mat[, k] <- colSums(groups)
        }
        return(score_mat - mu/sigma_mu^2)
    }

    #### Initializing####

    mu_ini <- as.matrix(t(wskm::ewkm(x, K)$centers)) + matrix(rnorm(p * K,
        mean = 0, sd = 0.01), p, K)

    sigma_mu <- 1  #max(abs(x))

    if (is.null(sigma_k)) {
        for (k in 1:K) {
            mat <- as.matrix(dist(rbind(x, mu_ini[, k])))
            sigma_k[k] <- mean(dist(t(mu_ini))) * sqrt(2)  #*mean(mat[n+1,1:n])#mean(dist(t(mu_ini)))
        }
    }

    gamma_ini <- matrix(0, gs, K)

    # what is the best way to initialize gamma here?
    for (g in 1:gs) {
        gamma_ini[g, c(1:Kls[g] + sum(Kls[1:g - 1]))] <- 1
    }

    X_mat <- X_mat_fun2(x, mu_ini)
    X_test_mat <- X_mat_fun2(x_test, mu_ini)

    X_star_mat <- X_star_fun2(X_mat, gamma_ini)

    cx <- crossprod(X_star_mat)
    theta_ini <- solve(cx + 0.01 * mean(diag(cx)) * diag(nrow(cx))) %*% colSums(X_star_mat *
        y)
    theta_0 <- rep(0, K)

    V <- V_func(gamma_ini, theta_ini)

    y_xv <- sum(sapply(1:gs, function(g) y_xv_func(y, X_mat, V, g)))

    sigmahat <- sd(lm(y ~ x + z)$residuals)

    theta <- theta_ini
    sigma <- sigmahat
    gammas <- gamma_ini
    mu <- mu_ini

    llhood <- llhood_fun(y, mu_ini)

    itr <- 0
    epsilon1 <- 0.1
    flag1 <- 0

    p_g <- rep(0.5, gs)

    # Define lists to store postburn samples thetals <- list() gammals
    # <- list() sigmals <- list() f1ls <- list() f2ls <- list() f3ls <-
    # list() #f4ls <- list() muls <- list() treatls <- list()
    tau21ls <- list()
    tau31ls <- list()

    Total_itr <- Total_itr
    burn <- burn  #Sample before itr==burn are in burn-in period.

    # eta <- 1 alpha <- 3 thr <- pgamma(1/sigmahat^2, alpha, eta)
    # while(thr > 0.10){ eta <- eta*0.95 thr <- pgamma(1/sigmahat^2,
    # alpha, eta) } print(eta)

    K1 <- max(rowSums(gammas))

    D_0inv <- diag(rep(C * (K1 + 1), K))
    while (itr < Total_itr) {

        start_time <- Sys.time()

        itr <- itr + 1

        if (itr == 2000) {
            K1 <- max(rowSums(gammas))

            D_0inv <- diag(rep(C * (K1 + 1), K))
        }


        ##### Gibbs Sampling#####
        mu0 <- (y - X_star_mat %*% theta)/sigma^2

        var0 <- (n/sigma^2 + D_0inv[1, 1])^(-1)

        muinter <- rnorm(1, mu0 * var0, sqrt(var0))

        yred <- y - muinter

        ### Sample theta conditional given gamma, sigma, y###
        ### theta_tilde <- solve(solve(D_0) + (1/sigma)*t(X_star_mat)
        ### %*% X_star_mat) %*% solve(D_0) %*% theta_0 +
        ### (1/sigma)*rowSums(X_star_mat *y)

        # D <- solve(solve(D_0) + (1/sigma)*x_star %*% x_star)

        D <- solve(D_0inv + (1/sigma^2) * crossprod(X_star_mat))

        D <- (D + t(D))/2

        theta_tilde <- D %*% ((D_0inv %*% theta_0) + (1/sigma^2) * (t(X_star_mat) %*%
            yred))


        # theta <- rnorm(10, theta_tilde, D)
        theta <- array(rmvnorm(1, mean = theta_tilde, sigma = D))

        # a <- 1/rgamma(1, 1/2, C*(K+1)) alpha <- 0.5 eta <- 1/a gen <-
        # rgamma(1, alpha + K/2+1/2, eta + sum(theta^2)/2+muinter^2/2)

        # a <- 1/rgamma(1, (K1+1)/2+1, sum(theta^2)/2+muinter^2/2) gen
        # <- D_0inv[1, 1] D_check <- dcauchy(sqrt(a),0, sqrt(C*(K1+1)),
        # log=T)-dcauchy(sqrt(D_0inv[1,1]),0, sqrt(C*(K1+1)), log=T)
        # D_check <- D_check + 3*(log(sqrt(a))-log(sqrt(gen))) u <-
        # runif(1) if(log(u) < D_check){ gen <- a } D_0inv <-
        # gen*diag(K)

        ### Sample gamma_j given gamma_-j, theta, sigma, y###
        if (itr < 1000) {
            if (itr%%1 == 0) {
                for (g in 1:gs) {
                  for (j in c(1:K)[-round(c(1:Kls[g] + sum(Kls[1:g - 1])))]) {
                    gamma1 <- gamma0 <- gammas
                    gamma1[g, j] <- 1
                    gamma0[g, j] <- 0
                    # V_star <- gamma1*theta
                    V_star <- V_func(gamma1, theta)
                    # V_2star <- gamma0*theta
                    V_2star <- V_func(gamma0, theta)

                    y_xvstar <- y_xv_func(yred, X_mat, V_star, g)
                    y_xv2star <- y_xv_func(yred, X_mat, V_2star, g)

                    # c_j <- p_j*exp(-1/(2*sigma^2)*y_xvstar) d_j <-
                    # (1-p_j)*exp(-1/(2*sigma^2)*y_xv2star)
                    d_jDivc_j <- (1 - p_g[g]) * exp(-1/(2 * sigma^2) * y_xv2star +
                      1/(2 * sigma^2) * y_xvstar)/p_g[g]
                    p_j_tilde <- 1/(1 + d_jDivc_j)
                    gammas[g, j] <- rbinom(1, 1, p_j_tilde)
                  }
                  p_g[g] <- rbeta(1, sum(gammas[g, ]) - Kls[g] + 1, K - sum(gammas[g,
                    ]) + 1)
                }
            }
        }

        if (itr >= 1000) {
            if (itr%%1 == 0) {
                for (g in 1:gs) {
                  for (j in c(1:K)) {
                    # [-round(c(1:Kls[g]+sum(Kls[1:g-1])))]
                    gamma1 <- gamma0 <- gammas
                    gamma1[g, j] <- 1
                    gamma0[g, j] <- 0
                    # V_star <- gamma1*theta
                    V_star <- V_func(gamma1, theta)
                    # V_2star <- gamma0*theta
                    V_2star <- V_func(gamma0, theta)

                    y_xvstar <- y_xv_func(yred, X_mat, V_star, g)
                    y_xv2star <- y_xv_func(yred, X_mat, V_2star, g)

                    # c_j <- p_j*exp(-1/(2*sigma^2)*y_xvstar) d_j <-
                    # (1-p_j)*exp(-1/(2*sigma^2)*y_xv2star)
                    d_jDivc_j <- (1 - p_g[g]) * exp(-1/(2 * sigma^2) * y_xv2star +
                      1/(2 * sigma^2) * y_xvstar)/p_g[g]
                    p_j_tilde <- 1/(1 + d_jDivc_j)
                    gammas[g, j] <- rbinom(1, 1, p_j_tilde)
                  }
                  p_g[g] <- rbeta(1, sum(gammas[g, ]) + 1, K - sum(gammas[g,
                    ]) + 1)
                }
            }
        }


        V <- V_func(gammas, theta)
        X_star_mat <- X_star_fun2(X_mat, gammas)

        y_xv <- sum(sapply(1:gs, function(g) y_xv_func(yred, X_mat, V, g)))

        a <- 1/rgamma(1, 1 + n/2, y_xv/2)
        gen <- sigma

        D_check <- dcauchy(sqrt(a), 0, sigmahat, log = T) - dcauchy(sigma,
            0, sigmahat, log = T)
        D_check <- D_check + 3 * (log(sqrt(a)) - log(sqrt(gen)))

        u <- runif(1)

        if (log(u) < D_check) {
            gen <- sqrt(a)
        }
        sigma <- gen

        #### Gradient MH#####
        if (itr > 0) {
            if (itr%%1 == 0) {
                dermu <- derivmu2(yred, x, V, mu, K, p, gs)  #, 0)
                # for random walk dermu <- 0

                temp <- mu + dermu * epsilon1/2 + matrix(rnorm(1, sd = sqrt(epsilon1)),
                  p, K)

                muc <- temp

                llhoodc <- llhood_fun(yred, muc)

                llhood <- llhood_fun(yred, mu)

                q1 <- -sum((mu - derivmu2(yred, x, V, temp, K, p, gs) * epsilon1/2 -
                  temp)^2/epsilon1)/2  #, 0)
                q2 <- -sum((temp - derivmu2(yred, x, V, mu, K, p, gs) * epsilon1/2 -
                  mu)^2/epsilon1)/2  #, 0)
                D_check <- llhoodc - llhood + q1 - q2
                u <- runif(1)

                if (log(u) < D_check) {
                  # u < exp(D) [=ratio of likelihoods]
                  # P(log(u)<D)=P(u<exp(D))=exp(D)
                  mu <- muc
                  llhood <- llhoodc
                  X_mat <- X_mat_fun2(x, mu)
                  X_test_mat <- X_mat_fun2(x_test, mu)
                  X_star_mat <- X_star_fun2(X_mat, gammas)
                  y_xv <- sum(sapply(1:gs, function(g) y_xv_func(yred, X_mat,
                    V, g)))
                  flag1 <- flag1 + 1  #To count number of accepts

                }
            }
        }

        #### Update the treatment list####

        f1est <- sy * (X_star_est(X_test_mat, gammas, 1) %*% theta + muinter) +
            my
        f2est <- sy * (X_star_est(X_test_mat, gammas, 2) %*% theta + muinter) +
            my
        f3est <- sy * (X_star_est(X_test_mat, gammas, 3) %*% theta + muinter) +
            my
        # f4est <- sy*X_star_est(X_test_mat, gammas, 4) %*% theta

        tau21est <- f2est - f1est
        tau31est <- f3est - f1est

        if (itr%%200 == 0) {
            if (itr > 0) {
                # Auto-tuning step
                ar <- flag1 * 1/(itr - 0)  #Acceptance rate
                if (ar < 0.45) {
                  # Acceptance is too low
                  epsilon1 <- epsilon1/2
                }
                if (ar > 0.7) {
                  # Acceptance is too high
                  epsilon1 <- epsilon1 * 2
                }
            }
            # print(range(mu))
            if (itr%%400 == 0) {
                # print(D_0inv[1,1])
                for (k in 1:K) {
                  mat <- as.matrix(dist(rbind(x, mu[, k])))
                }
            }

            K1 <- K
        }

        #### Store the posterior samples####
        if (itr > burn) {
            if ((itr - burn)%%skip == 0) {
                # thetals[[(itr-burn)/skip]] <- theta
                # gammals[[(itr-burn)/skip]] <- gammas
                # sigmals[[(itr-burn)/skip]] <- sigma
                # muls[[(itr-burn)/skip]] <- mu f1ls[[(itr-burn)/skip]]
                # <- f1est f2ls[[(itr-burn)/skip]] <- f2est
                # f3ls[[(itr-burn)/skip]] <- f3est
                # f4ls[[(itr-burn)/skip]] <- f4est
                # treatls[[(itr-burn)/skip]] <- trt_effect

                tau21ls[[(itr - burn)/skip]] <- tau21est
                tau31ls[[(itr - burn)/skip]] <- tau31est
            }
        }
        end_time <- Sys.time()
    }

    # thetas <- Reduce('+', thetals)/length(thetals) gammas <-
    # Reduce('+', gammals)/length(gammals) sigmas <- Reduce('+',
    # sigmals)/length(sigmals) mus <- Reduce('+', muls)/length(muls)
    # f1hat <- Reduce('+', f1ls)/length(f1ls) f2hat <- Reduce('+',
    # f2ls)/length(f2ls) f3hat <- Reduce('+', f3ls)/length(f3ls) #f4hat
    # <- Reduce('+', f4ls)/length(f4ls) f_list <- list(theta_est =
    # thetas, gamma_est = gammas, sigma_est = sigmas, mu_est = mus,
    # f1est = f1hat, f2est = f2hat, f3est = f3hat #f4est = f4hat )
    f_list2 <- list(tau21ls, tau31ls)

    return(f_list2)
}
