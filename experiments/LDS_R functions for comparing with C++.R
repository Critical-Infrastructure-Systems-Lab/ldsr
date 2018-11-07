learnLDS_R <- function(y, u, v, mu, init, niter = 1000, tol = 1e-5) {
    
    # Learn a single LDS model using EM starting from the initial conditions specified by init
    
    d <- nrow(u)
    theta <- list(A = init[1], 
                  B = matrix(init[2], 1, d),
                  C = init[3],
                  D = matrix(init[4], 1, d),
                  Q = 1,
                  R = 1,
                  mu1 = 0,
                  V1 = 1)
    
    lik <- rep(0, niter)
    fit <- Kalman_smoother(y, u, v, theta)
    lik[1] <- fit$lik
    for (i in 2:niter) {
        theta <- Mstep(y, u, v, fit)
        fit <- Kalman_smoother(y, u, v, theta)
        lik[i] <- fit$lik
        if (i > 2 && abs(lik[i] - lik[i-1]) < tol && abs(lik[i-1] - lik[i-2]) < tol) break
    }    
    
    # Final propagation
    # N <- ncol(u)
    # X2 <- rep(0, N)
    # X2[1] <- theta$mu1
    # for (t in 2:N) X2[t] <- theta$A * X2[t-1] + theta$B %*% u[, t-1]
    # Y2 <- theta$C * X2 + theta$D %*% v
    # fit$X2 <- X2
    # fit$Y2 <- Y2
    
    return(list(theta = theta, fit = fit, lik = lik[i]))
}

learnLDS_CRR <- function(y, u, v, mu, init, niter = 1000, tol = 1e-5) {
    
    # Learn a single LDS model using EM starting from the initial conditions specified by init
    
    d <- nrow(u)
    theta <- list(A = matrix(init[1]), 
                  B = matrix(init[2], 1, d),
                  C = matrix(init[3]),
                  D = matrix(init[4], 1, d),
                  Q = matrix(1),
                  R = matrix(1),
                  mu1 = matrix(0),
                  V1 = matrix(1))
    
    lik <- rep(0, niter)
    fit <- Kalman_smoother_C(y, u, v, theta)
    lik[1] <- fit$lik
    for (i in 2:niter) {
        theta <- lapply(Mstep(y, u, v, fit), as.matrix)
        fit <- Kalman_smoother_C(y, u, v, theta)
        lik[i] <- fit$lik
        if (i > 2 && abs(lik[i] - lik[i-1]) < tol && abs(lik[i-1] - lik[i-2]) < tol) break
    }    
    
    # Final propagation
    # N <- ncol(u)
    # X2 <- rep(0, N)
    # X2[1] <- theta$mu1
    # for (t in 2:N) X2[t] <- theta$A %*% X2[t-1] + theta$B %*% u[, t-1]
    # Y2 <- theta$C %*% X2 + theta$D %*% v
    # fit$X2 <- X2
    # fit$Y2 <- Y2
    
    return(list(theta = theta, fit = fit, lik = lik[i]))
}

learnLDS_CCR <- function(y, u, v, mu, init, niter = 1000, tol = 1e-5) {
    
    # Learn a single LDS model using EM starting from the initial conditions specified by init
    
    d <- nrow(u)
    theta <- list(A = matrix(init[1]), 
                  B = matrix(init[2], 1, d),
                  C = matrix(init[3]),
                  D = matrix(init[4], 1, d),
                  Q = matrix(1),
                  R = matrix(1),
                  mu1 = matrix(0),
                  V1 = matrix(1))
    
    lik <- rep(0, niter)
    fit <- Kalman_smoother_C(y, u, v, theta)
    lik[1] <- fit$lik
    for (i in 2:niter) {
        theta <- Mstep_C(y, u, v, fit)
        fit <- Kalman_smoother_C(y, u, v, theta)
        lik[i] <- fit$lik
        if (i > 2 && abs(lik[i] - lik[i-1]) < tol && abs(lik[i-1] - lik[i-2]) < tol) break
    }    
    
    # Final propagation
    # N <- ncol(u)
    # X2 <- rep(0, N)
    # X2[1] <- theta$mu1
    # for (t in 2:N) X2[t] <- theta$A %*% X2[t-1] + theta$B %*% u[, t-1]
    # Y2 <- theta$C %*% X2 + theta$D %*% v
    # fit$X2 <- X2
    # fit$Y2 <- Y2
    
    return(list(theta = theta, fit = fit, lik = lik[i]))
}

learnLDS_restart_R <- function(y, u, v, mu, init = NULL, num.restart, niter = 1000, tol = 1e-5, return.init = FALSE, parallel = FALSE) {
    
    # Learn multiple LDS models with randomized initial conditions (specified by num.restart) 
    # and select the one with the highest likelihood
    # Use a pecified grid of initial conditions instead of randomized ones if provided
    
    if (is.null(init)) init <- replicate(num.restart, 
                                         runif(4, min = c(0, -1, 0, -1), max = c(0.9, 1, 1, 1)), 
                                         simplify = F)
    
    if (parallel) {
        models <- foreach(init.val = init) %dopar% learnLDS_R(y, u, v, mu, init.val, niter, tol)
    } else {
        models <- lapply(init, function(init.val) learnLDS_R(y, u, v, mu, init.val, niter, tol))
    }
    
    liks <- sapply(models, '[[', 'lik')
    
    max.ind <- which.max(liks)
    if (return.init) {
        ans <- list(model = models[[max.ind]],
                    init = init[[max.ind]]) 
    } else {
        ans <- models[[which.max(liks)]] 
    }
    return(ans)
}