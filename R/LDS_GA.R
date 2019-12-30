#' Converts theta from a vector (as used in GA) to list (as used in Kalman smoothing)
#'
#' @param theta.vec a vector of parameter elements
#' @param d dimention of inputs
vec_to_list <- function(theta.vec, d) {
    d2 <- d*2
    list(A = matrix(theta.vec[1]),
         B = matrix(theta.vec[2:(d+1)],1,d),
         C = matrix(theta.vec[d+2]),
         D = matrix(theta.vec[(d+3):(d2+2)], 1, d),
         Q = matrix(theta.vec[d2+3]),
         R = matrix(theta.vec[d2+4]),
         mu1 = matrix(0),
         V1 = matrix(theta.vec[d2+5]))
}

#' Penalized likelihood objective function
#'
#' [Kalman_smoother] returns X and likelihood.
#' The penalized likelihood is the likelihood minus the sum-of-squares of the
#' measurement update. This is used as the fitness function in genetic algorihm.
#' @inheritParams Kalman_smoother
#' @param theta.vec a vector of parameter elements (i.e, the vectorized version of theta
#' in \code{Kalman_smoother})
#' @param lambda weight of the penalty
#' @return The penalized likelihood (a real number)
#' @export
penalized_likelihood <- function(y, u, v, theta.vec, lambda) {

    theta.list <- vec_to_list(theta.vec, nrow(u))
    ks.result <- Kalman_smoother(y, u, v, theta.list, stdlik = FALSE)

    # Calculate sum of squares of measurement update
    A <- theta.list$A
    B <- theta.list$B
    N <- ncol(u)
    X.tp1 <- ks.result$X[1, 2:N]
    X.t <- ks.result$X[1,1:(N-1)]
    ssq <- sum((X.tp1 - A %*% X.t - B %*% u[,1:(N-1)])^2)
    # X <- rep(0, N)
    # for (i in 2:N) X[i] <- A * X[i-1] + B %*% u[,i-1]
    # ssq <- sum((ks.result$X - X)^2)
    ks.result$lik - lambda * ssq
}

#' Learn a linear dynamical system using Genetic Algorithm.
#'
#' **Warning** This is an experimental feature. Use with care.
#' @inheritParams LDS_reconstruction
#' @param y Transformed and standardized streamflow
#'

LDS_GA <- function(y, u, v, lambda = 1, ub, lb, num.islands = 4, pop.per.island = 100, niter = 1000, parallel = TRUE) {

    # Run genetic algorithm
    d <- nrow(u)
    GA <- GA::gaisl(type = 'real-valued',
             fitness = function(theta) penalized_likelihood(y, u, v, theta, lambda),
             lower = lb,
             upper = ub,
             popSize = pop.per.island * num.islands,
             numIslands = num.islands,
             # migrationRate = 0.01,
             names = c('A', paste0('B', 1:d), 'C', paste0('D', 1:d), 'Q', 'R',   'V1'),
             run = 100,
             maxiter = niter,
             monitor = FALSE,
             optim = FALSE, # optim = TRUE gives ugly results
             parallel = parallel)
    # fit <- summary(GA)
    # best <- which.max(fit$fitnessValues)
    theta.vec <- GA@solution

    # Reconstruction with measurement update
    theta <- vec_to_list(theta.vec, d)
    ks.result <- Kalman_smoother(y, u, v, theta) # Return standardized lik so that it's comparable with EM
    list(theta = theta,
         fit = ks.result,
         lik = ks.result$lik,
         pl = GA@fitnessValue)
}

#' Learn LDS with L-BFGS-B
#'
#' **Warning** This is an experimental feature. Use with care.
#' @inheritParams LDS_GA
#' @inheritParams LDS_reconstruction
#' @export
LDS_BFGS <- function(y, u, v, lambda = 1, ub, lb, num.restarts = 100, parallel = TRUE) {

    d <- nrow(u)
    # Initial guess for L-BFGS-B
    par.list <- replicate(num.restarts,
                          runif(d + d + 5, min = lb, max = ub),
                          simplify = FALSE)
    if (parallel) {
        nbCores <- detectCores() - 1
        cl <- makeCluster(nbCores)
        registerDoParallel(cl)
        optim.result <- foreach(par = par.list, .packages = 'stats') %dopar%
            stats::optim(par,
                  function(theta) -penalized_likelihood(y, u, v, theta, lambda),
                  method = 'L-BFGS-B',
                  lower = lb,
                  upper = ub)
        stopCluster(cl)
    } else {
        optim.result <- lapply(par.list, function(par)
            stats::optim(par,
                  function(theta) -penalized_likelihood(y, u, v, theta, lambda),
                  method = 'L-BFGS-B',
                  lower = lb,
                  upper = ub))
    }

    optim.vals <- sapply(optim.result, '[[', 'value')
    best <- which.max(optim.vals)
    theta.vec <- optim.result[[best]]$par

    # Reconstruction with measurement update
    theta <- vec_to_list(theta.vec, d)
    ks.result <- Kalman_smoother(y, u, v, theta)
    list(theta = theta,
         fit = ks.result,
         lik = ks.result$lik,
         pl = optim.vals[best])

}
