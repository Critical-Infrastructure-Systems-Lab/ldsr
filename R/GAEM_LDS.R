#' Converts theta from a vector (as used in GA) to list (as used in Kalman smoothing)
#'
#' @param theta.vec a vector of parameter elements
#' @param d dimention of inputs
#' @keywords internal
vec_to_list <- function(theta.vec, d) {
    d2 <- d*2
    list(A = matrix(theta.vec[1]),
         B = matrix(theta.vec[2:(d+1)],1,d),
         C = matrix(theta.vec[d+2]),
         D = matrix(theta.vec[(d+3):(d2+2)], 1, d),
         Q = matrix(theta.vec[d2+3]),
         R = matrix(theta.vec[d2+4]),
         mu1 = matrix(theta.vec[d2+5]),
         V1 = matrix(theta.vec[d2+6]))
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
    ks.result <- Kalman_smoother(y, u, v, theta.list)

    # Calculate sum of squares of measurement update
    A <- theta.list$A
    B <- theta.list$B
    N <- ncol(u)
    X.tp1 <- ks.result$X[1, 2:N]
    X.t <- ks.result$X[1,1:(N-1)]
    ssq <- sum((X.tp1 - A %*% X.t - B %*% u[,1:(N-1)])^2)
    # Because likelihood is normalized, ssq needs to be normalized as well
    ks.result$lik - lambda * ssq / sum(is.na(y))
}

GAEM_LDS <- function(y, u, v, lambda = 1, num.restart = 10, niter = 100, pop.size = 100, parallel = TRUE) {

    # Upper and lower bounds
    lb <- c(0.1,   -1, -1, -1, -1, -1, 0.001, -1,
            0.001, -1, -1, -1, -1, -1, 0.001, -1,
            0.001, 0.001, 0, 0.001)
    ub <- c(0.999, -0.001, -0.001, -0.001, -0.001, -0.001, 1, -0.001,
            0.999, -0.001, -0.001, -0.001, -0.001, -0.001, 1, -0.001,
            10, 1, 0, 1)
    # Run genetic algorithm
    d <- nrow(u)
    GA <- GA::gaisl(type = 'real-valued',
             fitness = function(theta) penalized_likelihood(y, u, v, theta, lambda),
             lower = lb,
             upper = ub,
             popSize = pop.size,
             numIslands = num.restart,
             migrationRate = 0,
             names = c('A', paste0('B', 1:d), 'C', paste0('D', 1:d), 'Q', 'R', 'mu1', 'V1'),
             run = 100,
             maxiter = niter,
             monitor = FALSE,
             optim = TRUE,
             parallel = parallel)
    fit <- summary(GA)
    best <- which.max(fit$fitnessValues)
    theta.vec <- fit$solution[best,]

    # Reconstruction with measurement update
    theta <- vec_to_list(theta.vec, d)
    ks.result <- Kalman_smoother(y, u, v, theta)
    list(theta = theta,
         fit = ks.result,
         lik = ks.result$lik)
}
