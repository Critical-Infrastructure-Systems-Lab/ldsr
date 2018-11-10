#' Learn LDS model with multiple initial conditions
#'
#' The initial conditions can either be randomized (specifiled by num.restart) or provided beforehand.
#' @param Qa Observations: a data.table of annual streamflow with at least two columns: year and Qa.
#' @inheritParams LDS_EM_restart
#' @param method Either EM or GA (note: GA is experimental)
#' @param num.restarts if init is not given then num.restarts must be provided. In this case the function will randomize the initial value by sampling uniformly within the range for each parameters (A in \[0, 1\], B in \[-1, 1\], C in \[0, 1\] and D in \[-1, 1\]).
#' @return A list of the following elements
#' * rec: reconstruction results, a data.table with the following columns
#'     - year: calculated from Qa and the length of u
#'     - X: the estimated hidden state
#'     - Xl, Xu: lower and upper range for the 95\% confidence interval of X
#'     - Q: the reconstructed streamflow
#'     - Ql, Qu: lower and upper range for the 95\% confidence interval of Q
#' @export
LDS_reconstruction <- function(Qa, u, v, method = 'EM',
                               init = NULL, num.restarts = 100, return.init = TRUE,
                               lambda = 1, num.islands = 10, pop.size = 250,
                               niter = 1000, tol = 1e-5,
                               parallel = FALSE) {

    # Attach NA and make the y matrix
    mu <- mean(log(Qa$Qa), na.rm = T)
    n.paleo <- ncol(u) - nrow(Qa) # Number of years in the paleo period
    y <- t(c(rep(NA, n.paleo), log(Qa$Qa) - mu))

    if (method == 'EM') {
        # Randomize initial conditions if not given
        if (is.null(init)) {
            d <- nrow(u)
            init <- replicate(num.restart,
                              runif(1 + d + 1 + d,
                                    min = c(0, rep(-1, d), 0, rep(-1, d)),
                                    max = c(1, rep( 1, d), 1, rep( 1, d))),
                              simplify = F)
        } else {
            if (!is.list(init)) # learnLDS expects init is a list
                init <- list(init)
        }

        # Learn multiple models and select the best one. Run in parallel mode if the init list is long
        if (parallel && length(init) > 10) {
            results <- LDS_EM_restart(y, u, v, init, niter, tol, return.init, parallel = TRUE)
        } else {
            if (parallel)
                warning('Initial condition list is short, LDS_EM_restart() is run in sequential mode.')
            results <- LDS_EM_restart(y, u, v, init, niter, tol, return.init, parallel = FALSE)
        }
    } else {
        results <- LDS_GA(y, u, v, lambda, niter = niter, pop.size = 100, parallel = parallel)
    }
    # Construct 95% confidence intervals and return
    with(results, {
        X <- as.vector(fit$X)
        V <- as.vector(fit$V)
        CI.X <- 1.96*sqrt(V)
        Xl <- X - CI.X # Lower range for X
        Xu <- X + CI.X # Upper range for X
        Y <- as.vector(fit$Y + mu)
        CI.Y <- 1.96 * (as.vector(theta$C) * as.vector(fit$V) * as.vector(theta$C) +
                            as.vector(theta$R))
        Yl <- Y - CI.Y # Lower range for Y
        Yu <- Y + CI.Y # Upper range for Y
        # Transform Y to Q and put in a data.table
        rec <- data.table(year = (Qa[1, year] - n.paleo):(Qa[.N, year]),
                          X, Xl, Xu,
                          Q = exp(Y),
                          Ql = exp(Yl),
                          Qu = exp(Yu))

        ans <- list(rec = rec, theta = theta, lik = lik)
        if (return.init) ans$init <- results$init
        ans
    })
}

#' Cross validate LDS model
#'
#' @inheritParams LDS_reconstruction
#' @details Allows different experimental setups:
#'   * init: if given, all cross validation runs start with the same init, otherwise each cross validation run is learned using randomized restarts
#'   * Z: a list of n.reps elements, each is a vector of length k. If given, cross validation will be run on these points (useful when comparing LDS with another reconstruction method); otherwise, randomized cross validation points will be created
#' @export
cvLDS <- function(Qa, u, v, method = 'EM',
                  init = NULL, num.restart = 20,
                  lambda = 1,
                  k, n.reps = 100, niter = 1000, tol = 1e-5, Z = NULL,
                  reupdate = FALSE,
                  parallel = F) {

    mu <- mean(log(Qa$Qa), na.rm = T)
    n.paleo <- ncol(u) - nrow(Qa) # Number of years in the paleo period
    y <- t(c(rep(NA, n.paleo), log(Qa$Qa) - mu))

    obs.ind <- which(!is.na(y))
    n.obs <- length(obs.ind)

    if (missing(k)) k <- ceiling(n.obs / 10)
    if (is.null(Z)) {
        Z <- replicate(n.reps, sample(obs.ind, k), simplify = F)
        Z2 <- lapply(Z, '-', n.paleo)
    } else {
        Z2 <- Z
        Z <- lapply(Z2, '+', n.paleo)
    }

    one_CV <- function(omit, y, u, v, method, init, num.restart, lambda, niter, tol) {
        y2 <- y
        y2[omit] <- NA
        if (method == 'EM') {
            # Randomize initial conditions if not given
            if (is.null(init)) {
                d <- nrow(u)
                init <- replicate(num.restart,
                                  runif(1 + d + 1 + d,
                                        min = c(0, rep(-1, d), 0, rep(-1, d)),
                                        max = c(1, rep( 1, d), 1, rep( 1, d))),
                                  simplify = F)
            } else {
                if (!is.list(init)) # learnLDS expects init is a list
                    init <- list(init)
            }
            ans <- LDS_EM_restart(y2, u, v, init, niter, tol, return.init = FALSE, parallel = FALSE)
            if (reupdate) {
                ans$fit <- Kalman_smoother(y, u, v, ans$theta)
            }
            ans
        } else {
            LDS_GA(y2, u, v, lambda, num.restart = num.restart, niter = niter, parallel = FALSE)
        }
    }

    if (parallel)  {
        nbCores <- detectCores()
        cl <- makeCluster(nbCores)
        registerDoParallel(cl)
        cv.results <- foreach(omit = Z) %dopar%
            one_CV(omit, y, u, v, method, init, num.restart, lambda, niter, tol)
        stopCluster(cl)
    } else {
        cv.results <- lapply(Z, function(omit)
            one_CV(omit, y, u, v, method, init, num.restart, lambda, niter, tol))
    }

    fit <- lapply(cv.results, '[[', 'fit')
    Y <- lapply(fit, '[[', 'Y') %>% lapply(function(v) as.vector(v)[obs.ind])
    Q <- lapply(Y, function(v) exp(v + mu))
    all.Q <- rbindlist(lapply(1:n.reps, function(i)
        data.table(rep = i, year = Qa$year, Q = Q[[i]])))
    cvQ <- rbindlist(lapply(1:n.reps, function(i) {
        ind <- Z2[[i]]
        data.table(rep = i,
                   year = ind + Qa$year[1] - 1,
                   cvQ = Q[[i]][ind],
                   cvObs = Qa$Qa[ind])
    }))
    metrics.dist <- rbindlist(lapply(1:n.reps,
                                     function(i) calculate_metrics(Q[[i]], Qa$Qa, Z2[[i]])))

    return(list(metrics.dist = metrics.dist,
                cvQ = cvQ,
                Z2 = Z2))
}
