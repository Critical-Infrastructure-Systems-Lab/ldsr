#' Make a list of initial parameter values for the search
#'
#' If init is a vector, make it a list
#' If init is NULL, randomize it
#' @param init The init provided by the caller, can be NULL
#' @param d Dimension of input
#' @param num.restarts Number of randomized initial conditions
#' @return A list of iniial conditions
make_init <- function(init, d, num.restarts) {

    if (is.null(init)) {
        replicate(num.restarts,
                  runif(1 + d + 1 + d,
                        min = c(0, rep(-1, d), 0, rep(-1, d)),
                        max = c(1, rep( 1, d), 1, rep( 1, d))),
                        simplify = F)
    } else {
        # learnLDS expects init is a list
        if (!is.list(init)) list(init) else init
    }
}

#' Learn LDS model with multiple initial conditions
#'
#' The initial conditions can either be randomized (specifiled by num.restarts) or provided beforehand.
#' @param Qa Observations: a data.table of annual streamflow with at least two columns: year and Qa.
#' @inheritParams LDS_EM_restart
#' @param method Either EM or GA (note: GA is experimental)
#' @param num.restarts if init is not given then num.restarts must be provided. In this case the function
#' will randomize the initial value by sampling uniformly within the range for each parameters
#' (A in \[0, 1\], B in \[-1, 1\], C in \[0, 1\] and D in \[-1, 1\]).
#' @param ub Upper bounds, a vector whose length is the number of parameters
#' @param lb Lower bounds
#' @return A list of the following elements
#' * rec: reconstruction results, a data.table with the following columns
#'     - year: calculated from Qa and the length of u
#'     - X: the estimated hidden state
#'     - Xl, Xu: lower and upper range for the 95\% confidence interval of X
#'     - Q: the reconstructed streamflow
#'     - Ql, Qu: lower and upper range for the 95\% confidence interval of Q
#' * theta: model parameters
#' * lik: maximum likelihood
#' * init: the initial condition that resulted in the maximum likelihood (if )
#' @export
LDS_reconstruction <- function(Qa, u, v, method = 'EM',
                               init = NULL, num.restarts = 100, return.init = TRUE,
                               lambda = 1, ub, lb, num.islands = 4, pop.per.island = 250,
                               niter = 1000, tol = 1e-5,
                               parallel = TRUE) {

    # Attach NA and make the y matrix
    mu <- mean(log(Qa$Qa), na.rm = TRUE)
    n.paleo <- ncol(u) - nrow(Qa) # Number of years in the paleo period
    y <- t(c(rep(NA, n.paleo), log(Qa$Qa) - mu))

    results <- switch(method,
        EM = {
            init <- make_init(init, nrow(u), num.restarts)
            # To avoid unnecessary overhead, only run in parallel mode if the init list is long enough
            if (parallel && length(init) > 10) {
                LDS_EM_restart(y, u, v, init, niter, tol, return.init, parallel = TRUE)
            } else {
                if (parallel)
                    warning('Initial condition list is short, LDS_EM_restart() is run in sequential mode.')
                LDS_EM_restart(y, u, v, init, niter, tol, return.init, parallel = FALSE)
            }
        },
        GA = {
            if (missing(ub) || missing(lb)) stop("Upper and lower bounds must be provided.")
            LDS_GA(y, u, v, lambda, ub, lb, num.islands, pop.per.island, niter, parallel)
        },
        BFGS = {
            if (missing(ub) || missing(lb)) stop("Upper and lower bounds must be provided.")
            LDS_BFGS(y, u, v, lambda, ub, lb, num.restarts, parallel)
        },
        {
            stop("Method undefined. It has to be either EM, GA or BFGS.")
        }
    )

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
        if (method != 'EM') ans$pl <- results$pl
        ans
    })
}

#' Cross validate LDS model
#'
#' @inheritParams LDS_reconstruction
#' @param k Numer of data points to be left out in each cross validation-run.
#' @param CV.reps Number of cross-validation runs.
#' @param Z A list of CV.reps elements, each is a vector of length k.See Details.
#' @details Allows different experimental setups:
#'   * init: if given, all cross validation runs start with the same init, otherwise each cross validation run is learned using randomized restarts.
#'   * Z: If given, cross validation will be run on these points (useful when comparing LDS with another reconstruction method); otherwise, randomized cross validation points will be created.
#' @export
cvLDS <- function(Qa, u, v, method = 'EM',
                  k, CV.reps = 100, Z = NULL,
                  init = NULL, num.restarts = 100,
                  lambda = 1, ub, lb, num.islands = 4, pop.per.island = 100,
                  niter = 1000, tol = 1e-5, parallel = TRUE) {

    mu <- mean(log(Qa$Qa), na.rm = TRUE)
    n.paleo <- ncol(u) - nrow(Qa) # Number of years in the paleo period
    y <- t(c(rep(NA, n.paleo), log(Qa$Qa) - mu))

    obs.ind <- which(!is.na(y))
    n.obs <- length(obs.ind)

    if (missing(k)) k <- ceiling(n.obs / 10)
    if (is.null(Z)) {
        Z <- replicate(CV.reps, sample(obs.ind, k), simplify = F)
        Z2 <- lapply(Z, '-', n.paleo)
    } else {
        Z2 <- Z
        Z <- lapply(Z2, '+', n.paleo)
    }

    one_CV <- function(omit, method, y, u, v, init, num.restarts,
                       lambda, ub, lb, num.islands, pop.per.island, niter, tol) {
        y2 <- y
        y2[omit] <- NA
        # Parallel is run at the outer loop, i.e. for each cross-validation run
        switch(method,
            EM = {
                init <- make_init(init, nrow(u), num.restarts)
                LDS_EM_restart(y2, u, v, init, niter, tol, return.init = FALSE, parallel = FALSE)
            },
            GA = {
                if (missing(ub) || missing(lb)) stop("Upper and lower bounds must be provided.")
                LDS_GA(y2, u, v, lambda, ub, lb, num.islands, pop.per.island, niter, parallel = FALSE)
            },
            BFGS = {
                stop("BFGS has not been implemented. Please choose either EM or GA for now.")
            },
            {
                stop("Method undefined. It has to be either EM, GA or BFGS.")
            }
        )
    }

    if (parallel)  {
        nbCores <- detectCores()
        cl <- makeCluster(nbCores)
        registerDoParallel(cl)
        cv.results <- foreach(omit = Z) %dopar%
            one_CV(omit, method, y, u, v, init, num.restarts, lambda, ub, lb, num.islands, pop.per.island, niter, tol)
        stopCluster(cl)
    } else {
        cv.results <- lapply(Z, function(omit)
            one_CV(omit, method, y, u, v, init, num.restarts, lambda, ub, lb, num.islands, pop.per.island, niter, tol))
    }

    fit <- lapply(cv.results, '[[', 'fit')
    Y <- lapply(fit, '[[', 'Y') %>% lapply(function(v) as.vector(v)[obs.ind])
    Q <- lapply(Y, function(v) exp(v + mu))
    all.Q <- rbindlist(lapply(1:CV.reps, function(i)
        data.table(rep = i, year = Qa$year, Q = Q[[i]])))
    cvQ <- rbindlist(lapply(1:CV.reps, function(i) {
        ind <- Z2[[i]]
        data.table(rep = i,
                   year = ind + Qa$year[1] - 1,
                   cvQ = Q[[i]][ind],
                   cvObs = Qa$Qa[ind])
    }))
    metrics.dist <- rbindlist(lapply(1:CV.reps,
                                     function(i) calculate_metrics(Q[[i]], Qa$Qa, Z2[[i]])))

    return(list(metrics.dist = metrics.dist,
                cvQ = cvQ,
                Z2 = Z2))
}
