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
#' @param trans Transformatin of flow before model fitting, can be either "log", "boxcox" or "none".
#' @param num.restarts if init is not given then num.restarts must be provided. In this case the function
#' will randomize the initial value by sampling uniformly within the range for each parameters
#' (A in \[0, 1\], B in \[-1, 1\], C in \[0, 1\] and D in \[-1, 1\]).
#' @param ub Upper bounds, a vector whose length is the number of parameters
#' @param lb Lower bounds
#' @param return.raw If TRUE,
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
LDS_reconstruction <- function(Qa, u, v, method = 'EM', trans = 'log',
                               init = NULL, num.restarts = 100, return.init = TRUE,
                               lambda = 1, ub, lb, num.islands = 4, pop.per.island = 250,
                               niter = 1000, tol = 1e-5, return.raw = FALSE,
                               parallel = TRUE) {

    # Attach NA and make the y matrix
    if (trans == 'log')
        Q.trans <- log(Qa$Qa)
    else if (trans == 'boxcox') {
        lambda <- car::powerTransform(Qa$Qa ~ 1)$roundlam %>% 'names<-'('lambda')
        Q.trans <- car::bcPower(Qa$Qa, lambda)
    } else if (trans == 'none')
        Q.trans <- Qa$Qa
    else stop('Accepted transformations are "log", "boxcox" and "none"')

    mu <- mean(Q.trans, na.rm = TRUE)
    N <- ncol(u)
    n.paleo <- N - nrow(Qa) # Number of years in the paleo period
    y <- t(c(rep(NA, n.paleo), Q.trans - mu))

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
        rec <- switch(trans,
                      log = data.table(year = (Qa[1, year] - n.paleo):(Qa[.N, year]),
                                       X, Xl, Xu,
                                       Q = exp(Y),
                                       Ql = exp(Yl),
                                       Qu = exp(Yu)),
                      boxcox = if (lambda == 0)
                          data.table(year = (Qa[1, year] - n.paleo):(Qa[.N, year]),
                                     X, Xl, Xu,
                                     Q = exp(Y),
                                     Ql = exp(Yl),
                                     Qu = exp(Yu))
                      else
                          data.table(year = (Qa[1, year] - n.paleo):(Qa[.N, year]),
                                     X, Xl, Xu,
                                     Q = (Y*lambda + 1)^(1/lambda),
                                     Ql = (Yl*lambda + 1)^(1/lambda),
                                     Qu = (Yu*lambda + 1)^(1/lambda)),
                      data.table(year = (Qa[1, year] - n.paleo):(Qa[.N, year]),
                                 X, Xl, Xu,
                                 Q = Y,
                                 Ql = Yl,
                                 Qu = Yu)
        )

        ans <- list(rec = rec, theta = theta, lik = lik)

        # Return time series without measurement updates.
        # Note that we are unable to estimate confidence interval in this case
        if (return.raw) {
            X2 <- rep(0, N)
            V2 <- rep(0, N)
            V2[1] <- V[1]
            for (t in 2:N) {
                X2[t] <- theta$A %*% X2[t - 1] + theta$B %*% u[, t - 1]
                V2[t] <- A*V2[t]*A + theta$Q
            }
            CI.X2 <- 1.96*sqrt(V2)
            Xl2 <- X2 - CI.X2 # Lower range for X
            Xu2 <- X2 + CI.X2 # Upper range for X
            Y2 <- as.vector(theta$C) * X2 + as.vector(theta$D %*% u) + mu
            CI.Y2 <- 1.96 * (as.vector(theta$C) * V2 * as.vector(theta$C) + as.vector(theta$R))
            Yl2 <- Y2 - CI.Y2 # Lower range for Y
            Yu2 <- Y2 + CI.Y2 # Upper range for Y
            # Transform Y to Q and put in a data.table
            rec2 <- switch(trans,
                           log = data.table(year = (Qa[1, year] - n.paleo):(Qa[.N, year]),
                                            X = X2, Xl = Xl2, Xu = Xu2,
                                            Q = exp(Y2),
                                            Ql = exp(Yl2),
                                            Qu = exp(Yu2)),
                           boxcox = if (lambda == 0)
                               data.table(year = (Qa[1, year] - n.paleo):(Qa[.N, year]),
                                          X = X2, Xl = Xl2, Xu = Xu2,
                                          Q = exp(Y2),
                                          Ql = exp(Yl2),
                                          Qu = exp(Yu2))
                           else
                               data.table(year = (Qa[1, year] - n.paleo):(Qa[.N, year]),
                                          X = X2, Xl = Xl2, Xu = Xu2,
                                          Q = (Y2*lambda + 1)^(1/lambda),
                                          Ql = (Yl2*lambda + 1)^(1/lambda),
                                          Qu = (Yu2*lambda + 1)^(1/lambda)),
                           data.table(year = (Qa[1, year] - n.paleo):(Qa[.N, year]),
                                      X = X2, Xl = Xl2, Xu = Xu2,
                                      Q = Y2,
                                      Ql = Yl2,
                                      Qu = Yu2)
            )
            ans$rec2 <- rec2
        }

        if (return.init) ans$init <- results$init
        if (method != 'EM') ans$pl <- results$pl
        if (trans == 'boxcox') ans$lambda <- lambda
        ans
    })
}

#' Cross validate LDS model
#'
#' @inheritParams LDS_reconstruction
#' @param k Numer of data points to be left out in each cross validation-run.
#' @param CV.reps Number of cross-validation runs.
#' @param Z A list of CV.reps elements, each is a vector of length k.See Details.
#' @param metrics.space Whether performance metrics are calculated in the 'original' or 'transformed' space, or 'both'.
#' @details Allows different experimental setups:
#'   * init: if given, all cross validation runs start with the same init, otherwise each cross validation run is learned using randomized restarts.
#'   * Z: If given, cross validation will be run on these points (useful when comparing LDS with another reconstruction method); otherwise, randomized cross validation points will be created.
#' @export
cvLDS <- function(Qa, u, v, method = 'EM', trans = 'log',
                  k, CV.reps = 100, Z = NULL, metrics.space = 'original',
                  init = NULL, num.restarts = 100,
                  lambda = 1, ub, lb, num.islands = 4, pop.per.island = 100,
                  niter = 1000, tol = 1e-5, parallel = TRUE) {

    n.paleo <- ncol(u) - nrow(Qa) # Number of years in the paleo period
    inst.period <- (n.paleo + 1):ncol(u)
    if (trans == 'log')
        Q.trans <- log(Qa$Qa)
    else if (trans == 'boxcox') {
        lambda <- car::powerTransform(Qa$Qa ~ 1)$roundlam %>% 'names<-'('lambda')
        Q.trans <- car::bcPower(Qa$Qa, lambda)
    } else if (trans == 'none')
        Q.trans <- Qa$Qa
    else stop('Accepted transformations are "log", "boxcox" and "none"')

    if (!(metrics.space) %in% c('original', 'transformed', 'both'))
        stop('Metrics space undefined. Must be "original", "transformed", or "both"')

    mu <- mean(Q.trans, na.rm = TRUE)
    y <- t(c(rep(NA, n.paleo), Q.trans - mu))

    obs.ind <- which(!is.na(Qa$Qa))
    n.obs <- length(obs.ind)

    if (missing(k)) k <- ceiling(n.obs / 10)

    # Z: index in instrumental period
    if (is.null(Z)) Z <- replicate(CV.reps, sample(obs.ind, k), simplify = F)

    # ub and lb must be passed through to one_Cv, otherwise parallel run will produce an erroo.
    # TODO: change to ... in argument list. This is a temporary fix. With (...) we can past
    # different arguments to different calls of one_CV based on context.
    if (method == 'EM') ub <- lb <- NULL

    one_CV <- function(omit, method, y, u, v, init, num.restarts,
                       lambda, ub, lb, num.islands, pop.per.island, niter, tol) {

        y2 <- y
        y2[omit + n.paleo] <- NA
        # Parallel is run at the outer loop, i.e. for each cross-validation run
        ans <- switch(method,
                      EM = {
                          init <- make_init(init, nrow(u), num.restarts)
                          LDS_EM_restart(y2, u, v, init, niter, tol, return.init = FALSE, parallel = FALSE)
                      },
                      GA = {
                          if (missing(ub) || missing(lb)) stop("Upper and lower bounds must be provided.")
                          LDS_GA(y2, u, v, lambda, ub, lb, num.islands, pop.per.island, niter, parallel = FALSE)
                      },
                      BFGS = {
                          if (missing(ub) || missing(lb)) stop("Upper and lower bounds must be provided.")
                          LDS_BFGS(y, u, v, lambda, ub, lb, num.restarts, parallel)
                      },
                      stop("Method undefined. It has to be either EM, GA or BFGS.")
        )
        Y <- as.vector(ans$fit$Y + mu)[inst.period]
        Qa.hat <- switch(trans,
                         log = exp(Y),
                         boxcox = if (lambda == 0) exp(Y) else (Y*lambda + 1)^(1/lambda),
                         Y)
        metric <- switch(metrics.space,
                         original = calculate_metrics(Qa.hat, Qa$Qa, omit) %>% .[, space := 'original'],
                         trans = calculate_metrics(Y, Q.trans, omit) %>% .[, space := 'transformed'],
                         both = rbind(
                             calculate_metrics(Qa.hat, Qa$Qa, omit) %>% .[, space := 'original'],
                             calculate_metrics(Y, Q.trans, omit) %>% .[, space := 'transformed']
                         ))
        list(metric = metric, Ycv = Qa.hat)
    }

    if (parallel)  {
        nbCores <- detectCores() - 1
        cl <- makeCluster(nbCores)
        registerDoParallel(cl)
        cv.results <- foreach(omit = Z) %dopar%
            one_CV(omit, method, y, u, v, init, num.restarts,
                   lambda, ub, lb, num.islands, pop.per.island, niter, tol)
        stopCluster(cl)
    } else {
        cv.results <- lapply(Z, function(omit)
            one_CV(omit, method, y, u, v, init, num.restarts,
                   lambda, ub, lb, num.islands, pop.per.island, niter, tol))
    }

    metrics.dist <- rbindlist(lapply(cv.results, '[[', 'metric'))
    Ycv <- as.data.table(sapply(cv.results, '[[', 'Ycv'))
    Ycv$year <- Qa$year
    Ycv <- melt(Ycv, id.vars = 'year', variable.name = 'rep', value.name = 'Qa')

    return(list(metrics = {
        if (metrics.space != 'both') colMeans(metrics.dist[, 1:5])
        else metrics.dist[, lapply(.SD, mean), by = space]
    },
    metrics.dist = metrics.dist,
    Ycv = Ycv,
    Z = Z))
}
