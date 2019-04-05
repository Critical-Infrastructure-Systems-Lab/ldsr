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

#' Learn LDS model with multiple initial conditions. This function can be called directly or wrapped in [LDS_ensemble]
#'
#' The initial conditions can either be randomized (specifiled by num.restarts) or provided beforehand.
#' @param Qa Observations: a data.table of annual streamflow with at least two columns: year and Qa.
#' @inheritParams LDS_EM_restart
#' @param start.year Starting year of the paleo period. `start.year + ncol(u) - 1` will determine the last year of the study horizon, which must be greater than or equal to the last year in `Qa`.
#' @param method Either EM or GA (note: GA is experimental)
#' @param trans Transformatin of flow before model fitting, can be either "log", "boxcox" or "none".
#' @param num.restarts if init is not given then num.restarts must be provided. In this case the function
#' will randomize the initial value by sampling uniformly within the range for each parameters
#' (A in \[0, 1\], B in \[-1, 1\], C in \[0, 1\] and D in \[-1, 1\]).
#' @param ub Upper bounds, a vector whose length is the number of parameters
#' @param lb Lower bounds
#' @param return.raw If TRUE, state and streamflow estimates without measurement updates will be returned.
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
LDS_reconstruction <- function(Qa, u, v, start.year, method = 'EM', trans = 'log',
                               init = NULL, num.restarts = 100, return.init = TRUE,
                               lambda = 1, ub, lb, num.islands = 4, pop.per.island = 250,
                               niter = 1000, tol = 1e-5, return.raw = FALSE,
                               parallel = TRUE, all.cores = FALSE) {

  if (method != 'EM' && (is.null(ub) || is.null(lb)))
    stop("For GA and BFGS methods, upper and lower bounds of parameters must be provided.")

  if (trans == 'log')
    Q.trans <- log(Qa$Qa)
  else if (trans == 'boxcox') {
    lambda <- car::powerTransform(Qa$Qa ~ 1)$roundlam %>% 'names<-'('lambda')
    Q.trans <- car::bcPower(Qa$Qa, lambda)
  } else if (trans == 'none')
    Q.trans <- Qa$Qa
  else stop('Accepted transformations are "log", "boxcox" and "none"')

  if (!identical(dim(u), dim(v))) stop('Dimensions of u and v must be the same.')

  N <- ncol(u)
  end.year <- start.year + N - 1
  if (end.year < Qa[.N, year])
    stop('The last year of u is earlier than the last year of the instrumental period.')

  # Attach NA and make the y matrix
  mu <- mean(Q.trans, na.rm = TRUE)
  y <- t(c(rep(NA, Qa[1, year] - start.year), # Before the instrumental period
           Q.trans - mu,     # Instrumental period
           rep(NA, end.year - Qa[.N, year])))

  results <-
    switch(method,
           EM = {
             init <- make_init(init, nrow(u), num.restarts)
             # To avoid unnecessary overhead, only run in parallel mode if the init list is long enough
             if (parallel && length(init) > 10) {
               LDS_EM_restart(y, u, v, init, niter, tol, return.init, parallel = TRUE, all.cores = all.cores)
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

  tidy_rec <- function(X, V, Y, theta, trans, years) {
    CI.X <- 1.96*sqrt(V)
    Xl <- X - CI.X # Lower range for X
    Xu <- X + CI.X # Upper range for X

    CI.Y <- 1.96 * (as.vector(theta$C) * V * as.vector(theta$C) + as.vector(theta$R))
    Yl <- Y - CI.Y # Lower range for Y
    Yu <- Y + CI.Y # Upper range for Y

    # Transform Y to Q and put in a data.table
    X.out <- data.table(year = years, X, Xl, Xu)
    Q.out <- switch(trans,
                    log = data.table(Q = exp(Y),
                                     Ql = exp(Yl),
                                     Qu = exp(Yu)),
                    boxcox = if (lambda == 0)
                      data.table(Q = exp(Y),
                                 Ql = exp(Yl),
                                 Qu = exp(Yu))
                    else
                      data.table(Q = (Y*lambda + 1)^(1/lambda),
                                 Ql = (Yl*lambda + 1)^(1/lambda),
                                 Qu = (Yu*lambda + 1)^(1/lambda)),
                    # no transformation
                    data.table(Q = Y,
                               Ql = Yl,
                               Qu = Yu)
    )
    cbind(X.out, Q.out)
  }
  # Construct 95% confidence intervals and return
  years <- start.year : end.year
  with(results, {

    rec <- tidy_rec(as.vector(fit$X),
                    as.vector(fit$V),
                    as.vector(fit$Y + mu),
                    theta, trans, years)
    ans <- list(rec = rec, theta = theta, lik = lik)
    if (return.raw) {
      # Return time series without measurement updates.
      # Note that we are unable to estimate confidence interval in this case
      X2 <- rep(0, N)
      V2 <- rep(0, N)
      V2[1] <- fit$V[1]
      for (t in 2:N) {
        X2[t] <- theta$A %*% X2[t - 1] + theta$B %*% u[, t - 1]
        V2[t] <- theta$A*V2[t]*theta$A + theta$Q
      }
      Y2 <- as.vector(theta$C) * X2 + as.vector(theta$D %*% u) + mu
      ans$rec2 <- tidy_rec(X2, V2, Y2, theta, trans, years)
    }

    if (return.init) ans$init <- results$init
    if (method != 'EM') ans$pl <- results$pl
    if (trans == 'boxcox') ans$lambda <- lambda
    ans
  })
}

#' Ensemble reconstruction. This is a wrapper for [LDS_reconstruction].
#'
#' @param u.list List of u matrices
#' @param v.list List of v matrices#'
#' @inheritParams LDS_reconstruction
#' @return A list of models in the ensemble
#' @export
LDS_ensemble <- function(Qa, u.list, v.list, start.year, method = 'EM', trans = 'log',
                         init = NULL, num.restarts = 100,
                         lambda = 1, ub = NULL, lb = NULL, num.islands = 4, pop.per.island = 250,
                         niter = 1000, tol = 1e-5, return.raw = FALSE,
                         parallel = TRUE, all.cores = FALSE) {

  if (length(u.list) != length(v.list))
    stop("Length of u and v lists must be the same.")

  if (parallel) {
    nbCores <- detectCores() - { if (all.cores) 0 else 1 }
    cl <- makeCluster(nbCores)
    registerDoParallel(cl)
    ensemble.full <- foreach(i = 1:length(u.list)) %dopar%
      LDS_reconstruction(Qa, u.list[[i]], v.list[[i]], start.year, method, trans,
                         init, num.restarts, return.init = FALSE,
                         lambda, ub, lb, num.islands, pop.per.island,
                         niter, tol, return.raw, parallel = FALSE)
    stopCluster(cl)
  } else {
    ensemble.full <- mapply(LDS_reconstruction, u = u.list, v = v.list,
                            MoreArgs = list(
                              Qa = Qa, start.year = start.year, method = method, trans = trans,
                              init = init, num.restarts = num.restarts, return.init = FALSE,
                              lambda = lambda, ub = ub, lb = lb,
                              num.islands = num.islands, pop.per.island = pop.per.island,
                              niter = niter, tol = tol, return.raw = return.raw,
                              parallel = FALSE
                            ),
                            SIMPLIFY = FALSE)
  }

  ensemble <- ensemble.full %>%
    lapply( '[[', 'rec') %>%
    rbindlist() %>%
    .[, member := 1:.N, by = year]

  rec <- ensemble[, .(Q = mean(Q)), by = year]
  Xstd <- ensemble[, .(year, X = X / sd(X)), by = member
              ][, .(X = mean(X)), by = year
              ][, X]
  rec[, X := Xstd]
  list(rec = rec[],
       ensemble = ensemble,
       theta = lapply(ensemble.full, '[[', 'theta'))

}

#' Cross validate LDS model. This is a wrapper for [LDS_reconstruction]
#'
#' @inheritParams LDS_reconstruction
#' @param k Numer of data points to be left out in each cross validation-run.
#' @param CV.reps Number of cross-validation runs.
#' @param Z A list of CV.reps elements, each is a vector of length k. See Details.
#' @param metrics.space Whether performance metrics are calculated in the 'original' or 'transformed' space, or 'both'.
#' @details Allows different experimental setups:
#'   * init: if given, all cross validation runs start with the same init, otherwise each cross validation run is learned using randomized restarts.
#'   * Z: If given, cross validation will be run on these points (useful when comparing LDS with another reconstruction method); otherwise, randomized cross validation points will be created.
#' @export
cvLDS <- function(Qa, u, v, start.year, method = 'EM', trans = 'log',
                  k, CV.reps = 100, Z = NULL,
                  init = NULL, num.restarts = 100,
                  lambda = 1, ub, lb, num.islands = 4, pop.per.island = 100,
                  niter = 1000, tol = 1e-5, parallel = TRUE, all.cores = FALSE) {

  obs.ind <- which(!is.na(Qa$Qa))
  n.obs <- length(obs.ind)
  if (missing(k)) k <- ceiling(n.obs / 10)
  # Z: index in instrumental period
  if (is.null(Z)) {
    Z <- replicate(CV.reps, sample(obs.ind, k), simplify = FALSE)
  } else {
    CV.reps <- length(Z)
  }
  # ub and lb must be passed through to one_Cv, otherwise parallel run will produce an erroo.
  if (method == 'EM') ub <- lb <- NULL

  one_CV <- function(omit, Qa, u, v, start.year, method, trans, init, num.restarts,
                     lambda, ub, lb, num.islands, pop.per.island, niter, tol) {

    # Don't change Qa so that we can still calculate metrics later
    Qa2 <- copy(Qa) %>% .[omit, Qa := NA]
    # Parallel is run at the outer loop, i.e. for each cross-validation run
    ans <- LDS_reconstruction(Qa2, u, v, start.year, method, trans, init, num.restarts, return.init = FALSE,
                              lambda, ub, lb, num.islands, pop.per.island, niter, tol,
                              parallel = FALSE)

    Qa.hat <- ans$rec[year %in% Qa$year, Q]
    metric <- calculate_metrics(Qa.hat, Qa$Qa, omit)

    list(metric = metric, Ycv = Qa.hat)
  }

  if (parallel)  {
    nbCores <- detectCores() - { if (all.cores) 0 else 1 }
    cl <- makeCluster(nbCores)
    registerDoParallel(cl)
    cv.results <- foreach(omit = Z) %dopar%
      one_CV(omit, Qa, u, v, start.year, method, trans, init, num.restarts,
             lambda, ub, lb, num.islands, pop.per.island, niter, tol)
    stopCluster(cl)
  } else {
    cv.results <- lapply(Z, function(omit)
      one_CV(omit, Qa, u, v, start.year, method, trans, init, num.restarts,
             lambda, ub, lb, num.islands, pop.per.island, niter, tol))
  }

  metrics.dist <- rbindlist(lapply(cv.results, '[[', 'metric'))
  Ycv <- as.data.table(sapply(cv.results, '[[', 'Ycv'))
  Ycv$year <- Qa$year
  Ycv <- melt(Ycv, id.vars = 'year', variable.name = 'rep', value.name = 'Qa')

  return(list(metrics = colMeans(metrics.dist[, 1:5]),
              metrics.dist = metrics.dist,
              Ycv = Ycv,
              Z = Z))
}

#' Cross validatdion for ensemble reconstruction. This is a wrapper for [cvLDS].
#'
#' @param u.list List of u matrices
#' @param v.list List of v matrices#'
#' @inheritParams cvLDS
#' @return same as cvLDS
#' @export
cvLDS_ensemble <- function(Qa, u.list, v.list, start.year, method = 'EM', trans = 'log',
                           k, CV.reps = 100, Z = NULL,
                           init = NULL, num.restarts = 100,
                           lambda = 1, ub, lb, num.islands = 4, pop.per.island = 100,
                           niter = 1000, tol = 1e-5, parallel = TRUE, all.cores = FALSE) {

  if (length(u.list) != length(v.list))
    stop("Length of u and v lists must be the same.")

  obs.ind <- which(!is.na(Qa$Qa))
  n.obs <- length(obs.ind)
  if (missing(k)) k <- ceiling(n.obs / 10)
  # Z: index in instrumental period
  if (is.null(Z)) {
    Z <- replicate(CV.reps, sample(obs.ind, k), simplify = FALSE)
  } else {
    CV.reps <- length(Z)
  }
  # ub and lb must be passed through to one_Cv, otherwise parallel run will produce an erroo.
  if (method == 'EM') ub <- lb <- NULL

  one_CV_ensemble <- function(omit, Qa, u.list, v.list, start.year, method, trans, init, num.restarts,
                              lambda, ub, lb, num.islands, pop.per.island, niter, tol) {

    # Don't change Qa so that we can still calculate metrics later
    Qa2 <- copy(Qa) %>% .[omit, Qa := NA]
    # Parallel is run at the outer loop, i.e. for each cross-validation run
    ans <- LDS_ensemble(Qa2, u.list, v.list, start.year, method, trans, init, num.restarts,
                        lambda, ub, lb, num.islands, pop.per.island, niter, tol,
                        return.raw = FALSE, parallel = FALSE)

    Qa.hat <- ans$rec[year %in% Qa$year, Q]
    metric <- calculate_metrics(Qa.hat, Qa$Qa, omit)

    list(metric = metric, Ycv = Qa.hat)
  }

  if (parallel)  {
    nbCores <- detectCores() - { if (all.cores) 0 else 1 }
    cl <- makeCluster(nbCores)
    registerDoParallel(cl)
    cv.results <- foreach(omit = Z) %dopar%
      one_CV_ensemble(omit, Qa, u.list, v.list, start.year, method, trans, init, num.restarts,
                      lambda, ub, lb, num.islands, pop.per.island, niter, tol)
    stopCluster(cl)
  } else {
    cv.results <- lapply(Z, function(omit)
      one_CV_ensemble(omit, Qa, u.list, v.list, start.year, method, trans, init, num.restarts,
                      lambda, ub, lb, num.islands, pop.per.island, niter, tol))
  }

  metrics.dist <- rbindlist(lapply(cv.results, '[[', 'metric'))
  Ycv <- as.data.table(sapply(cv.results, '[[', 'Ycv'))
  Ycv$year <- Qa$year
  Ycv <- melt(Ycv, id.vars = 'year', variable.name = 'rep', value.name = 'Qa')

  return(list(metrics = colMeans(metrics.dist[, 1:5]),
              metrics.dist = metrics.dist,
              Ycv = Ycv,
              Z = Z))
}
