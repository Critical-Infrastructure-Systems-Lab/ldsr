#' Learn LDS model with multiple initial conditions
#'
#' The initial conditions can either be randomized (specifiled by num.restart) or provided beforehand.
#' @param Qa Observations: a data frame of annual streamflow with at least two columns: year and Qa.
#' @param u Input matrix for the state equation (m_u rows, T columns)
#' @param v Input matrix for the output equation (m_v rows, T columns)
#' @param init A list of initial conditions, each element is a vector of length 4, the initial values for A, B, C and D. The initial values for Q and R are always 1, and mu_1 is 0 and V_1 is 1.
#' @param num.restart if init is not given then num.restart must be provided. In this case the function will randomize the initial value by sampling uniformly within the range for each parameters (A in \[0, 1\], B in \[-1, 1\], C in \[0, 1\] and D in \[-1, 1\]).
#' @param niter Maximum number of iterations, default 1000
#' @param tol Tolerance for likelihood convergence, default 1e-5. Note that the log-likelihood is normalized by dividing by the number of observations.
#' @param return.init Indicate whether the initial condition that results in the highest log-likelihood is returned. Default is TRUE.
#' @param parallel If TRUE, the computation is done in parallel using all available cores (using the doParallel backend). If FALSE, the computation is done serially.
#' @export
learnLDS_restart <- function(Qa, u, v, init = NULL, num.restart,
                             niter = 1000, tol = 1e-5, return.init = TRUE,
                             parallel = FALSE) {

    if (is.null(init))
        init <- replicate(num.restart,
                          runif(4, min = c(0, -1, 0, -1), max = c(0.9, 1, 1, 1)),
                          simplify = F)
    # Attach NA and make the y matrix
    mu <- mean(log(Qa$Qa), na.rm = T)
    n.paleo <- ncol(u) - nrow(Qa) # Number of years in the paleo period
    y <- t(c(rep(NA, n.paleo), log(Qa$Qa) - mu))

    if (parallel) {
        nbCores <- detectCores()
        cl <- makeCluster(nbCores)
        registerDoParallel(cl)
        models <- foreach(init.val = init) %dopar% learnLDS(y, u, v, init.val, niter, tol)
        stopCluster(cl)
    } else {
        models <- lapply(init, function(init.val) learnLDS(y, u, v, init.val, niter, tol))
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

cvLDS <- function(Qa, u, v, init = NULL, num.restart = 20, k, n.reps = 100, niter = 1000, tol = 1e-5, Z = NULL, parallel = F) {

    # Qa: dataframe (year, Qa)
    # u : d x years matrix
    # allows different experimental setups:
    #   + init: if given, all cross validation runs start with the same init,
    #           otherwise each cross validation run is learned using randomized restarts
    #   + Z: a list of n.reps elements, each is a vector of length k. If given, cross validation will be run on these points
    #        (useful when comparing LDS with another reconstruction method),
    #        Otherwise randomized cross validation points will be created

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
        z2 <- Z
        Z <- lapply(Z2, '+', n.paleo)
    }

    if (parallel)  {
        cv.results <- foreach::foreach(omit = Z) %dopar% {
            y2 <- y
            y2[omit] <- NA
            if (is.null(init)) {
                learnLDS_restart(y2, u, v, mu, num.restart, niter, tol)
            } else {
                learnLDS(y2, u, v, mu, init, niter, tol)
            }
        }
    } else {
        cv.results <- lapply(Z, function(omit) {
            y2 <- y
            y2[omit] <- NA
            if (is.null(init)) {
                learnLDS_restart(y2, u, v, mu, num.restart, niter, tol)
            } else {
                learnLDS(y2, u, v, mu, init, niter, tol)
            }
        })
    }

    fit <- lapply(cv.results, '[[', 'fit')
    Y <- lapply(fit, '[[', 'Y') %>% lapply(function(v) as.vector(v)[obs.ind])
    Y2 <- lapply(fit, '[[', 'Y2') %>% lapply(function(v) as.vector(v)[obs.ind])
    Q <- lapply(Y, function(v) exp(v + mu))
    Q2 <- lapply(Y2, function(v) exp(v + mu))

    all.Q <- rbindlist(lapply(1:n.reps, function(i) data.table(rep = i,
                                                            year = Qa$year,
                                                            Q = Q[[i]],
                                                            Q2 = Q2[[i]])))
    cvQ <- rbindlist(lapply(1:n.reps, function(i) {

        ind <- Z2[[i]]
        data.table(rep = i,
                   year = ind + Qa$year[1] - 1,
                   cvQ = Q[[i]][ind],
                   cvQ2 = Q2[[i]][ind],
                   cvObs = Qa$Qa[ind])
    }
    ))
    metrics.dist <- rbindlist(lapply(1:n.reps,
                                     function(i) calculate_metrics(Q[[i]], Qa$Qa, Z2[[i]])))
    metrics.dist2 <- rbindlist(lapply(1:n.reps,
                                     function(i) calculate_metrics(Q2[[i]], Qa$Qa, Z2[[i]])))

    return(list(metrics.dist = metrics.dist,
                metrics.dist2 = metrics.dist2,
                cvQ = cvQ,
                Z2 = Z2))
}

plot_cv <- function(cv.result, case.name) {

    metrics.dist <- cv.result$metrics.dist
    metrics.dist2 <- cv.result$metrics.dist2
    all.Q <- cv.result$all.Q

    # TODO: fix this, maybe copy from another folder
    #cvQ <-

    metrics.dist[, rep := 1:.N]
    metrics.dist2[, rep := 1:.N]

    p <- ggplot(all.Q) +
        geom_line(aes(year, Q, colour = 'LDS'), size = 0.1) +
        geom_line(aes(year, Qa, colour = 'Inst'), data = Qa, size = 0.1) +
        geom_point(aes(year, cvQ, colour = 'LDS'), data = cvQ, size = 0.1) +
        geom_point(aes(year, cvObs, colour = 'Inst'), data = cvQ, size = 0.1) +
        geom_text(aes(x = 1922, y = 4000, label = paste0('CE = ', round(CE, 2))),
                  data = metrics.dist,
                  hjust = 'left') +
        facet_wrap(~rep)
    ggsave(paste0('cv_results_', case.name, 'with_update.pdf'),
                  p,
                  width = 420, height = 297, units = 'mm')

    q <- ggplot(all.Q) +
        geom_line(aes(year, Q2, colour = 'LDS'), size = 0.1) +
        geom_line(aes(year, Qa, colour = 'Inst'), data = Qa, size = 0.1) +
        geom_point(aes(year, cvQ2, colour = 'LDS'), data = cvQ, size = 0.1) +
        geom_point(aes(year, cvObs, colour = 'Inst'), data = cvQ, size = 0.1) +
        geom_text(aes(x = 1922, y = 4000, label = paste0('CE = ', round(CE, 2))),
                  data = metrics.dist,
                  hjust = 'left') +
        facet_wrap(~rep)
    ggsave(paste0('cv_results_', case.name, 'without_update.pdf'),
           q,
           width = 420, height = 297, units = 'mm')
}
