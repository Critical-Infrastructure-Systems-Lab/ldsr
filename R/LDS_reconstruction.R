learnLDS_restart <- function(y, u, v, mu, init = NULL, num.restart,
                             niter = 1000, tol = 1e-5, return.init = FALSE,
                             parallel = FALSE) {

    # Learn multiple LDS models with randomized initial conditions (specified by num.restart)
    # and select the one with the highest likelihood
    # Use a pecified grid of initial conditions instead of randomized ones if provided

    if (is.null(init))
        init <- replicate(num.restart,
                          runif(4, min = c(0, -1, 0, -1), max = c(0.9, 1, 1, 1)),
                          simplify = F)

    if (parallel) {
        models <- foreach(init.val = init) %dopar%
            learnLDS(y, u, v, mu, init.val, niter, tol)
    } else {
        models <- lapply(init, function(init.val)
            learnLDS(y, u, v, mu, init.val, niter, tol))
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
        cv.results <- foreach(omit = Z) %dopar% {
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
