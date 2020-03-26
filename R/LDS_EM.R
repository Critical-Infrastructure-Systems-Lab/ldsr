#' Learn LDS model with multiple initial conditions
#'
#' This is the backend computation for [LDS_reconstruction].
#' @inheritParams LDS_EM
#' @param niter Maximum number of iterations, default 1000
#' @param tol Tolerance for likelihood convergence, default 1e-5. Note that the log-likelihood is normalized by dividing by the number of observations.
#' @param return.init Indicate whether the initial condition that results in the highest
#' @return a list as produced by [LDS_EM]. If return.init is true, a vector of initial condition is included in the list as well.
#' @export
LDS_EM_restart <- function(y, u, v, init, niter = 1000, tol = 1e-5, return.init = TRUE) {

    # To prevent global variable error in R CMD check
    init.val <- NULL
    models <- foreach(init.val = init, .packages = 'ldsr') %dopar% LDS_EM(y, u, v, init.val, niter, tol)

    # Select the model with highest likelihood
    # Only select models with C > 0 for physical interpretation (if possible)
    liks <- sapply(models, '[[', 'lik')
    all.C <- lapply(lapply(models, '[[', 'theta'), '[[', 'C')
    pos.C <- which(all.C > 0) # Positions of positive C
    if (length(pos.C) > 0) {
        max.ind <- which(liks == max(liks[pos.C]))
    } else {
        max.ind <- which.max(liks)
    }
    ans <- models[[max.ind]]
    if (return.init) ans$init <- init[[max.ind]]

    ans
}
