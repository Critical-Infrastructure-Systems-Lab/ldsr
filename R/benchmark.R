#' Reconstruction with Principal Component Regression
#'
#' Reconstruct streamflow using the given principal components and backward selection.
#' Since this is the benchmark for the LDS model, and for simplicity, learning and cross-validation are performed in one go, instead of with two separate functions as in the case of LDS.
#'
#' @inheritParams cvLDS
#' @param pc A data frame of principal components with first column called years (the years in) the
#' study period, followed by one column for each principal components.
#' @return A list of reconstruction and cross-validation results
#' * rec: reconstruction
#' * metrics: performance metrics
#' * metrics.dist: distribution of performance metrics across all cross-validation runs.
#' * Ycv: the predicted streamflow in each cross validation run
#' * Z: the cross-validation points
#' @export
PCR_reconstruction <- function(Qa, pc, k, CV.reps = 100, Z = NULL) {

    # Function to calculate cross validation metrics
    cv <- function(z) {
        # Leave-k-out cross validation with the indices of the k data points to be left out supplied in z
        # Returns a vector of performance metrics

        fit <- lm(log(Qa) ~ . , data = df2, subset = setdiff(1:N, z))
        Qa.hat <- exp(predict(fit, newdata = df))

        return(list(metric = calculate_metrics(Qa.hat, df2$Qa, z),
                    Ycv = Qa.hat))
    }

    # Main model  -------------------------------------------------------------------
    Qa <- as.data.table(Qa)
    pc <- as.data.table(pc)
    # Remove observed value after the end of the paleorecord and missing values
    y <- Qa[year <= max(pc$year) & (!is.na(Qa))]
    mu <- mean(y$Qa)

    # Important to remove year otherwise it'll be come a predictor
    years <- pc$year
    df <- merge(pc, Qa, by = 'year')
    df$year <- NULL
    pc$year <- NULL

    # Model fitting
    fit <- step(lm(log(Qa) ~ . , data = df), direction = 'backward', trace = 0)
    selected <- names(fit$model)[-1]   # First element is intercept
    if (length(selected) == 0) {
        warning("Backward stepwise returned empty model; use all variables instead.")
        fit <- lm(log(Qa) ~ . , data = df)
        selected <- names(fit$model)[-1]   # First element is intercept
    }
    rec <- data.table(exp(predict(fit, newdata = pc, interval = 'confidence')))
    colnames(rec) <- c('Q', 'Ql', 'Qu')
    rec$year <- years

    df2 <- df[, ..selected]
    df2$Qa <- df$Qa

    # Cross validation ------------------------------------------------------------
    N <- nrow(y)
    if (missing(k)) k <- ceiling(N/10) # Default 10%
    if (is.null(Z)) Z <- replicate(CV.reps, sample(1:N, size = k, replace = FALSE))
    if (is.list(Z)) Z <- matrix(unlist(Z), ncol = CV.reps)
    if (k > 1) { # Leave-k-out, k > 1
        cv_res <- apply(Z, 2, cv)
        metrics.dist <- rbindlist(lapply(cv_res, '[[', 'metric'))
        Ycv <- as.data.table(sapply(cv_res, '[[', 'Ycv'))
        Ycv$year <- Qa$year
        Ycv <- melt(Ycv, id.vars = 'year', variable.name = 'rep', value.name = 'Qa')
    } else {
        # Leave-one-out
        cv_res <- rbindlist(lapply(1:N, cv))
    }

    return(list(
        rec = rec,
        coeffs = fit$coefficients,
        sigma = summary(fit)$sigma,
        selected = selected,
        metrics.dist = metrics.dist,
        metrics = colMeans(metrics.dist),
        Ycv = Ycv,
        Z = t(Z)
    ))
}
