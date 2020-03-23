#' Reconstruction with Principal Component Regression
#'
#' Reconstruct streamflow using the given principal compone11nts and backward selection.
#' Since this is the benchmark for the LDS model, and for simplicity, learning and cross-validation are performed in one go, instead of with two separate functions as in the case of LDS.
#'
#' @param Qa Observations: a data.frame of annual streamflow with at least two columns: year and Qa.
#' @param pc A data.frame, one colulmn for each principal component
#' @param start.year Starting year of `pc`, i.e, the first year of the paleo period. `start.year + nrow(pc) - 1` will determine the last year of the study horizon, which must be greater than or equal to the last year in `Qa`.
#' @param transform Flow transformation, either "log", "boxcox" or "none".
#' @return A list of reconstruction results
#' * rec: reconstruction
#' * selected: a vector of selected principal components
#' @export
PCR_reconstruction <- function(Qa, pc, start.year, transform = 'log') {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)
  pc <- as.data.table(pc)
  end.year <- start.year + nrow(pc) - 1
  if (end.year < Qa[.N, year])
    stop('The last year of pc is earlier than the last year of the instrumental period.')

  y <- Qa$Qa
  if (transform == 'log') {
    y <- log(y)
  } else if (transform == 'boxcox') {
    bc <- MASS::boxcox(y ~ 1, plotit = FALSE)
    lambda <- bc$x[which.max(bc$y)]
    y <- (y^lambda - 1) / lambda
  } else if (transform != 'none') stop('Accepted transformations are "log", "boxcox" and "none" only. If you need another transformation, please do so first, and then supplied the transformed variable in Qa, with transform = "none".')

  years <- start.year:end.year
  df <- pc[years %in% Qa$year][, y := y]

  # Reconstruction ----------------------------------------------------------
  fit <- lm(y ~ ., data = df, na.action = na.omit)
  rec <- predict(fit, newdata = pc, interval = 'prediction')
  if (transform == 'log') {
    rec <- exp(rec)
  } else {
    if (transform == 'boxcox') rec <- (rec*lambda + 1)^(1/lambda)
  }
  rec <- data.table(rec)
  setnames(rec, c('Q', 'Ql', 'Qu'))
  rec$year <- years

  list(rec = rec,
       coeffs = fit$coefficients,
       sigma = summary(fit)$sigma)
}

#' Cross validation of PCR reconstruction.
#'
#' @inheritParams PCR_reconstruction
#' @param Z A list of cross-validation folds. If `NULL`, will be created with `make_Z()` with default settings. Users are advised to use `make_Z()` to create the cross-validation folds beforehand. See [make_Z] for details.
#' @return A list of cross validation results
#' * metrics.dist: distribution of performance metrics across all cross-validation runs; a matrix, one column for each metric, with column names.
#' * metrics: average performance metrics; a named vector.
#' * obs: the (transformed) observations, a data.table with two columns (year, y)
#' * Ycv: the predicted streamflow in each cross validation run; a matrix, one column for each cross-validation run
#' * Z: the cross-validation folds
#' @export
cvPCR <- function(Qa, pc, start.year, transform = 'log', Z = NULL) {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)
  pc <- as.data.table(pc)
  end.year <- start.year + nrow(pc) - 1
  if (end.year < Qa[.N, year])
    stop('The last year of pc is earlier than the last year of the instrumental period.')

  y <- Qa$Qa
  if (transform == 'log') {
    y <- log(y)
  } else if (transform == 'boxcox') {
    bc <- MASS::boxcox(y ~ 1, plotit = FALSE)
    lambda <- bc$x[which.max(bc$y)]
    y <- (y^lambda - 1) / lambda
  } else if (transform != 'none') stop('Accepted transformations are "log", "boxcox" and "none" only. If you need another transformation, please do so first, and then supplied the transformed variable in Qa, with transform = "none".')

  years <- start.year:end.year
  df <- pc[years %in% Qa$year][, y := y]

  if (is.null(Z)) {
    Z <- make_Z(y)
  } else {
    if (!is.list(Z)) stop("Please provide the cross-validation folds (Z) in a list.")
  }
  # Cross validation ------------------------------------------------------------
  cv <- function(z) {
    # Leave-k-out cross validation with the indices of the k data points to be left out supplied in z
    # Returns a vector of performance metrics, which is calculated on the transformed time series.
    fit <- lm(y ~ ., data = df[-z])
    yhat <- predict(fit, newdata = df)
    list(metric = calculate_metrics(yhat, y, z),
         Ycv = yhat)
  }

  cv_res <- lapply(Z, cv)
  metrics.dist <- sapply(cv_res, '[[', 'metric')
  Ycv <- sapply(cv_res, '[[', 'Ycv')
  dimnames(Ycv) <- NULL
  list(metrics.dist = t(metrics.dist),
       metrics = rowMeans(metrics.dist),
       obs = data.table(year = Qa$year, y = y),
       Ycv = Ycv,
       Z = Z) # Retain Z so that we can plot the CV points when analyzing CV results
}


#' Ensemble PCR reconstruction
#'
#' Same as `PCR` but using ensemble model.  This function builds a reconstruction for each ensemble member then calculates the ensemble average.
#' @inheritParams PCR_reconstruction
#' @param pc.list A list, each element is a set of principal component as in `PCR_reconstruction`'s `pc`
#' @return A list with two elements
#' * rec: the ensemble averaged reconstructed streamflow
#' * ensemble: the ensemble of reconstructions; a matrix, one column for each ensemble member
#' @export
PCR_ensemble <- function(Qa, pc.list, start.year, transform = 'log') {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)
  end.year <- start.year + nrow(pc.list[[1]]) - 1
  if (end.year < Qa[.N, year])
    stop('The last year of the principal components is earlier than the last year of the instrumental period.')

  y <- Qa$Qa
  if (transform == 'log') {
    y <- log(y)
  } else if (transform == 'boxcox') {
    bc <- MASS::boxcox(y ~ 1, plotit = FALSE)
    lambda <- bc$x[which.max(bc$y)]
    y <- (y^lambda - 1) / lambda
  } else if (transform != 'none') stop('Accepted transformations are "log", "boxcox" and "none" only. If you need another transformation, please do so first, and then supplied the transformed variable in Qa, with transform = "none".')

  years <- start.year:end.year
  instPeriod <- which(years %in% Qa$year)

  # Reconstruction ----------------------------------------------------------
  ensemble <- sapply(pc.list,
                     function(pc) {
                       pc <- as.data.table(pc)
                       df <- pc[instPeriod][, y := y]
                       fit <- lm(y ~ ., data = df, na.action = na.omit)
                       predict(fit, newdata = pc)
                     })
  Q <- rowMeans(ensemble)
  if (transform == 'log') {
    Q <- exp(Q)
  } else {
    if (transform == 'boxcox') Q <- (Q*lambda + 1)^(1/lambda)
  }
  rec <- data.table(year = years, Q = Q)
  list(rec = rec, ensemble = ensemble)
}

#' Cross validation of PCR_ensemble reconstruction.
#'
#' @inheritParams cvPCR
#' @inheritParams PCR_ensemble
#' @return A list of cross validation results
#' * metrics.dist: distribution of performance metrics across all cross-validation runs; a matrix, one column for each metric, with column names.
#' * metrics: average performance metrics; a named vector.
#' * obs: the (transformed) observations, a data.table with two columns (year, y)
#' * Ycv: the predicted streamflow in each cross validation run; a matrix, one column for each cross-validation run
#' * Z: the cross-validation folds
#' @export
cvPCR_ensemble <- function(Qa, pc.list, start.year, transform = 'log', Z = NULL) {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)
  end.year <- start.year + nrow(pc.list[[1]]) - 1
  if (end.year < Qa[.N, year])
    stop('The last year of the principal components is earlier than the last year of the instrumental period.')

  y <- Qa$Qa
  if (transform == 'log') {
    y <- log(y)
  } else if (transform == 'boxcox') {
    bc <- MASS::boxcox(y ~ 1, plotit = FALSE)
    lambda <- bc$x[which.max(bc$y)]
    y <- (y^lambda - 1) / lambda
  } else if (transform != 'none') stop('Accepted transformations are "log", "boxcox" and "none" only. If you need another transformation, please do so first, and then supplied the transformed variable in Qa, with transform = "none".')

  years <- start.year:end.year
  instPeriod <- which(years %in% Qa$year)
  df.list <- lapply(pc.list, function(pc) as.data.table(pc)[instPeriod][, y := y][])

  if (is.null(Z)) {
    Z <- make_Z(y)
  } else {
    if (!is.list(Z)) stop("Please provide the cross-validation folds (Z) in a list.")
  }

  # Function to calculate cross validation metrics
  cv <- function(z) {
    ensemble <- sapply(df.list,
                       function(df) {
                         fit <- lm(y ~ ., data = df[-z], na.action = na.omit)
                         predict(fit, newdata = df)
                       })
    yhat <- rowMeans(ensemble)
    list(metric = calculate_metrics(yhat, y, z),
         Ycv = yhat)
  }

  # Cross validation ------------------------------------------------------------
  cv_res <- lapply(Z, cv)
  metrics.dist <- sapply(cv_res, '[[', 'metric')
  Ycv <- sapply(cv_res, '[[', 'Ycv')
  dimnames(Ycv) <- NULL
  list(metrics.dist = t(metrics.dist),
       metrics = rowMeans(metrics.dist),
       obs = data.table(year = Qa$year, y = y),
       Ycv = Ycv,
       Z = Z)
}


#' Select the best reconstruction
#'
#' Select the best reconstruction from an ensemble.
#' @inheritParams cvPCR_ensemble
#' @param agg.type Type of ensemble aggregate. There are 2 options: 'best member': the member with the best performance score is used; 'best overall': if the ensemble average is better than the best member, it will be used, otherwise the best member will be used.
#' @param criterion The performance criterion to be used.
#' @param return.all.metrics Logical, if TRUE, all members' performance scores (and the ensemble average's score, if `agg.type == 'best overall'`) are returned.
#' @return A list of two elements:
#' * choice: The index of the selection. If the ensemble is selected, returns 0.
#' * cv: the cross-validation results of the choice, see [cvPCR] for details.
#' * all.metrics: all members' scores, and if `agg.type == 'best overall'`, the ensemble average's scores as well, in the last column.
#' @export
PCR_ensemble_selection <- function(Qa, pc.list, start.year, transform = 'log', Z = NULL,
                                   agg.type = c('best member', 'best overall'),
                                   criterion = c('RE', 'CE', 'nRMSE', 'KGE'),
                                   return.all.metrics = FALSE) {
  if (!(agg.type %in% c('best member', 'best overall')))
    stop('agg.type must be either "best member" or "best overall".')

  if (!(criterion %in% c('RE', 'CE', 'nRMSE', 'KGE')))
    stop('Criterion must be either RE, CE, nRMSE or KGE')

  memberCV <- lapply(pc.list, function(pc) cvPCR(Qa, pc, start.year, transform, Z))
  memberMetrics <- sapply(memberCV, '[[', 'metrics')
  bestMember <- which.max(memberMetrics[criterion, ])

  if (agg.type == 'best member') {
    ans <- list(choice = bestMember, cv = memberCV[[bestMember]])
    if (return.all.metrics) ans$all.metrics <- memberMetrics
  } else {
    ensembleCV <- cvPCR_ensemble(Qa, pc.list, start.year, transform, Z)
    ensembleMetrics <- ensembleCV$metrics
    ans <- if (ensembleMetrics[criterion] > memberMetrics[criterion, bestMember]) {
      list(choice = 0, cv = ensembleCV)
    } else {
      list(choice = bestMember, cv = memberCV[[bestMember]])
    }
    if (return.all.metrics) ans$all.metrics <- cbind(memberMetrics, ensembleMetrics)
  }
  ans
}
