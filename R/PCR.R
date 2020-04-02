#' Principal Component Regression Reconstruction
#'
#' Reconstruction with principal component linear regression.
#' @param Qa Observations: a data.frame of annual streamflow with at least two columns: year and Qa.
#' @param pc For a single model: a data.frame, one column for each principal component. For an ensemble reconstruction: a list, each element is a data.frame of principal components.
#' @param start.year Starting year of the climate proxies, i.e, the first year of the paleo period. `start.year + nrow(pc) - 1` will determine the last year of the study horizon, which must be greater than or equal to the last year in `Qa`.
#' @param transform Flow transformation, either "log", "boxcox" or "none".
#' @return A list of reconstruction results, with the following elements:
#' ## For a single-model reconstruction:
#' * rec: reconstructed streamflow with 95% prediction interval; a data.table with four columns: year, Q, Ql (lower bound), and Qu (upper bound).
#' * coeffs: the regression coefficients.
#' * sigma: the residual standard deviation.
#' ## For an ensemble reconstruction:
#' * rec: the ensemble average reconstruction; a data.table with two columns: year and Q.
#' * ensemble: a list of ensemble members, each element is reconstructed from one element of `pc` and is itself a list of three elements: Q (a vector of reconstructed flow), coeffs and sigma.
#' Note that for ensemble reconstruction, `ldsr` does not provide uncertainty estimates. It is up to the user to do so, for example, using ensemble spread.
#' @export
PCR_reconstruction <- function(Qa, pc, start.year, transform = 'log') {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)
  single <- is.data.frame(pc)
  N <- if (single) nrow(pc) else nrow(pc[[1]])
  end.year <- start.year + N - 1
  if (end.year < Qa[.N, year])
    stop('The last year of the principal components is earlier than the last year of the instrumental period.')

  y <- Qa$Qa
  if (transform == 'log') {
    y <- log(y)
  } else if (transform == 'boxcox') {
    bc <- MASS::boxcox(y ~ 1, plotit = FALSE)
    lambda <- bc$x[which.max(bc$y)]
    y <- (y^lambda - 1) / lambda
  } else if (transform != 'none') stop('Accepted transformations are "log", "boxcox" and "none" only. If you need another transformation, please do so first, and then supplied the transformed variable in Qa, and set transform = "none".')

  years <- start.year:end.year
  instPeriod <- which(years %in% Qa$year)

  # Reconstruction ----------------------------------------------------------
  if (single) {
    pc <- as.data.table(pc)
    df <- pc[years %in% Qa$year][, y := y]
    fit <- stats::lm(y ~ ., data = df, na.action = na.omit)
    rec <- stats::predict(fit, newdata = pc, interval = 'prediction')
    if (transform == 'log') {
      rec <- exp(rec)
    } else {
      if (transform == 'boxcox') rec <- (rec*lambda + 1)^(1/lambda)
    }
    rec <- data.table(rec)
    setnames(rec, c('Q', 'Ql', 'Qu'))
    rec$year <- years
    ans <- list(rec = rec,
                coeffs = fit$coefficients,
                sigma = stats::sigma(fit))
  } else {
    ensemble <- lapply(pc,
                       function(pcX) {
                         setDT(pcX)
                         df <- pcX[instPeriod][, y := y]
                         fit <- lm(y ~ ., data = df, na.action = na.omit)
                         list(Q = stats::predict(fit, newdata = pcX),
                              coeffs = fit$coefficients,
                              sigma = stats::sigma(fit))
                       })
    Q <- rowMeans(sapply(ensemble, '[[', 'Q'))
    if (transform == 'log') {
      Q <- exp(Q)
    } else {
      if (transform == 'boxcox') Q <- (Q*lambda + 1)^(1/lambda)
    }
    rec <- data.table(year = years, Q = Q)
    ans <- list(rec = rec, ensemble = ensemble)
  }
  ans
}

#' One cross-validation run
#'
#' Make one prediction for one cross-validation run. This is a subroutine to be called by other cross-validation functions.
#' @param df The training data, with one column named y, the (transformed) observations. and other columns the predictors.
#' @param z A vector of left-out points
#' @return A vector of prediction.
#' @export
one_pcr_cv <- function(df, z) {
  fit <- lm(y ~ ., data = df[-z])
  predict(fit, newdata = df)
}

#' Cross validation of PCR reconstruction.
#'
#' @inheritParams PCR_reconstruction
#' @param Z A list of cross-validation folds. If `NULL`, will be created with `make_Z()` with default settings. Users are advised to use `make_Z()` to create the cross-validation folds beforehand. See [make_Z] for details.
#' @param metric.space Either "transformed" or "original", the space to calculate the performance metrics.
#' @return A list of cross validation results
#' * metrics.dist: distribution of performance metrics across all cross-validation runs; a matrix, one column for each metric, with column names.
#' * metrics: average performance metrics; a named vector.
#' * obs: the (transformed) observations, a data.table with two columns (year, y)
#' * Ycv: the predicted streamflow in each cross validation run; a matrix, one column for each cross-validation run
#' * Z: the cross-validation folds
#' @export
cvPCR <- function(Qa, pc, start.year, transform = 'log', Z = NULL, metric.space = 'transformed') {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)
  if (is.data.frame(pc)) {
    N <- nrow(pc)
    pc.list <- list(pc)
  } else {
    N <- nrow(pc[[1]])
    pc.list <- pc
  }
  end.year <- start.year + N - 1
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
  instPeriod <- which(years %in% Qa$year)
  df.list <- lapply(pc.list, function(pc) as.data.table(pc)[instPeriod][, y := y][])

  if (is.null(Z)) {
    Z <- make_Z(y)
  } else {
    if (!is.list(Z)) stop("Please provide the cross-validation folds (Z) in a list.")
  }

  # Cross validation ------------------------------------------------------------
  Ycv <- lapply(Z, function(z) rowMeans(sapply(df.list, one_pcr_cv, z = z)))

  if (metric.space == 'original') {
    if (transform == 'log') {
      Ycv <- lapply(Ycv, exp)
    } else if (transform == 'boxcox') {
      Ycv <- lapply(Ycv, function(x) (x*lambda + 1)^(1/lambda))
    }
    target <- Qa$Qa
  } else {
    target <- y
  }
  # doing mapply is a lot faster than working on data.table
  metrics.dist <- mapply(calculate_metrics, sim = Ycv, z = Z, MoreArgs = list(obs = target))
  metrics.dist <- data.table(t(metrics.dist))
  metrics <- metrics.dist[, lapply(.SD, mean)]

  names(Ycv) <- seq_along(Ycv)
  setDT(Ycv, check.names = FALSE)
  Ycv[, year := Qa$year]

  list(metrics.dist = metrics.dist,
       metrics = metrics,
       target = data.table(year = Qa$year, y = target),
       Ycv = melt(Ycv, id.vars = 'year', variable.name = 'rep', value.name = 'Y'),
       Z = Z) # Retain Z so that we can plot the CV points when analyzing CV results
}

#' Select the best reconstruction
#'
#' Select the best reconstruction from an ensemble.
#' @inheritParams cvPCR
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
    ensembleCV <- cvPCR(Qa, pc.list, start.year, transform, Z)
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
