#' Reconstruction with Principal Component Regression
#'
#' Reconstruct streamflow using the given principal compone11nts and backward selection.
#' Since this is the benchmark for the LDS model, and for simplicity, learning and cross-validation are performed in one go, instead of with two separate functions as in the case of LDS.
#'
#' @inheritParams LDS_reconstruction
#' @param pc A data frame, one colulmn for each principal component
#' @param start.year Starting year of `pc`, i.e, the first year of the paleo period. `start.year + nrow(pc) - 1` will determine the last year of the study horizon, which must be greater than or equal to the last year in `Qa`.
#' @param transform Flow transformation, either "log", "boxcox" or "none".
#' @param stepwise If `TRUE`, backward stepwise selection will be performed. Otherwise, all PCs will be used.
#' @return A list of reconstruction results
#' * rec: reconstruction
#' * selected: a vector of selected principal components
#' @export
PCR_reconstruction <- function(Qa, pc, start.year, transform = 'log') {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)
  pc <- as.data.table(pc)
  N <- nrow(pc)
  end.year <- start.year + N - 1
  if (end.year < Qa[.N, year])
    stop('Error in PCR_reconstruction(): the last year of pc is earlier than the last year of the instrumental period.')

  if (transform == 'log') {
    y <- log(Qa$Qa)
  } else if (transform == 'boxcox') {
    lambda <- car::powerTransform(Qa$Qa ~ 1)$roundlam %>% 'names<-'('lambda')
    y <- car::bcPower(Qa$Qa, lambda)
  } else if (transform == 'none') {
    y <- Qa$Qa
  } else stop('Error in PCR_reconstruction(): accepted transformations are "log", "boxcox" and "none" only.')

  # Training  -------------------------------------------------------------------

  years <- start.year:end.year
  df <- cbind(pc[years %in% Qa$year], y)
  fit <- lm(y ~ . , data = df, na.action = na.omit)

  # Reconstruction ----------------------------------------------------------

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
#' @inheritParams cvLDS
#' @param pc Selected principal components. No input selection is done during cross-validation.
#' @return A list of cross validation results
#' * metrics: performance metrics
#' * metrics.dist: distribution of performance metrics across all cross-validation runs.
#' * Ycv: the predicted streamflow in each cross validation run
#' * Z: the cross-validation points
#' @export

cvPCR <- function(Qa, pc, start.year, transform = 'log', k, CV.reps = 100, Z = NULL) {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)
  pc <- as.data.table(pc)
  N <- nrow(pc)
  end.year <- start.year + N - 1
  if (end.year < Qa[.N, year])
    stop('Error in PCR_reconstruction(): the last year of pc is earlier than the last year of the instrumental period.')

  if (transform == 'log') {
    y <- log(Qa$Qa)
  } else if (transform == 'boxcox') {
    lambda <- car::powerTransform(Qa$Qa ~ 1)$roundlam
    names(lambda) <-'lambda'
    y <- car::bcPower(Qa$Qa, lambda)
  } else if (transform == 'none') {
    y <- Qa$Qa
  } else stop('Error in PCR_reconstruction(): accepted transformations are "log", "boxcox" and "none" only.')

  years <- start.year:end.year
  df <- cbind(pc[years %in% Qa$year], y)

  # Leave-k-out cross validation with the indices of the k data points to be left out supplied in z
  # Returns a vector of performance metrics
  cv <- function(z, df) {
    fit <- lm(y ~ . , data = df[-z])
    rec <- predict(fit, newdata = df)
    if (transform == 'log') {
      rec <- exp(rec)
    } else {
      if (transform == 'boxcox') rec <- (rec*lambda + 1)^(1/lambda)
    }
    list(metric = calculate_metrics(rec, Qa$Qa, z),
         Ycv = rec)
  }

  # Cross validation ------------------------------------------------------------
  obsInd <- which(!is.na(Qa$Qa))
  nObs <- length(obsInd)

  if (is.null(Z)) {
    if (missing(k)) k <- ceiling(nObs/10) # Default 10%
    Z <- if (k == 1) {
      split(obsInd, obsInd) # Leave-one-out
    } else {
      replicate(CV.reps, sample(obsInd, size = k, replace = FALSE), simplify = FALSE)
    }
  } else {
    if (!is.list(Z)) stop("Error in cvPCR: please provide the cross-validation points (Z) in a list.")
  }
  CV.reps <- length(Z)
  cv_res <- lapply(Z, cv, df = df)
  metrics.dist <- sapply(cv_res, '[[', 'metric')
  Ycv <- cbind(Qa$year, sapply(cv_res, '[[', 'Ycv'))
  colnames(Ycv) <- c('year', paste0('cv', 1:CV.reps))
  list(metrics.dist = as.data.table(t(metrics.dist)),
       metrics = rowMeans(metrics.dist),
       Ycv = Ycv,
       Z = Z) # Retain Z so that we can plot the CV points when analyzing CV results
}


#' Same as `PCR` but using ensemble model
#'
#' @inheritParams PCR_reconstruction
#' @param pc.list A list, each element is a set of principal component as in `PCR_reconstruction`'s `pc`
#' @param agg.type Type of ensemble aggregate. There are 3 options: 'average': the ensemble average is returned; 'best member': the member with the best performance score is used; 'best overall': if the ensemble average is better than the best member, it will be used, otherwise the best member will be used.
#' @param criterion The performance criterion to be used.
#' @export
PCR_ensemble <- function(Qa, pc.list, start.year, transform = 'log',
                         agg.type = c('average', 'best member', 'best overall'),
                         criterion = c('RE', 'CE', 'nRMSE', 'KGE')) {

  # Non-standard call issue in R CMD check
  # Q <- NULL

   ensemble <- lapply(pc.list, PCR_reconstruction,
                     Qa = Qa, start.year = start.year, transform = transform)

  rec <- rbindlist(lapply(ensemble, '[[', 'rec'))[, list(Qa = mean(Q)), by = year]


  list(rec = rec,
       ensemble = ensemble)
}

#' Cross validation of PCR_ensemble reconstruction.
#'
#' @inheritParams cvPCR
#' @inheritParams PCR_ensemble
#' @return A list of cross validation results
#' * metrics: performance metrics
#' * metrics.dist: distribution of performance metrics across all cross-validation runs.
#' * Ycv: the predicted streamflow in each cross validation run
#' * Z: the cros
#' @export
cvPCR_ensemble <- function(Qa, pc.list, start.year, transform = 'log',
                           k, CV.reps = 100, Z = NULL) {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)

  # Function to calculate cross validation metricsO
  cv <- function(z, Qa, pc.list) {
    # Leave-k-out cross validation with the indices of the k data points to be left out supplied in z
    # Returns a vector of performance metrics
    instYears <- Qa$year
    Qa2 <- copy(Qa)
    Qa2[z, Qa := NA]
    recResults <- PCR_ensemble(Qa2, pc.list, start.year, transform)
    rec <- recResults$rec[year %in% instYears, Qa]

    list(metric = calculate_metrics(rec, Qa$Qa, z),
         Ycv = rec)
  }

  # Cross validation ------------------------------------------------------------
  obsInd <- which(!is.na(Qa$Qa))
  nObs <- length(obsInd)

  if (is.null(Z)) {
    if (missing(k)) k <- ceiling(nObs/10) # Default 10%
    Z <- if (k == 1) {
      split(obsInd, obsInd) # Leave-one-out
    } else {
      replicate(CV.reps, sample(obsInd, size = k, replace = FALSE), simplify = FALSE)
    }
  } else {
    if (!is.list(Z)) stop("Error in cvPCR_ensemble: please provide the cross-validation points (Z) in a list.")
  }
  cv_res <- lapply(Z, cv, Qa = Qa, pc.list = pc.list)
  metrics.dist <- rbindlist(lapply(cv_res, '[[', 'metric'))
  Ycv <- as.data.table(sapply(cv_res, '[[', 'Ycv'))
  Ycv$year <- Qa$year
  Ycv <- melt(Ycv, id.vars = 'year', variable.name = 'rep', value.name = 'Qa')

  list(metrics.dist = metrics.dist,
       metrics = colMeans(metrics.dist),
       Ycv = Ycv)
}
