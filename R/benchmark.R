#' Reconstruction with Principal Component Regression
#'
#' Reconstruct streamflow using the given principal compone11nts and backward selection.
#' Since this is the benchmark for the LDS model, and for simplicity, learning and cross-validation are performed in one go, instead of with two separate functions as in the case of LDS.
#'
#' @inheritParams LDS_reconstruction
#' @param pc A data frame, one colulmn for each principal component
#' @param start.year Starting year of `pc`, i.e, the first year of the paleo period. `start.year + nrow(pc) - 1` will determine the last year of the study horizon, which must be greater than or equal to the last year in `Qa`.
#' @param stepwise If `TRUE`, backward stepwise selection will be performed. Otherwise, all PCs will be used.
#' @return A list of reconstruction results
#' * rec: reconstruction
#' * selected: a vector of selected principal components
#' @export
PCR_reconstruction <- function(Qa, pc, start.year, trans = 'log', stepwise = TRUE) {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)
  pc <- as.data.table(pc)
  N <- nrow(pc)
  end.year <- start.year + N - 1
  if (end.year < Qa[.N, year])
    stop('The last year of pc is earlier than the last year of the instrumental period.')

  if (trans == 'log') {
    y <- log(Qa$Qa)
  } else if (trans == 'boxcox') {
    lambda <- car::powerTransform(Qa$Qa ~ 1)$roundlam %>% 'names<-'('lambda')
    y <- car::bcPower(Qa$Qa, lambda)
  } else if (trans == 'none') {
    y <- Qa$Qa
  } else stop('Accepted transformations are "log", "boxcox" and "none"')

  indStart <- Qa[1, year] - start.year + 1
  indEnd <- N - (end.year - Qa[.N, year])
  df <- cbind(pc[indStart:indEnd], y)

  # Main model  -------------------------------------------------------------------

  # Model fitting
  if (stepwise) {
    fit <- tryCatch(
      step(lm(y ~ . , data = df), direction = 'backward', trace = 0),
      error = function(e) {
        if (substr(e$message, 1, 16) == 'AIC is -infinity') {
          warning('Backward selection finds AIC = -Infinity, only PC1 is used.')
          lm(y ~ PC1 , data = df)
        }
      }
    )
    selected <- names(fit$model)[-1]   # First element is intercept
    if (length(selected) == 0) {
      warning('Backward selection returned empty model; model selection is skipped.')
      fit <- lm(y ~ . , data = df)
      selected <- names(fit$model)[-1]   # First element is intercept
    }
  } else {
    fit <- lm(y ~ . , data = df)
    selected <- names(fit$model)[-1]
  }

  rec <- predict(fit, newdata = pc, interval = 'confidence')
  if (trans == 'log') {
    rec <- exp(rec)
  } else {
    if (trans == 'boxcox') rec <- (rec*lambda + 1)^(1/lambda)
  }
  rec <- data.table(rec)
  setnames(rec, c('Q', 'Ql', 'Qu'))
  rec$year <- start.year:end.year

  list(rec = rec,
       coeffs = fit$coefficients,
       sigma = summary(fit)$sigma,
       selected = selected)
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

cvPCR <- function(Qa, pc, start.year, trans = 'log', k, CV.reps = 100, Z = NULL) {

  # Preprocessing ------------------------------------------------------------------
  Qa <- as.data.table(Qa)
  pc <- as.data.table(pc)
  N <- nrow(pc)
  n <- nrow(Qa)
  end.year <- start.year + N - 1
  if (end.year < Qa[.N, year])
    stop('The last year of pc is earlier than the last year of the instrumental period.')

  if (trans == 'log') {
    y <- log(Qa$Qa)
  } else if (trans == 'boxcox') {
    lambda <- car::powerTransform(Qa$Qa ~ 1)$roundlam %>% 'names<-'('lambda')
    y <- car::bcPower(Qa$Qa, lambda)
  } else if (trans == 'none') {
    y <- Qa$Qa
  } else stop('Accepted transformations are "log", "boxcox" and "none"')

  indStart <- Qa[1, year] - start.year + 1
  indEnd <- N - (end.year - Qa[.N, year])
  df <- cbind(pc[indStart:indEnd], y)

  # Function to calculate cross validation metricsO
  cv <- function(z, df) {
    # Leave-k-out cross validation with the indices of the k data points to be left out supplied in z
    # Returns a vector of performance metrics
    fit <- lm(y ~ . , data = df, subset = setdiff(1:n, z))
    rec <- predict(fit, newdata = df)
    if (trans == 'log') {
      rec <- exp(rec)
    } else {
      if (trans == 'boxcox') rec <- (rec*lambda + 1)^(1/lambda)
    }

    list(metric = calculate_metrics(rec, Qa$Qa, z),
         Ycv = rec)
  }

  # Cross validation ------------------------------------------------------------
  nObs <- Qa[!is.na(Qa), .N]

  if (missing(k)) k <- ceiling(nObs/10) # Default 10%

  if (is.null(Z)) {
    Z <- replicate(CV.reps, sample(1:nObs, size = k, replace = FALSE), simplify = FALSE)
  } else {
    if (!is.list(Z)) stop("Please provide the cross-validation points in a list.")
  }
  cv_res <- if (k > 1) { # Leave-k-out, k > 1
    lapply(Z, cv, df = df)
  } else { # Leave-one-out
    lapply(1:n, cv, df = df)
  }
  metrics.dist <- rbindlist(lapply(cv_res, '[[', 'metric'))
  Ycv <- as.data.table(sapply(cv_res, '[[', 'Ycv'))
  Ycv$year <- Qa$year
  Ycv <- melt(Ycv, id.vars = 'year', variable.name = 'rep', value.name = 'Qa')

  list(metrics.dist = metrics.dist,
       metrics = colMeans(metrics.dist),
       Ycv = Ycv,
       Z = Z)
}


#' Same as `PCR` but using ensemble model
#'
#' @inheritParams PCR_reconstruction
#' @param pc.list A list, each element is a set of principal component as in `PCR_reconstruction`'s `pc`
#' @export
PCR_ensemble <- function(Qa, pc.list, start.year, trans = 'log', stepwise = TRUE) {

  # Non-standard call issue in R CMD check
  Q <- NULL

  ensemble <- lapply(pc.list, function(pc)
    PCR_reconstruction(Qa, pc, start.year, trans, stepwise))

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
cvPCR_ensemble <- function(Qa, pc.list, start.year, trans = 'log',
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
    recResults <- PCR_ensemble(Qa2, pc.list, start.year, trans, stepwise = FALSE)
    rec <- recResults$rec[year %in% instYears, Qa]

    list(metric = calculate_metrics(rec, Qa$Qa, z),
         Ycv = rec)
  }

  # Cross validation ------------------------------------------------------------
  nObs <- Qa[!is.na(Qa), .N]

  if (missing(k)) k <- ceiling(nObs/10) # Default 10%

  if (is.null(Z)) {
    Z <- replicate(CV.reps, sample(1:nObs, size = k, replace = FALSE), simplify = FALSE)
  } else {
    if (!is.list(Z)) stop("Please provide the cross-validation points in a list.")
  }
  cv_res <- if (k > 1) { # Leave-k-out, k > 1
    lapply(Z, cv, Qa = Qa, pc.list = pc.list)
  } else { # Leave-one-out
    lapply(1:n, cv, Qa = Qa, pc.list = pc.list)
  }
  metrics.dist <- rbindlist(lapply(cv_res, '[[', 'metric'))
  Ycv <- as.data.table(sapply(cv_res, '[[', 'Ycv'))
  Ycv$year <- Qa$year
  Ycv <- melt(Ycv, id.vars = 'year', variable.name = 'rep', value.name = 'Qa')

  list(metrics.dist = metrics.dist,
       metrics = colMeans(metrics.dist),
       Ycv = Ycv,
       Z = Z)
}
