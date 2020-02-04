#' Reconstruction with Principal Component Regression
#'
#' Reconstruct streamflow using the given principal components and backward selection.
#' Since this is the benchmark for the LDS model, and for simplicity, learning and cross-validation are performed in one go, instead of with two separate functions as in the case of LDS.
#'
#' @inheritParams cvLDS
#' @param pc A data frame of principal components with first column called years (the years in) the
#' study period, followed by one column for each principal components.
#' @param stepwise If `TRUE`, backward stepwise selection will be performed. Otherwise, all PCs will be used.
#' @return A list of reconstruction and cross-validation results
#' * rec: reconstruction
#' * metrics: performance metrics
#' * metrics.dist: distribution of performance metrics across all cross-validation runs.
#' * Ycv: the predicted streamflow in each cross validation run
#' * Z: the cross-validation points
#' @export
PCR_reconstruction <- function(Qa, pc, stepwise = TRUE, trans = 'log', k, CV.reps = 100, Z = NULL) {

  # Function to calculate cross validation metrics
  cv <- function(z) {
    # Leave-k-out cross validation with the indices of the k data points to be left out supplied in z
    # Returns a vector of performance metrics
    fit <- lm(log(Qa) ~ . , data = df2, subset = setdiff(1:N, z))
    Qa.hat <- exp(predict(fit, newdata = df2))

    return(list(metric = calculate_metrics(Qa.hat, df2$Qa, z), Ycv = Qa.hat))
  }

  # Main model  -------------------------------------------------------------------
  Qa <- as.data.table(Qa)[year <= max(pc$year)]
  pc <- as.data.table(pc)

  # Important to remove year otherwise it'll be come a predictor
  years <- pc$year
  df <- merge(pc, Qa, by = 'year')
  df$year <- NULL
  pc$year <- NULL
  if (trans == 'log') {
    df[, Qa := log(Qa)]
  } else if (trans == 'boxcox') {
    lambda <- car::powerTransform(Qa$Qa ~ 1)$roundlam %>% 'names<-'('lambda')
    df[, Qa := car::bcPower(Qa, lambda)]
  } else { # Has to be 'none' here
    if (trans != 'none') stop('Accepted transformations are "log", "boxcox" and "none".')
  }

  # Model fitting
  if (stepwise) {
    fit <- tryCatch(
      step(lm(Qa ~ . , data = df), direction = 'backward', trace = 0),
      error = function(e) {
        if (substr(e$message, 1, 16) == 'AIC is -infinity') {
          warning('Backward selection finds AIC = -Infinity, only PC1 is used.')
          lm(log(Qa) ~ PC1 , data = df)
        }
      }
    )

    selected <- names(fit$model)[-1]   # First element is intercept
    if (length(selected) == 0) {
      warning('Backward selection returned empty model; model selection is skipped.')
      fit <- lm(Qa ~ . , data = df)
      selected <- names(fit$model)[-1]   # First element is intercept
    }
  } else {
    fit <- lm(Qa ~ . , data = df)
    selected <- names(fit$model)[-1]
  }

  rec <- data.table(exp(predict(fit, newdata = pc, interval = 'confidence')))
  colnames(rec) <- c('Q', 'Ql', 'Qu')
  rec$year <- years

  df2 <- df[, selected, with = FALSE]
  df2$Qa <- df$Qa

  # Cross validation ------------------------------------------------------------
  N <- Qa[!is.na(Qa), .N]
  if (missing(k)) k <- ceiling(N/10) # Default 10%
  if (is.null(Z)) {
    Z <- replicate(CV.reps, sample(1:N, size = k, replace = FALSE))
  } else {
    if (is.list(Z)) {
      CV.reps <- length(Z)
      Z <- matrix(unlist(Z), ncol = CV.reps)
    } else {
      CV.reps <- ncol(Z)
    }
  }
  cv_res <- if (k > 1) { # Leave-k-out, k > 1
    apply(Z, 2, cv)
  } else { # Leave-one-out
    lapply(1:N, cv)
  }
  metrics.dist <- rbindlist(lapply(cv_res, '[[', 'metric'))
  Ycv <- as.data.table(sapply(cv_res, '[[', 'Ycv'))
  Ycv$year <- Qa$year
  Ycv <- melt(Ycv, id.vars = 'year', variable.name = 'rep', value.name = 'Qa')

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

#' Same as `PCR` but using ensemble model
#'
#' @inheritParams PCR_reconstruction
#' @param pc.list A list, each element is a set of principal component as in `PCR_reconstruction`'s `pc`
#' @export
PCR_ensemble <- function(Qa, pc.list, stepwise = TRUE, trans = 'log', k, CV.reps = 100, Z = NULL) {

  # Non-standard call issue in R CMD check
  Q <- NULL

  # Function to calculate cross validation metrics
  oneCV <- function(z, df) {
    # Leave-k-out cross validation with the indices of the k data points to be left out supplied in z
    # Returns a vector of performance metrics

    fit <- lm(Qa ~ . , data = df, subset = setdiff(1:N, z))
    Qa.hat <- exp(predict(fit, newdata = df))
    return(list(metric = calculate_metrics(Qa.hat, df$Qa, z), Ycv = Qa.hat))
  }

  # Main model  -------------------------------------------------------------------
  Qa <- as.data.table(Qa)[year <= max(pc.list[[1]]$year)]
  pc.list <- lapply(pc.list, as.data.table)

  # Important to remove year otherwise it'll be come a predictor
  ensemble <- lapply(pc.list, function(pc) {
    years <- pc$year
    df <- merge(pc, Qa, by = 'year')
    df$year <- NULL
    pc$year <- NULL
    if (trans == 'log') {
      df[, Qa := log(Qa)]
    } else if (trans == 'boxcox') {
      lambda <- car::powerTransform(Qa$Qa ~ 1)$roundlam %>% 'names<-'('lambda')
      df[, Qa := car::bcPower(Qa, lambda)]
    } else {
      if (trans != 'none') stop('Accepted transformations are "log", "boxcox" and "none".')
    }

    # Model fitting
    if (stepwise) {
      fit <- tryCatch(
        step(lm(Qa ~ . , data = df), direction = 'backward', trace = 0),
        error = function(e) {
          if (substr(e$message, 1, 16) == 'AIC is -infinity') {
            warning('Backward selection finds AIC = -Infinity, only PC1 is used.')
            lm(log(Qa) ~ PC1 , data = df)
          }
        }
      )

      selected <- names(fit$model)[-1]   # First element is intercept
      if (length(selected) == 0) {
        warning('Backward selection returned empty model; model selection is skipped.')
        fit <- lm(Qa ~ . , data = df)
        selected <- names(fit$model)[-1]   # First element is intercept
      }
    } else {
      fit <- lm(Qa ~ . , data = df)
      selected <- names(fit$model)[-1]
    }
    rec <- data.table(exp(predict(fit, newdata = pc, interval = 'confidence')))
    colnames(rec) <- c('Q', 'Ql', 'Qu')
    rec$year <- years

    outCols <- c(selected, 'Qa')
    list(rec = rec,
         coeffs = fit$coefficients,
         sigma = summary(fit)$sigma,
         selectedPC = df[, outCols, with = FALSE])

  })

  rec <- rbindlist(lapply(ensemble, '[[', 'rec'))[, list(Qa = mean(Q)), by = year]

  # Cross validation ------------------------------------------------------------
  N <- Qa[!is.na(Qa), .N]
  if (missing(k)) k <- ceiling(N/10) # Default 10%
  if (is.null(Z)) {
    Z <- replicate(CV.reps, sample(1:N, size = k, replace = FALSE))
  } else {
    if (is.list(Z)) { # LDS takes Z as a list but for legacy this function takes Z as a matrix
      CV.reps <- length(Z)
      Z <- matrix(unlist(Z), ncol = CV.reps)
    } else {
      CV.reps <- ncol(Z)
    }
  }
  metrics.dist <- sapply(1:length(pc.list), function(i) {

    df <- ensemble[[i]]$selectedPC
    cv_res <- if (k > 1) { # Leave-k-out, k > 1
      apply(Z, 2, oneCV, df = df)
    } else { # Leave-one-out
      lapply(1:N, oneCV, df = df)
    }
    lapply(cv_res, '[[', 'metric') %>% rbindlist() %>% colMeans()
  }) %>%
    t() %>%
    as.data.table()

  return(list(
    rec = rec,
    ensemble = ensemble,
    metrics.dist = metrics.dist,
    metrics = colMeans(metrics.dist),
    Z = t(Z)
  ))
}
