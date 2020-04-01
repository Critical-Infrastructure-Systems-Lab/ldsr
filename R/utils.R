calendar_to_water_year <- function(df, lag) {

    # Lag: the difference between month 1 of the calendar year and month 1 of the water year
    # eg: if the water year starts in April then lag = 3
    # Lag is negative if the water year begins at the end of the calendar year
    # eg: if the water year starts in December then lag = -1

    if (lag > 0) {
        df$year <- as.integer(ifelse(df$month %in% 1:lag, df$year-1, df$year))
        df$month <- as.integer(ifelse(df$month %in% 1:lag, df$month+12-lag, df$month-lag))
        # Remove incomplete years at the beginning and end
        df <- df[(lag+1) : (.N - 12 + lag)]
    } else if (lag < 0) {
        lag <- -lag
        df$month <- as.integer(ifelse((df$month+lag) %in% (lag+1):12, df$month+lag, df$month+lag-12))
        df$year <- as.integer(ifelse((df$month+lag) %in% (lag+1):12, df$year, df$year+1))
    }
    return(df)
}

#' Water to calendar year
#'
#' Converting a sub-annual data frame from water year to calendar year. Users must pay attention to the notion of the water year, whether it uses the year in which it starts (e.g., water year 2000 is from April 2000 to March 2001) or the year in which it ends (e.g., water year 2000 is from September 1999 to October 2000).
#' @param df A data frame following water year with at least two columns
#' \describe{
#'     \item{year}{The water year}
#'     \item{months}{The months in the water year, with `month = 1` indicating the first month of the water year}
#' }
#' @param dt The time shift required to change from water year to calendar year. For example, if the water year 2000 starts from April 2000 then `delta.t = 3` (caldendar year is 3 months a head of water year). On the other hand, if the water year 2000 starts from October 1999 then `delta.t = -3` (calendear year is 3 months behind water year).
#' @param keep.all If `FALSE` (the default), only complete years are retained. Otherwise, the whole data frame is retained.
#' @return A new data frame where the year and month is now shifted to calendar year, with `month = 1` indicating January.
#' @details For example, in Thailand a water year is from April to March, i.e, month 1 of 1921 denotes April 1921, and month 12 of 1921 denotes March 1922. After shifting, the period of April to December 1921 are discarded, and the data frame starts from January 1922.
#' @export
water_to_calendar_year <- function(df, dt, keep.all = FALSE) {

    df <- as.data.table(df)
    dt <- as.integer(dt)
    if (dt > 0) {
        k <- which(df$month + dt > 12)
        df[k, ':='(month = month + dt - 12L,
                   year = year + 1L)]
        df[-k, month := month + dt]
        if (keep.all) df[] else df[(13-dt):(.N - dt)]
    } else {
        dt <- -dt
        k <- which(df$month - dt < 1)
        df[k, ':='(month = month - dt + 12L,
                   year = year - 1L)]
        df[-k, month := month - dt]
        # delta.t has changed sign
        if (keep.all) df[] else df[(dt + 1):(.N - (12 - dt))]
    }
}

#' Transform the estimates before calculating metrics
#'
#' If you already ran the cross-validation on transformed output and now wanted to calculate performance on the back-transformed one, or vice-versa, you don't have to rerun the whole cross-validation, but just need to transform or back-transform the cross-validation Ycv. This function helps you do that.
#' @param cv Cross-validation output as produced by cvLDS or cvPCR
#' @param transform Either "log", "exp", "boxcox" or "inv_boxcox"
#' @param lambda Lambda value used in Box-Cox or inverse Box-Cox
#' @return A new cv object wit hthe new metrics
#' @export
metrics_with_transform <- function(cv, transform, lambda = NULL) {
  Ycv <- copy(cv$Ycv)
  Ycv[,
      Y := switch(transform,
                 log = log(Y),
                 exp = exp(Y),
                 boxcox = {if (lambda == 0) log(Y) else (Y^lambda - 1) / lambda},
                 inv_boxcox = {if (lambda == 0) exp(Y) else (Y*lambda + 1)^(1/lambda)})
      ]
  obs <- copy(cv$obs)
  obs[,
      y := switch(transform,
                 log = log(y),
                 exp = exp(y),
                 boxcox = {if (lambda == 0) log(y) else (y^lambda - 1) / lambda},
                 inv_boxcox = {if (lambda == 0) exp(Y) else (y*lambda + 1)^(1/lambda)})
      ]

  metrics.dist <- Ycv[, data.table(t(calculate_metrics(Y, obs$y, cv$Z[[rep]]))), by = rep]
  metrics.dist[, rep := NULL]
  metrics <- metrics.dist[, lapply(.SD, mean)]

  list(metrics.dist = metrics.dist,
       metrics = metrics,
       target = obs,
       Ycv = Ycv,
       Z = Z)
}

#' Reconstruction metrics
#'
#' Calculate reconstruction metrics from the instrumental period
#' @param sim A vector of reconstruction output for instrumental period
#' @param obs A vector of all observations
#' @param z A vector of left out indices in cross validation
#' @param norm.fun The function (unquoted name) used to calculate the normalizing constant. Default is `mean()`, but other functions such as `sd()` can also be used. THe function must take a vector as input and return a scalar as output, and must have an argument `na.rm = TRUE`.
#' @return A named vector of performance metrics
#' @export
calculate_metrics <- function(sim, obs, z, norm.fun = mean) {
    train.sim <- sim[-z]
    train.sim <- train.sim[!is.na(train.sim)]
    train.obs <- obs[-z]
    train.obs <- train.obs[!is.na(train.obs)]
    val.sim <- sim[z]
    val.obs <- obs[z]

    c(R2    = NSE(train.sim, train.obs), # Use the NSE form of R2
      RE    = RE(val.sim, val.obs, mean(train.obs)),
      CE    = NSE(val.sim, val.obs),
      nRMSE = nRMSE(val.sim, val.obs, norm.fun(obs, na.rm = TRUE)),
      KGE   = KGE(val.sim, val.obs)
    )
}

#' Make cross-validation folds.
#'
#' Make a list of cross-validation folds. Each element of the list is a vector of the cross-validation points for one cross-validation run.
#' @param obs Vector of observations.
#' @param nRuns Number of repetitions.
#' @param frac Fraction of left-out points. For leave-one-out, use `frac = 1`, otherwise use any value less than 1. Default is 0.1 (leave-10%-out).
#' @param contiguous Logical. If `TRUE`, the default, the left-out points are made in contiguous blocks; otherwise, they are scattered randomly.
#' @export
make_Z <- function(obs, nRuns = 30, frac = 0.1, contiguous = TRUE) {
  obsInd <- which(!is.na(obs))
  if (frac == 1) {
    split(obsInd, obsInd)
  } else {
    n <- length(obsInd)
    k <- floor(n * frac) # leave-k-out
    if (contiguous) {
      maxInd <- n - k # Highest possible position in of obsInd
      if (maxInd < nRuns) { # Not enough samples, reduce k
        maxInd <- nRuns
        k <- n - nRuns
      }
      lapply(sort(sample(1:maxInd, nRuns)), function(x) obsInd[x:(x + k)])
    } else {
      replicate(nRuns, sort(sample(obsInd, k, replace = FALSE)), simplify = FALSE)
    }
  }
}
