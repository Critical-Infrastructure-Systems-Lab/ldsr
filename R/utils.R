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
