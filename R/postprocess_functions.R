#' Calculate some metrics from reconstructed and observed streamflow
#'
#' @param Qa.hat reconstructed
#' @param Qa observed
#' @param z indices of left-out values in cross validation
#' @export
calculate_metrics <- function(Qa.hat, Qa, z) {

    data.table(
        R2    = 1 - var(Qa.hat[-z] - Qa[-z], na.rm = TRUE) / var(Qa[-z], na.rm = TRUE),
        RE    = 1 - ssq(Qa.hat[z], Qa[z]) / ssq(rep(mean(Qa[-z], na.rm = TRUE), length(z)), Qa[z]),
        CE    = NSE(Qa.hat[z], Qa[z]),
        nRMSE = rmse(Qa.hat[z], Qa[z]) / mean(Qa, na.rm = TRUE),
        KGE   = KGE(Qa.hat[z], Qa[z])
    )
}
