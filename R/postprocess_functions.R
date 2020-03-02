#' Calculate some metrics from reconstructed and observed streamflow
#'
#' @param Qa.hat reconstructed
#' @param Qa observed
#' @param z indices of left-out values in cross validation
#' @export
calculate_metrics <- function(Qa.hat, Qa, z) {

    data.table(
        R2    = 1 - var(Qa.hat[-z] - Qa[-z], na.rm = TRUE) / var(Qa[-z], na.rm = TRUE),
        RE    = 1 - hydroGOF::ssq(Qa.hat[z], Qa[z]) / hydroGOF::ssq(rep(mean(Qa[-z], na.rm = TRUE), length(z)), Qa[z]),
        CE    = hydroGOF::NSE(Qa.hat[z], Qa[z]),
        nRMSE = hydroGOF::rmse(Qa.hat[z], Qa[z]) / mean(Qa, na.rm = TRUE),
        KGE   = hydroGOF::KGE(Qa.hat[z], Qa[z])
    )
}
