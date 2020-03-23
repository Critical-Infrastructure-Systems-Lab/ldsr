#' ldsr: Paleo reconstruction with linear dynamical systems
#'
#' Reconstruct paleoclimate variables (streamflow, precipitation, etc...) with linear dynamical systems
#' @docType package
#' @name ldsr
#' @import data.table
#' @import stats
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib ldsr, .registration = TRUE
#' @keywords internal
"_PACKAGE"
NULL

.onUnload <- function(libpath) {
    library.dynam.unload("ldsr", libpath)
}
