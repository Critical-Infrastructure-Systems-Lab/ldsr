#' ldsr: Paleo reconstruction with linear dynamical systems
#'
#' Reconstruct paleoclimate variables (streamflow, precipitation, etc...) with linear dynamical systems
#' @docType package
#' @name ldsr
#' @importFrom magrittr %>%
#' @importFrom data.table rbindlist
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom Rcpp evalCpp
#' @useDynLib ldsr, .registration = TRUE
NULL

.onUnload <- function(libpath) {
    library.dynam.unload("ldsr", libpath)
}
