#' ldsr: Paleo reconstruction with linear dynamical systems
#'
#' Reconstruct paleoclimate variables (streamflow, precipitation, etc...) with linear dynamical systems
#' @docType package
#' @name ldsr
#' @import data.table
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom hydroGOF ssq NSE KGE rmse
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom stats runif var
#' @useDynLib ldsr, .registration = TRUE
#' @keywords internal
"_PACKAGE"
NULL

.onUnload <- function(libpath) {
    library.dynam.unload("ldsr", libpath)
}
