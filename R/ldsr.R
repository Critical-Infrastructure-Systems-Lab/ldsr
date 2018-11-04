#' ldsr: Paleo reconstruction with linear dynamical systems
#'
#' Reconstruct paleoclimate variables (streamflow, precipitation, etc...) with linear dynamical systems
#' @docType package
#' @name ldsr
#' @importFrom magrittr %>%
#' @importFrom data.table rbindlist
NULL

.onUnload <- function(libpath) {
    library.dynam.unload("ldsr", libpath)

}
