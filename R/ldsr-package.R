#' @docType package
#' @name ldsr
#' @import data.table
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom dplR tbrm
#' @useDynLib ldsr, .registration = TRUE
#' @keywords internal
"_PACKAGE"
NULL

.onUnload <- function(libpath) {
    library.dynam.unload("ldsr", libpath)
}
