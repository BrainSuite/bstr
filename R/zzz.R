#' @importFrom grDevices col2rgb colorRamp colorRampPalette dev.new
#' @importFrom graphics axis plot rect
#' @importFrom methods .valueClassTest new
#' @importFrom stats formula model.matrix p.adjust pf pt
#' @importFrom utils read.csv read.table write.csv
NULL

.onLoad <- function(libname, pkgname) {

  if(!is_brainsute_installed()) {
    setup(quiet = FALSE, raise_error = FALSE)
    invisible()
  }
}
