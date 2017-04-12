#' log10 transform
#' @param values numeric vector
#' @return numeric vector containing the log10 transformed \code{values}
#' @export
log10_transform <- function(values) {
  eps <- .Machine$double.eps
  sgn <- (values + eps)/abs(values + eps)
  logvalues <- - 1*sgn*log10(abs(values + eps))
  return(logvalues)
}

# TBD
image_to_shape <- function(imagefile, shapefile, outputshapefile, resample) {

}
