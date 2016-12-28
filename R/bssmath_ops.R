# math ops for BSS
#' @export
log10_transform <- function(values) {
  eps <- .Machine$double.eps
  sgn <- (values + eps)/abs(values + eps)
  logvalues <- - 1*sgn*log10(abs(values + eps))
  return(logvalues)
}

#' @export
image_to_shape <- function(imagefile, shapefile, outputshapefile, resample) {

}
