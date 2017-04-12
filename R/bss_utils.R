#' Check if file exists
#'
#'
#' @param filename Name of file.
#' @param raise_error logical; if \code{TRUE}, stops the execution if file does not exist. The default
#' value is \code{FALSE}, in which case the function returns {FALSE} without stopping the execution.
#' @param errmesg character string of optional error message
#'
#' @export
check_file_exists <- function(filename, raise_error=FALSE, errmesg=NULL) {
  errmesg <- if (is.null(errmesg)) sprintf('File %s does not exist.', filename) else errmesg
  if (!identical(filename, character(0))) {
    if ( file.exists(filename) )
      return(TRUE)
    else {
      if (raise_error)
        stop(errmesg, call. = FALSE)
      return(FALSE)
    }
  }
  else {
    return(FALSE)
  }
}

