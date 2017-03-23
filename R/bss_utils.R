#' Utility functions for bss

#' @export
check_file_exists <- function(filename, raise_error=FALSE, errmesg=NULL) {
  errmesg <- if (is.null(errmesg)) sprintf('File %s does not exist.', filename) else errmesg
  if ( file.exists(filename) )
    return(TRUE)
  else {
    if (raise_error)
      stop(errmesg, call. = FALSE)
    else
      return(FALSE)
  }
}
