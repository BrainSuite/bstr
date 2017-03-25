.onLoad <- function(libname, pkgname) {

  if(!is_brainsute_installed()) {
    packageStartupMessage(paste('Running setup...', sep = ""), appendLF = FALSE)
    setup(quiet = FALSE, raise_error = FALSE)
  }
}
