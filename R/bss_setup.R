#' script for reading/writing setup configuration file (bssr.ini) for bss
#' Usually this should be executed by the user immediately after installing bssr
#'
#' @export
setup <- function(brainsuite_path = NULL, quiet = FALSE, raise_error = TRUE) {
  bssr_ini_file <- get_bssr_ini_path()
  bs_settings <- ini::read.ini(bssr_ini_file)

  if (is.null(brainsuite_path)) { # The user didn't specify the BrainSuite location
    message('Finding BrainSuite installation paths...', appendLF = FALSE)
    # Find BrainSuite installation paths automatically
    brainsuite_path <- get_brainsuite_install_path()
  }

  # Check if BrainSuite atlas files are present in the user specified location
  if (check_bs_atlas_exists(brainsuite_path, quiet = quiet, raise_error = raise_error)) {
    # At this point, a valid brainsuite_path should exist
    # Write it to the bssr.ini file
    bs_settings$path$brainsuite_path <- brainsuite_path
    message(bssr_ini_file, appendLF = TRUE)
    ini::write.ini(bs_settings, bssr_ini_file)
    message('bssr setup is complete.', appendLF = TRUE)
  }
  else
    message(paste('Warning: bssr setup is not complete.\n',
                  'After making sure BrainSuite is installed, please run bssr::setup("/path/to/brainsuite/") manually.', sep = ""), appendLF = TRUE)
}

get_os <- function() {
  if ( !is.null(Sys.info()) ) {
    if(Sys.info()['sysname'] == 'Darwin')
      return('macOS')
  }
  else if ( grepl("darwin", R.version$os) )
    return('macOS')

  if (.Platform$OS.type == 'unix')
    return('unix')

  if (.Platform$OS.type == 'windows')
    return('windows')
}

get_brainsuite_path_on_macOS <- function() {
  bs_paths <- sort(list.files('/Applications', 'BrainSuite', full.names = TRUE), decreasing = TRUE)[1]
  if (!is.na(bs_paths[1])) {
    # Check if the required files exist
    if (check_bs_atlas_exists(bs_paths[1], quiet = TRUE, raise_error = FALSE)) return(bs_paths[1]) else return("")
  }
  else
    return("")
}

get_brainsuite_path_on_unix <- function() {

  # Search /opt first
  bs_paths <- sort(list.files('/opt', 'BrainSuite', full.names = TRUE), decreasing = TRUE)[1]
  if (is.na(bs_paths[1])) {
    # Search home directory first
    bs_paths <- sort(list.files(path.expand('~'), 'BrainSuite', full.names = TRUE), decreasing = TRUE)[1]
    if (!is.na(bs_paths[1])) {
      # Check if the required files exist

    }
    else
      return(NULL)

  }
}

get_brainsuite_path_on_windows <- function() {

  message('get_brainsuite_path_on_windows not implemented.', appendLF = TRUE)

}

check_bs_atlas_exists <- function(brainsuite_path, quiet=FALSE, raise_error = TRUE) {

  if (!quiet) message('Finding BrainSuite atlas file paths...', appendLF = FALSE)
  for (i in bs_atlas_files ) {
    errmesg <- sprintf('Atlas file %s does not exist. \nPlease check if BrainSuite is installed correctly.', i)
    if (!check_file_exists(file.path(brainsuite_path, i), raise_error = raise_error,
                          errmesg = errmesg)) {
      if (!quiet) message(errmesg, appendLF = TRUE)
      return(FALSE)
    }
  }
  if (!quiet) message('Done.', appendLF = TRUE)
  return(TRUE)
}


#' Check if the BrainSuite installation is valid.
#' This function is called from .onLoad() when the package is loaded.
#' Opens bssr.ini and checks if all the paths are valid.
#'
#' @export
is_brainsute_installed <- function() {
  bssr_ini_file <- get_bssr_ini_path()
  if (check_file_exists(bssr_ini_file, raise_error = FALSE)) {
    bs_settings <- ini::read.ini(bssr_ini_file)
    if (check_bs_atlas_exists(bs_settings$path$brainsuite_path, quiet = TRUE, raise_error = FALSE))
      return(TRUE)
    else
      return(FALSE)
  }
  else
    return(FALSE)
}

#' Retrieve bssr.ini path in the package
#'
#' @export
get_bssr_ini_path <- function() {
  bssr_ini_file <- system.file("extdata", "bssr.ini", package = 'bssr')
  if (check_file_exists(bssr_ini_file, raise_error = TRUE)) return(bssr_ini_file) else return("")
}


#' Retrieve BrainSuite installation path
#'
#' @export
get_brainsuite_install_path <- function() {
  switch(get_os(),
         macOS = {brainsuite_path <- get_brainsuite_path_on_macOS()},
         unix = {brainsuite_path <- get_brainsuite_path_on_unix()},
         windows = {brainsuite_path <- get_brainsuite_path_on_windows()}
  )
  return(brainsuite_path)
}

get_brainsuite_path_from_bssr_ini <- function() {
  bssr_ini_file <- get_bssr_ini_path()
  bs_settings <- ini::read.ini(bssr_ini_file)
  return(bs_settings$path$brainsuite_path)
}

#' Check if a given BrainSuite path is valid
#'
#' @export
is_valid_brainsute_install_path <- function(brainsuite_path="") {
  return(check_bs_atlas_exists(brainsuite_path, quiet = TRUE, raise_error = FALSE))
}
