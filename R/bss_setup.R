# BrainSuite Statistics Toolbox in R (bssr)
# Copyright (C) 2017 The Regents of the University of California
# Creator: Shantanu H. Joshi, Department of Neurology, Ahmanson Lovelace Brain Mapping Center, UCLA
#
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation; version 2.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License version 2 for more details.
#
# You should have received a copy of the GNU General Public License along with this program;
# if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

#' Setup script for bssr
#'
#' This script reads from and writes to the setup configuration file (bssr.ini) for bss.
#' Usually this will be called automatically when the package is installed and loaded
#' for the first time. Optionally, it can be executed by the user immediately after installing bssr.
#'
#' @param brainsuite_path path to the BrainSuite installation
#' @param quiet logical; if \code{FALSE} does not display messages to the user
#' @param raise_error logical; if \code{TRUE}, stops the execution if file does not exist. The default
#' value is \code{FALSE}, in which case the function returns {FALSE} without stopping the execution.
#' @export
setup <- function(brainsuite_path = NULL, quiet = FALSE, raise_error = TRUE) {
  bssr_ini_file <- get_bssr_ini_path()
  bs_settings <- ini::read.ini(bssr_ini_file)

  if (is.null(brainsuite_path)) { # The user didn't specify the BrainSuite location
    message('Finding BrainSuite installation paths...', appendLF = FALSE)
    # Find BrainSuite installation paths automatically
    brainsuite_path <- get_brainsuite_install_path(quiet, raise_error)
  }

  # Check if BrainSuite atlas files and binaries are present in the user specified location
  if (check_bs_atlas_binaries_exist(brainsuite_path, quiet = quiet, raise_error = raise_error)) {
    # At this point, a valid brainsuite_path should exist
    # Write it to the bssr.ini file
    bs_settings$path$brainsuite_path <- brainsuite_path
    message(bssr_ini_file, appendLF = TRUE)
    ini::write.ini(bs_settings, bssr_ini_file)
    message('bssr setup is complete.', appendLF = TRUE)
  }

  else
    message(paste('bssr setup is not complete.\n',
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

get_brainsuite_path_on_macOS <- function(quiet = TRUE, raise_error = FALSE) {

  # Search /Applications first
  bs_opt_paths <- dir('/Applications', 'BrainSuite', full.names = TRUE)
  bs_opt_paths <- bs_opt_paths[dir.exists(bs_opt_paths)]
  # Search ~ next
  bs_home_paths <- dir('~', 'BrainSuite', full.names = TRUE)
  bs_home_paths <- bs_home_paths[dir.exists(bs_home_paths)]
  bs_paths <- sort(c(bs_opt_paths, bs_home_paths), decreasing = TRUE)

  # Test bs_paths for valid installations
  valid_bs_path = ""
  for (path in bs_paths) {
    if (check_bs_atlas_binaries_exist(path, quiet = FALSE, raise_error = FALSE)) {
      valid_bs_path = path
      break
    }
  }
  return(valid_bs_path)
}

get_brainsuite_path_on_unix <- function(quiet = TRUE, raise_error = FALSE) {

  # Search /opt first
  bs_opt_paths <- dir('/opt', 'BrainSuite', full.names = TRUE)
  bs_opt_paths <- bs_opt_paths[dir.exists(bs_opt_paths)]
  # Search ~ next
  bs_home_paths <- dir('~', 'BrainSuite', full.names = TRUE)
  bs_home_paths <- bs_home_paths[dir.exists(bs_home_paths)]
  bs_paths <- sort(c(bs_opt_paths, bs_home_paths), decreasing = TRUE)

  # Test bs_paths for valid installations
  valid_bs_path = ""
  for (path in bs_paths) {
    if (check_bs_atlas_binaries_exist(path, quiet = FALSE, raise_error = FALSE)) {
      valid_bs_path = path
      break
    }
  }
  return(valid_bs_path)
}

get_brainsuite_path_on_windows <- function(quiet = TRUE, raise_error = FALSE) {

  # Search C:/Program Files
  bs_paths <- dir('C:/Program Files', 'BrainSuite', full.names = TRUE)
  bs_paths <- bs_paths[dir.exists(bs_paths)]
  bs_paths <- sort(bs_paths, decreasing = TRUE)

  # Test bs_paths for valid installations
  valid_bs_path = ""
  for (path in bs_paths) {
    if (check_bs_atlas_binaries_exist(path, quiet = FALSE, raise_error = FALSE)) {
      valid_bs_path = path
      break
    }
  }
  return(valid_bs_path)
}

check_bs_atlas_binaries_exist <- function(brainsuite_path, quiet=FALSE, raise_error = TRUE) {

  if (!quiet) message('Finding BrainSuite atlas file paths...', appendLF = FALSE)
  for (i in bs_atlas_files ) {
    errmesg <- sprintf('Atlas file %s does not exist. \nPlease check if BrainSuite is installed correctly.', file.path(brainsuite_path, i))
    if (!check_file_exists(file.path(brainsuite_path, i), raise_error = raise_error,
                          errmesg = errmesg)) {
      if (!quiet) message(errmesg, appendLF = TRUE)
      return(FALSE)
    }
  }
  for (i in bs_binary_files ) {
    errmesg <- sprintf('Binary file %s does not exist. \nPlease check if BrainSuite is installed correctly.', file.path(brainsuite_path, i))
    if (!check_file_exists(file.path(brainsuite_path, i), raise_error = raise_error,
                           errmesg = errmesg)) {
      if (!quiet) message(errmesg, appendLF = TRUE)
      return(FALSE)
    }
  }
  if (!quiet) message(sprintf('Done. Valid BrainSuite installation found at %s', brainsuite_path), appendLF = TRUE)
  return(TRUE)
}


#' Check if BrainSuite is installed.
#'
#' Check if the BrainSuite installation is valid by verifying if the appropriate
#' atlas files and data exist. This function is called from \code{\link{.onLoad}}
#' when the package is loaded. It opens \code{bssr.ini} and checks if all the
#' paths are valid.
#' @param  quiet boolean specifying whether warnings/messages should be displayed
#' @param  raise_error boolean specifying whether an exception should be raised
#'
#' @export
is_brainsute_installed <- function(quiet = FALSE, raise_error = FALSE) {
  bssr_ini_file <- get_bssr_ini_path()
  if (check_file_exists(bssr_ini_file, raise_error = raise_error)) {
    bs_settings <- ini::read.ini(bssr_ini_file)
    if (check_bs_atlas_binaries_exist(bs_settings$path$brainsuite_path, quiet = quiet, raise_error = raise_error))
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

#' Retrieve label description path in the package
#'
#' @export
get_labeldesc_path <- function() {
  labeldesc_file <- file.path(get_brainsuite_install_path(), "labeldesc", "brainsuite_labeldescriptions_30March2018.xml")
  if (check_file_exists(labeldesc_file, raise_error = TRUE)) return(labeldesc_file) else return("")
}

#' Retrieve BrainSuite installation path
#' @param  quiet boolean specifying whether warnings/messages should be displayed
#' @param  raise_error boolean specifying whether an exception should be raised
#'
#' @export
get_brainsuite_install_path <- function(quiet = TRUE, raise_error = FALSE) {
  brainsuite_path_bssr_ini <- get_brainsuite_path_from_bssr_ini()
  if (is_valid_brainsute_install_path(brainsuite_path_bssr_ini))
    return(brainsuite_path_bssr_ini)
  switch(get_os(),
         macOS = {brainsuite_path <- get_brainsuite_path_on_macOS(quiet, raise_error)},
         unix = {brainsuite_path <- get_brainsuite_path_on_unix(quiet, raise_error)},
         windows = {brainsuite_path <- get_brainsuite_path_on_windows(quiet, raise_error)}
  )
  return(brainsuite_path)
}

#' Retrieve BrainSuite installation path from bssr.ini
#'
#' @export
get_brainsuite_path_from_bssr_ini <- function() {
  bssr_ini_file <- get_bssr_ini_path()
  bs_settings <- ini::read.ini(bssr_ini_file)
  return(bs_settings$path$brainsuite_path)
}

#' Check if a given BrainSuite path is valid
#'
#' @param brainsuite_path path to the BrainSuite installation
#' @export
is_valid_brainsute_install_path <- function(brainsuite_path="") {
  return(check_bs_atlas_binaries_exist(brainsuite_path, quiet = TRUE, raise_error = FALSE))
}
