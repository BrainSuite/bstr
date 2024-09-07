# BrainSuite Statistics Toolbox in R (bstr)
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

#' Check if file exists
#'
#'
#' @param filename Name of file.
#' @param raise_error logical; if `TRUE`, stops the execution if file does not exist. The default
#' value is `FALSE`, in which case the function returns {FALSE} without stopping the execution.
#' @param errmesg character string of optional error message
#'
#' @export
check_file_exists <- function(filename, raise_error=FALSE, errmesg=NULL) {
  errmesg <- if (is.null(errmesg)) sprintf('File or path %s does not exist.', filename) else errmesg
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
    if (raise_error)
      stop(errmesg, call. = FALSE)
    else
      return(FALSE)
  }
}

#' Check if multiple files exists all at once
#'
#' @param filelist Name of file.
#' @param errmesg character string of optional error message
#'
#' @export
check_multiple_files_exists <- function(filelist, errmesg=NULL) {
  if ( !all(file.exists(filelist)) ) {
    message('Following subjects have missing files')
    print(filelist[which(!file.exists(filelist))], row.names = FALSE)
    stop('\nCheck if svreg was run succesfully on all the subjects.',
         call. = FALSE)
  }
  return (TRUE)
}

#' Delete and recreate directory
#'
#' @param dir Name of directory.
#'
#' @export
delete_and_recreate_dir <- function(dir){
  unlink(dir,recursive=TRUE)
  dir.create(dir)
}

#' Convert a 3D coordinate from an image to a linear vector index
#' Equivalent to sub2ind from MATLAB(R)
#' Does not perform error checking for bounds
#' @param coord vector of 3D coordinates
#' @param dimI vector of image dimensions
#' @export
coord2ind <- function(coord, dimI) {
  return (coord[1] + (coord[2]-1)*dimI[1] + (coord[3]-1)*dimI[1]*dimI[2])
}

#' Convert a linear array index to a 3D coordinate in an image
#' Equivalent to ind2sub from MATLAB(R)
#' Does not perform error checking for bounds
#' @param idx numeric index in a linear vector
#' @param dimI vector of image dimensions
#' @export
ind2coord <- function(idx, dimI) {
  return (as.vector(arrayInd(idx, dimI)))
}
