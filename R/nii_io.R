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

read_nii_images_for_all_subjects_old <- function(nii_filelist, attrib_siz, mask_idx = NULL) {
  data_matrix <- vapply(nii_filelist, function(i) {
    cat(sprintf('%s\n', i))
    as.vector(RNifti::readNifti(i))
  }, FUN.VALUE = numeric(attrib_siz))
  colnames(data_matrix) <- NULL
  if ( !is.null(mask_idx) )
    data_matrix <- data_matrix[mask_idx,]
  return(t(data_matrix))
}

read_nii_images_for_all_subjects <- function(nii_filelist, attrib_siz, mask_idx = NULL) {

  if (is.null(mask_idx))
    mask_idx <- 1:attrib_siz

  data_matrix <- matrix(NA_real_, nrow = length(nii_filelist), ncol = length(mask_idx))
  for (ii in seq_along(nii_filelist)) {
    x <- RNifti::readNifti(nii_filelist[ii])
    data_matrix[ii, ] <- x[mask_idx]
    cat(sprintf('Loaded %s\n', nii_filelist[ii]))
  }

  colnames(data_matrix) <- NULL
  gc()
  return(data_matrix)
}

