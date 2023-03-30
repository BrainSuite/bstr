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

#' Load ROI data from subjects (BIDS compatible)
#'
#' Loads the ROI measure for statistical analysis.
#' @param subjects_dir character string for subject directory.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' @param roiids vector of numeric label identifiers for the region of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. Should either be "gmthickness", "gmvolume", or "wmvolume".
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @export
#'
bstr_load_roi_data <- function(subjects_dir, csv, roiids, roimeas, exclude_col) {

  if (!dir.exists(subjects_dir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", subjects_dir), call. = FALSE)
  }

  demographics <- read_demographics(csv, exclude_col)
  if("File_roi" %in% colnames(demographics)){
    warning(sprintf("The file %s already contains a File_roi column.\nWill overwrite this column.\n", csv), call.=FALSE)
  }

  demographics$subjID <- as.character(demographics$subjID)
  if (roimeas %in% roi_types$svreg_roi_types || roimeas %in% roi_types$swm_roi_types) {
    roi_bids_filelist <- file.path(subjects_dir, demographics$subjID, 'anat',
                                   sprintf('%s%s%s', demographics$subjID, '_T1w', bs_file_formats$roi_txt))
  }
  else if (roimeas %in% roi_types$wm_roi_types) {
    roi_bids_filelist <- file.path(subjects_dir, demographics$subjID, 'dwi',
                                   sprintf('%s%s%s', demographics$subjID, '_T1w', bs_file_formats$wm_roi_txt))
  }

  if(all(file.exists(roi_bids_filelist)))
    additional_dir <- 'anat'
  else
    additional_dir <- '/'
  if (roimeas %in% roi_types$svreg_roi_types || roimeas %in% roi_types$swm_roi_types) {
    roiwise_file_list <- lapply(demographics$subjID,
                              function(i) {
                                Sys.glob(file.path(subjects_dir, i, additional_dir, "*roiwise.stats.txt"))[1]
                              })
  }
  else if (roimeas %in% roi_types$wm_roi_types) {
    roiwise_file_list <- lapply(demographics$subjID,
                              function(i) {
                                Sys.glob(file.path(subjects_dir, i, additional_dir, "*wm.stats.tsv"))[1]
                              })
  }


  if (any(is.na(roiwise_file_list))) {
    stop(sprintf("Subject %s is either missing the roiwise.stats.txt or the wm.stats.tsv file depending on the ROI measure used.\n",
                 demographics$subjID[which(is.na(roiwise_file_list))]), call. = FALSE)
  }

  demographics$File_roi <- unlist(roiwise_file_list)

  roi_data_frame <- read_roi_data_for_all_subjects(demographics$File_roi, demographics,roiids, roimeas)

  # Put outputted data frame together with demographics data frame
  combined_roidata_and_demographics <- cbind(demographics[,-which(names(demographics) == "File_roi")], roi_data_frame)

  bstr_data <- list(df=combined_roidata_and_demographics, roiids=roiids, roimeas=roimeas)

  # Finally also include the command to load the data
  pasted_roiids <- paste(roiids,collapse = ", ")
  bstr_data$load_data_command <- sprintf("bstr_data <- bstr_load_roi_data( '%s', '%s', c( %s), '%s', '%s') ",
                                        subjects_dir, csv, pasted_roiids, roimeas, exclude_col)

  return(bstr_data)
}



