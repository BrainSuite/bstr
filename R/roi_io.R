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

#' Read ROI stats file
#'
#' Reads the BrainSuite ROI stats file saved for each subject after executing SVREG.
#' @param roiwise_txt_filename filename for the ROIwise stats in the individual subject directory
#' @param roiids vector of numeric label identifiers for the region of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. Should either be "gmthickness",
#' "gmvolume", or "wmvolume".
#' @export
read_roistats_txt <- function(roiwise_txt_filename, roiids, roimeas = 'gmthickness') {

  measure_dict <- list(gmthickness = "Mean_Thickness(mm)",
                       gmvolume = "GM_Volume(mm^3)",
                       area = "Cortical_Area_pial(mm^2)",
                       wmvolume = "WM_Volume(mm^3)",
                       swmFA = "swmFA",
                       swmMD = "swmMD",
                       swmRD = "swmRD",
                       swmAD = "swmAD",
                       axial_WM = "axial_WM",
                       FA_WM = "FA_WM",
                       FRT_GFA_WM = "FRT_GFA_WM",
                       L2_WM = "L2_WM",
                       L3_WM = "L3_WM",
                       mADC_WM = "mADC_WM",
                       MD_WM = "MD_WM",
                       radial_WM = "radial_WM"
  )

  if (!file.exists(roiwise_txt_filename)) {
    stop(sprintf("File name %s does not exist.\n", roiwise_txt_filename), call. = FALSE)
  }
  # Check if file is a roiwise.stats.txt file.
  file_con <- file(roiwise_txt_filename,open="r")
  on.exit(close(file_con))
  file_contents <- readLines(file_con)
  if (!substr(file_contents[1],1,6) == "ROI_ID") {
    stop(sprintf('The file %s is not a valid roiwise.stats.txt file.\n', roiwise_txt_filename), call. = FALSE)
  }
  roiwise_stats <- read.table(roiwise_txt_filename, header = TRUE, check.names = FALSE)

  # Check if all the inputted roiids are in the file
  roiids_not_in_file <- !(roiids %in% roiwise_stats$ROI_ID)
  if (any(roiids_not_in_file == TRUE)) {
    incorrect_roiids <- roiids[roiids_not_in_file]
    stop(sprintf('ROI ID %d not found in the roiwise stats file %s.\nPlease check if %d is a valid ROI.\n',
                 incorrect_roiids, roiwise_txt_filename, incorrect_roiids), call. = FALSE)
  }

  roiids_from_file <- sapply(roiids, function(x) {
    roiwise_stats[roiwise_stats$ROI_ID == x,as.character(measure_dict[roimeas])]
    })

  return(roiids_from_file)

}

#' Reads the BrainSuite ROI data file saved for each subject and returns the filenames
#' @param roi_filelist filenames for the ROIwise stats in the individual subject directory
#' @param demographics data frame containing the demographics of all subjects
#' @param roiids vector of numeric label identifiers for the region of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. Should either be "gmthickness",
#' "gmvolume", or "wmvolume".
#' @export
read_roi_data_for_all_subjects <- function(roi_filelist, demographics,roiids, roimeas = 'gmthickness') {

  roi_data_frame <- data.frame(matrix(nrow=length(roi_filelist),ncol=length(roiids)))
  rownames(roi_data_frame) <- demographics$subjID
  for (i in 1:length(roiids)){
    if (!(roiids[i] %in% read_label_desc()$roiid)){
      stop(sprintf('ROI ID %d not found. This is not a valid ROI for %s.', roiids[i], roimeas))
    }
    colnames(roi_data_frame)[i] <- paste0(as.character(get_roi_tag(read_label_desc(),roiids[i])), "(",roiids[i],")")
  }


  for (x in 1:length(roiids)) {
    roi_data_frame[,x] <- sapply(roi_filelist, function(i,roiid_column) {
      read_roistats_txt(i, roiids[x], roimeas)
    },roiid_column=roiids[x])
  }

  # Add file path column to data frame
  roi_data_frame$ROI_filelist <- roi_filelist

  return(roi_data_frame)
}


