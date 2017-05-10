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

#' Read ROI stats file
#'
#' Reads the BrainSuite ROI stats file saved for each subject after executing SVREG.
#' @param roiwise_txt_filename filename for the ROIwise stats in the individual subject directory
#' @param roiid numeric label identifier for the region of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. Should either be "gmthickness",
#' "gmvolume", or "wmvolume".
#' @export
read_roistats_txt <- function(roiwise_txt_filename, roiid, roimeas = 'gmthickness') {

  measure_dict <- list(gmthickness = "Mean_Thickness(mm)",
                       gmvolume = "GM_Volume(mm^3)",
                       area = "Cortical_Area_pial(mm^2)",
                       swmFA = "swmFA",
                       swmMD = "swmMD",
                       swmRD = "swmRD",
                       swmAD = "swmAD"
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
  # rownames(roiwise_stats) <- roiwise_stats[,1]
  if ( !(roiid %in% roiwise_stats$ROI_ID) ) {
    stop(sprintf('ROI ID %d not found in the roiwise stats file %s.\nPlease check if %d is a valid ROI.\n',
                 roiid, roiwise_txt_filename, roiid), call. = FALSE)
  }
  return (roiwise_stats[roiwise_stats$ROI_ID == roiid, as.character(measure_dict[roimeas])])
}

read_roi_data_for_all_subjects <- function(roi_filelist, roiid, roimeas = 'gmthickness') {

  roi_data <- lapply(roi_filelist, function(i) {
    read_roistats_txt(i, roiid, roimeas)
  })
  return(unlist(roi_data))
}


