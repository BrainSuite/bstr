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

#' Load ROI data
#'
#' Loads the ROI measure for statistical analysis.
#' @param subjects_dir character string for subject directory.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' @param roiid numeric label identifier for the region of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. Should either be "gmthickness", "gmvolume", or "wmvolume".
#' @param outdir output directory name to save results for ROI analysis.
#' @export
#'
bss_load_roi_data <- function(subjects_dir, csv, roiids, roimeas, outdir=NULL) {

  if (!dir.exists(subjects_dir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", subjects_dir), call. = FALSE)
  }
  if (is.null(outdir)) {
    # Outputted directory (if not specified) is of the same format (csv or tsv) as the inputted file
    if (substr(csv,nchar(csv)-2,nchar(csv)) == "csv"){
      outdir <- paste0(substr(csv,1,nchar(csv)-4), "_roidata.csv")
      specific_separator = ','
    } else if (substr(csv,nchar(csv)-2,nchar(csv)) == "tsv"){
      outdir <- paste0(substr(csv,1,nchar(csv)-4), "_roidata.tsv")
      specific_separator = '\t'
    }
    cat(sprintf('Output directory is not specified. Using %s to save outputs.\n', outdir))
  }
  else {
    dir.create(file.path(outdir), showWarnings = FALSE)
  }
  demographics <- read_demographics(csv)
  if("File_roi" %in% colnames(demographics)){
    warning(sprintf("The file %s already contains a File_roi column.\nWill overwrite this column.\n", csv), call.=FALSE)
  }

  demographics$subjID <- as.character(demographics$subjID)
  # for subjID in demographic_data['subjID'].astype('str')]
  roiwise_file_list <- lapply(demographics$subjID,
                              function(i) {
                                Sys.glob(file.path(subjects_dir, i, "*roiwise.stats.txt"))[1]
                              })
  if (any(is.na(roiwise_file_list))) {
    stop(sprintf("Subject %s is missing the roiwise.stats.txt file.\n",
                 demographics$subjID[which(is.na(roiwise_file_list))]), call. = FALSE)
  }

  demographics$File_roi <- unlist(roiwise_file_list)

  roi_data_frame <- read_roi_data_for_all_subjects(demographics$File_roi, roiids, roimeas)

  # Put outputted data frame together with demographics data frame
  combined_roidata_and_demographics <- cbind(demographics[,-which(names(demographics) == "File_roi")], roi_data_frame)
  write.table(combined_roidata_and_demographics, outdir, sep = specific_separator,row.names = FALSE)

  bss_data <- list(df=combined_roidata_and_demographics, outdir=outdir, roiids=roiids, roimeas=roimeas, csv=outdir)

  # Finally also include the command to load the data
  pasted_roiids <- paste(roiids,collapse = ", ")
  bss_data$load_data_command <- sprintf("bss_data <- bss_load_roi_data( '%s', '%s', c( %s), '%s', outdir = '%s') ",
                                        subjects_dir, csv, pasted_roiids, roimeas, outdir)

  return(bss_data)
}



