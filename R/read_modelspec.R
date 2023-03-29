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

#' Read modelspec.ini file
#'
#' The modelspec file specifies the subject directory, paths to the atlas files,
#' the model specification, whether linear regerssion, correlation etc.
#' The modelspec file was used in the previous versions of bstr.
#' Currently, the user does not have to create this explicitly.
#' This functionality is still kept in case a future need arises to generate these automatically.
#' @param modelspecfile path to the modelspec file
#' @export
read_modelspec <- function(modelspecfile) {

  mspec <- check_modelspec_validity(modelspecfile)

  # Open the svreg.log file and get the atlas file name
  if ( identical(tools::file_ext(mspec$csv), 'csv') ) {
    demo <- read_demographics(mspec$csv)
  }
  # The first column has to contain subject IDs which are same as subject directories
  first_subjid <- demo[[1]][1]
  # Get atlas names from log files.
  svreg_log_file <- file.path(mspec$subjdir, first_subjid, sprintf('%s.svreg.log', first_subjid))
  if ( !file.exists(svreg_log_file) ) {
    stop(sprintf('Subject %s does not contain the svreg log file %s.\nPlease check if this is a valid subject directory and if svreg was run on all subjects.', first_subjid, svreg_log_file), call. = FALSE)
  }
  bs_atlas_path <- get_brainsuite_atlas_path_from_logfile(svreg_log_file)

  # Get atlas files
  switch(mspec$stats$type,
         sba = { mspec <- get_sba_atlas_files(mspec, bs_atlas_path, svreg_log_file) },
         tbm = { mspec <- get_tbm_atlas_files(mspec, bs_atlas_path, svreg_log_file) },
         croi = { mspec <- get_roi_specs(mspec) }
         )

  return(mspec)
}

#' Check that modelspec.ini file is valid
#'
#' The modelspec file specifies the subject directory, paths to the atlas files,
#' the model specification, whether linear regerssion, correlation etc.
#' The modelspec file was used in the previous versions of bstr.
#' Currently, the user does not have to create this explicitly.
#' This functionality is still kept in case a future need arises to generate these automatically.
#' @param modelspecfile path to the modelspec file
#'
check_modelspec_validity <- function(modelspecfile) {

  check_file_exists(modelspecfile, raise_error = TRUE)

  mspec <- ini::read.ini(modelspecfile)
  mspec$mspec_file <- modelspecfile

  if (is.null(mspec$subject))
    stop(sprintf('The modelspec file %s does not contain a [subject] section.', modelspecfile), call. = FALSE)

  if (is.null(mspec$stats))
    stop(sprintf('The modelspec file %s does not contain a [stats] section.', modelspecfile), call. = FALSE)

  # Check for existence of [subject] sub-fields
  if (is.null(mspec$subject$subjdir))
    stop('Section subjdir under [subject] not found. It should point to the top level directory that contains individual subjects.', call. = FALSE)
  else
    mspec$subjdir <- mspec$subject$subjdir

  if (is.null(mspec$subject$demographics))
    stop('Section demographics under [subject] not found. It should point to the csv/xls demographics file.', call. = FALSE)
  else
    mspec$csv <- mspec$subject$demographics

  if (!is.null(mspec$subject$smooth)) {
    if ( is.na(suppressWarnings(as.numeric(mspec$subject$smooth))) )
      stop('Smoothing level is not numeric.', call. = FALSE)
    mspec$smooth <- as.numeric(mspec$subject$smooth)
    message(sprintf('Using smoothing level %2.1f.', mspec$smooth))
  }
  else
    message('Using no smoothing.')

  valid_stats_fields <- c("type", "main_effect", "covariates", "corr_var")
  stats_fields <- names(mspec$stats)
  if ( !all(stats_fields %in% valid_stats_fields) )
    stop(sprintf('[stats] contains invalid subsection(s). Valid subsections are %s.',
                 paste(valid_stats_fields, collapse = ', ')), call. = FALSE)

  # Check for existence of [stats] fields
  if (is.null(mspec$stats$type))
    stop('Section type under [stats] not found. It should be either sba, tbm or croi.', call. = FALSE)
  else {
    # Check if type is sba, tbm or roi
    if (! (identical(mspec$stats$type, 'tbm') || identical(mspec$stats$type, 'sba') || identical(mspec$stats$type, 'croi')))
      stop('Section type under [stats] should be either sba, tbm or croi.', call. = FALSE)
    mspec$type <- mspec$stats$type
  }

  # Check if main_effect, covariates, and/or corr_var are present
  main_effect_present <- !is.null(mspec$stats$main_effect)
  covariates_present <- !is.null(mspec$stats$covariates)
  corr_var_present <- !is.null(mspec$stats$corr_var)

  # Either main_effect and covariates or corr_var must be present
  condition1 <- main_effect_present & covariates_present & !corr_var_present
  condition2 <- !main_effect_present & !covariates_present & corr_var_present
  if ( !(condition1 || condition2) )
    stop('Either both main_effect and covariates *or* only corr_var must be present.
         All three cannot be specified at this time. If you include the main_effect, you also need to include covariates.', call. = FALSE)

  mspec$main_effect <- mspec$stats$main_effect
  mspec$covariates <- mspec$stats$covariates
  mspec$corr_var <- mspec$stats$corr_var

  return (mspec)
}

#' Get the cortical surface modelspec.ini file
#'
#' The modelspec file specifies the subject directory, paths to the atlas files,
#' the model specification, whether linear regerssion, correlation etc.
#' The modelspec file was used in the previous versions of bstr.
#' Currently, the user does not have to create this explicitly.
#' This functionality is still kept in case a future need arises to generate these automatically.
#' @param mspec path to the modelspec file
#' @param bs_atlas_path path to the atlas file
#' @param svreg_log_file file containing svreg log output
#'
get_sba_atlas_files <- function(mspec, bs_atlas_path, svreg_log_file) {

  lh_atlas_file <- file.path(dirname(bs_atlas_path), bs_file_formats$surf_atlas_left)
  rh_atlas_file <- file.path(dirname(bs_atlas_path), bs_file_formats$surf_atlas_right)

  if ( file.exists(lh_atlas_file) &&  file.exists(rh_atlas_file)) {
    mspec$lh_surf_atlas <- lh_atlas_file
    mspec$rh_surf_atlas <- rh_atlas_file
  }
  else {
    message(sprintf('Atlas file %s or %s in the log file %s do not exist. Will try atlas in the modelspec file.',
                    lh_atlas_file, rh_atlas_file, svreg_log_file))
    if (is.null(mspec$subject$lh_atlas) || is.null(mspec$subject$rh_atlas))
      stop('One or more atlas files are not specified in the modelspec file. Please specify lh_atlas and rh_atlas in the modelspeec file', call. = FALSE)
    else {
      if (file.exists(mspec$subject$lh_atlas) && file.exists(mspec$subject$rh_atlas)) {
        message(sprintf('Using atlas files %s and %s for left and right hemispheres.', mspec$subject$lh_atlas, mspec$subject$rh_atlas))
        mspec$lh_surf_atlas <- mspec$subject$lh_atlas
        mspec$rh_surf_atlas <- mspec$subject$rh_atlas
      }
      else
        stop(sprintf('Atlas files %s or %s do not exist', mspec$subject$lh_atlas, mspec$subject$rh_atlas), call. = FALSE)
    }
  }
  return(mspec)
}

#' Get the tensor-based morphometry modelspec.ini file
#'
#' The modelspec file specifies the subject directory, paths to the atlas files,
#' the model specification, whether linear regerssion, correlation etc.
#' The modelspec file was used in the previous versions of bstr.
#' Currently, the user does not have to create this explicitly.
#' This functionality is still kept in case a future need arises to generate these automatically.
#' @param mspec path to the modelspec file
#' @param bs_atlas_path path to the atlas file
#' @param svreg_log_file file containing svreg log output
#'
get_tbm_atlas_files <- function(mspec, bs_atlas_path, svreg_log_file) {

  atlas_file <- file.path(dirname(bs_atlas_path), bs_file_formats$nii_atlas)
  if ( file.exists(atlas_file) )
    mspec$nii_atlas <- atlas_file
  else {
    message(sprintf('Atlas file %s in the log file %s does not exist. Will try atlas in the modelspec file...',
                    atlas_file, svreg_log_file))
    if (is.null(mspec$subject$atlas))
      stop('Atlas file is not specified in the modelspec file. It needs to be either obtained from svreg.log or specified explicitly in the modelspec.ini.', call. = FALSE)
    else {
      if (file.exists(mspec$subject$atlas))
        mspec$nii_atlas <- mspec$subject$atlas
      else
        stop(sprintf('Atlas file %s does not exist.', mspec$subject$atlas), call. = FALSE)
    }
  }

  if (is.null(mspec$subject$maskfile)) {
    message('Mask file not specified in the modelspec file. Will try the atlas directory to get the maskfile path.')
    maskfile <- file.path(dirname(bs_atlas_path), bs_file_formats$nii_maskfile)
    if (file.exists(maskfile)) {
      mspec$maskfile <- maskfile
      message(sprintf('Using maskfile %s.', maskfile))
    }
    else {
      mspec$maskfile <- NULL
      message(sprintf('Could not find mask file %s. Will skip masking the voxels for statistical analysis.', maskfile))
    }
  } else {
    if (file.exists(mspec$subject$maskfile))
      mspec$maskfile <- mspec$subject$maskfile
    else {
      mspec$maskfile <- NULL
      message(sprintf('Could not find mask file %s. Will skip masking the voxels for statistical analysis.', mspec$subject$maskfile))
    }
  }
  return(mspec)
}

#' Get the ROI modelspec.ini file
#'
#' The modelspec file specifies the subject directory, paths to the atlas files,
#' the model specification, whether linear regerssion, correlation etc.
#' The modelspec file was used in the previous versions of bstr.
#' Currently, the user does not have to create this explicitly.
#' This functionality is still kept in case a future need arises to generate these automatically.
#' @param mspec path to the modelspec file
#'
get_roi_specs <- function(mspec) {

  if (is.null(mspec$subject$roiid))
    stop('[subject] does not contain the roiid= field. Please specify a roiid or a comma separated list for multiple roiids.', call. = FALSE)
  else {
    mspec$roiid <- as.numeric(unlist(strsplit(mspec$subject$roiid, ',', fixed = TRUE)))
  }
  if (is.null(mspec$subject$roimeasure))
    stop('[subject] does not contain the roimeasure= field. It should either be gmthickness, gmvolume or area.', call. = FALSE)
  else {
    if (is.element(mspec$subject$roimeasure, c('gmthickness', 'gmvolume', 'area')) )
      mspec$roimeasure <- mspec$subject$roimeasure
    else
      stop(sprintf('Incorrect roimeasure=%s. roimeasure should either be gmthickness, gmvolume or area.', mspec$subject$roimeasure), call. = FALSE)
  }

  return(mspec)
}
