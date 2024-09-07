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

#' List of file formats used in BrainSuite
#' @export
## TODO: use closures for this in the future
bs_file_formats <- list(
  jacdet = '*.svreg.inv.map.jacdet.*.nii.gz',
  roi_txt = '.roiwise.stats.txt',
  wm_roi_txt = '.wm.stats.tsv',
  svreg_log = '*.svreg.log',
  surf_atlas_left = 'mri.left.mid.cortex.dfs',
  surf_atlas_right = 'mri.right.mid.cortex.dfs',
  nii_atlas = 'mri.bfc.nii.gz',
  nii_maskfile = 'mri.cerebrum.mask.nii.gz'
)

#' List ofROI types
#' @export
roi_types <- list(
  svreg_roi_types = c("gmthickness", "gmvolume", "area", "wmvolume"),
  swm_roi_types = c("swmFA", "swmMD", "swmRD", "swmAD"),
  wm_roi_types = c("FA", "FRT_GFA",	"L2",	"L3",	"MD",	"axial",	"mADC",	"radial")

)

#' List of atlas files used in BrainSuite
#' @export
bs_atlas_files <- list(
  atlas_BS1_tbm = 'svreg/BrainSuiteAtlas1/mri.bfc.nii.gz',
  atlas_BS1_mask_tbm = 'svreg/BrainSuiteAtlas1/mri.mask.nii.gz',
  atlas_BS1_mask_dba = 'svreg/BrainSuiteAtlas1/mri.cortex.dewisp.mask.nii.gz',
  lh_atlas_BS1_sba = 'svreg/BrainSuiteAtlas1/mri.left.mid.cortex.dfs',
  rh_atlas_BS1_sba = 'svreg/BrainSuiteAtlas1/mri.right.mid.cortex.dfs',

  atlas_BCIDNI_tbm = 'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.bfc.nii.gz',
  atlas_BCIDNI_mask_tbm = 'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.mask.nii.gz',
  atlas_BCIDNI_mask_dba = 'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.cortex.dewisp.mask.nii.gz',
  lh_atlas_BCIDNI_sba = 'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.left.mid.cortex.dfs',
  rh_atlas_BCIDNI_sba = 'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.right.mid.cortex.dfs'
)

#' List of binaries used in BrainSuite analysis
#' @export
bs_binary_files <- list(
  clustermap = 'bin/clustermap',
  statmap = 'bin/statmap'
)

#' List of suffixes for atlas files used in BrainSuite
bs_atlas_files_suffix <- list(

  atlas_custom_suffix_tbm = 'bfc.nii.gz',
  atlas_custom_mask_suffix_tbm = 'mask.nii.gz',
  atlas_custom_suffix_dba = 'bfc.nii.gz',
  atlas_custom_mask_suffix_dba = 'wm.mask.nii.gz'
)

analysis_type_list <- list(
  sba = 'sba',
  tbm = 'tbm',
  roi = 'roi',
  dba = 'dba',
  nca = 'nca'
)

bs_data_types <- list(
  surface = '.dfs',
  nifti_image = '.nii.gz'
)

#' List of statistical and other overlay types used in BrainSuite
bs_stat_overlays <- list(
  log_pvalues_adjusted = "log_pvalues_adjusted",
  tvalues_adjusted = "tvalues_adjusted",
  log_pvalues = "log_pvalues",
  tvalues = "tvalues",
  pvalues = "pvalues",
  corr_values = "corr_values",
  corr_values_masked_adjusted = "corr_values_masked_adjusted"
)

#' Mapping of stat overlays to label strings
bs_stat_overlays_mapping_to_label <- list(
  log_pvalues_adjusted = "adjusted log(p-values)",
  tvalues_adjusted = "adjusted t-values",
  log_pvalues = "log(p-values)",
  tvalues = "t-values",
  pvalues = "p-values",
  corr_values = "correlations",
  corr_values_masked_adjusted = "masked adjusted correlations"
)

#' Mapping of model_type to readable string
bs_model_type_to_readable_string <- list(
  bstr_anova = "ANOVA",
  bstr_lm = "Linear Fixed-effects Model",
  bstr_lmer = "Linear Mixed-effects Model",
  bstr_corr = "Correlation",
  pairedttest = "Paired T-test",
  unpairedttest = "Unpaired T-test"
)

#' Mapping of analysis_type to readable string
bs_analysis_type_to_readable_string <- list(
  tbm = "TBM",
  sba = "SBA",
  roi = "ROI",
  dba = "DBA"
)

create_rmd_report_title_str <- function(bstr_data, bstr_model) {

  if(bstr_data@analysis_type == "dba") {
    title_string <- paste(bstr:::bs_analysis_type_to_readable_string[[bstr_data@analysis_type]], " ", bstr_data@measure, " ",
                          bstr:::bs_model_type_to_readable_string[[bstr_model@model_type]], '-- ', sep="")
  }
  else {
    title_string <- paste(bstr:::bs_analysis_type_to_readable_string[[bstr_data@analysis_type]], " ",
                          bstr:::bs_model_type_to_readable_string[[bstr_model@model_type]], '-- ', sep="")
  }



  switch(bstr_model@model_type,
         bstr_anova = {
           title_string <- paste(title_string, "Main effect of", bstr_model@main_effect, "with covariates",
                 bstr_model@covariates)
         },
         bstr_lm = {
           title_string <- paste(title_string, "Main effect of", bstr_model@main_effect, "with covariates",
                 bstr_model@covariates)
         },
         bstr_lmer = {
           title_string <- paste(title_string, "Fixed effect of", bstr_model@main_effect, "with covariates",
                 bstr_model@covariates, ", Random effect of", bstr_model@group_var)
         },
         bstr_corr = {
           title_string <- paste(title_string, "with", bstr_model@corr_var)
         },
         pairedttest = {
           title_string <- paste(title_string, "between samples for", bstr_model@group_var)
         },
         unpairedttest = {
           title_string <- paste(title_string, "between samples for", bstr_model@group_var)
         }
  )

  return(title_string)
}

#' Returns the filename for the designated image
#' @param outdir string denoting the output directory
#' @param voxelcoord all outputted voxelcoordinates
#' @param overlay_name string denoting the name of the overlay
#' @param brain_sector_index numeric value denoting which brain slice is desired
#' @param voxelcoord_index numeric value denoting desired voxel coordinate
#'
get_render_image_filename <- function(outdir,voxelcoord, overlay_name, brain_sector_index, voxelcoord_index) {
  view_order <- c("sag","cor","ax")
  return(paste0(outdir, "/png_images_crosshairs/", view_order[brain_sector_index], voxelcoord[[voxelcoord_index]][brain_sector_index],"_",overlay_name,"_cluster",voxelcoord_index,".png"))
}
#' Generate the designated cortical surface file
#' @param hemi string denoting which hemisphere of the brain is desired
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
bs_surface_file_string <- function(hemi="left", smooth = 0) {

  if (smooth != 0)
    return(paste('atlas.pvc-thickness_0-6mm.', sprintf('smooth%2.1fmm.', smooth), hemi, '.mid.cortex.dfs', sep = ''))
  else
    return(paste('atlas.pvc-thickness_0-6mm.', hemi, '.mid.cortex.dfs', sep = ''))
}
#' Generate the designated tensor-based morphometry file
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
bs_volume_jacobian_file_string <- function(smooth = 0) {

  if (smooth != 0)
    return(paste('%s.svreg.inv.jacobian.', sprintf('smooth%2.1fmm.nii.gz', smooth), sep = ''))
  else
    return('%s.svreg.inv.jacobian.nii.gz')
}

#' Generate the designated BIDS compatible (T1w) tensor-based morphometry file
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
bs_BIDS_volume_jacobian_file_string <- function(smooth = 0) {

  if (smooth != 0)
    return(paste('%s_T1w.svreg.inv.jacobian.', sprintf('smooth%2.1fmm.nii.gz', smooth), sep = ''))
  else
    return('%s_T1w.svreg.inv.jacobian.nii.gz')
}

#' Generate the designated diffusion surface file
#' @param measure string designating the type of diffusion measure used
#' @param smooth numeric value designating the smoothing used (default is 0)
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#'
bs_diffusion_file_string <- function(measure = "FA", smooth = 0, eddy = TRUE) {

  valid_diffusion_measures <- c('FA', 'MD', 'axial', 'radial', 'mADC', 'FRT_GFA')
  if (! measure %in% valid_diffusion_measures)
    stop(sprintf('Invalid diffusion measure: %s. Valid measures are %s.', measure, paste(valid_diffusion_measures, collapse = ', ')))

  if (eddy == TRUE && smooth != 0)
    return(paste0('%s.dwi.RAS.correct.atlas.', measure, sprintf('.smooth%2.1fmm.nii.gz', smooth)))
  if (eddy == FALSE && smooth != 0)
    return(paste0('%s.dwi.RAS.atlas.', measure, sprintf('.smooth%2.1fmm.nii.gz', smooth)))
  if (eddy == TRUE && smooth == 0)
    return(paste0('%s.dwi.RAS.correct.atlas.', measure, '.nii.gz'))
  if (eddy == FALSE && smooth == 0)
    return(paste0('%s.dwi.RAS.atlas.', measure, '.nii.gz'))
}

#' Generate the BIDS designated diffusion surface file
#' @param measure string designating the type of diffusion measure used
#' @param smooth numeric value designating the smoothing used (default is 0)
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#'
bs_BIDS_diffusion_file_string <- function(measure = "FA", smooth = 0, eddy = TRUE) {

  valid_diffusion_measures <- c('FA', 'MD', 'axial', 'radial', 'mADC', 'FRT_GFA')
  if (! measure %in% valid_diffusion_measures)
    stop(sprintf('Invalid diffusion measure: %s. Valid measures are %s.', measure, paste(valid_diffusion_measures, collapse = ', ')))

  if (eddy == TRUE && smooth != 0)
    return(paste0('%s_dwi.dwi.RAS.correct.atlas.', measure, sprintf('.smooth%2.1fmm.nii.gz', smooth)))
  if (eddy == FALSE && smooth != 0)
    return(paste0('%s_dwi.dwi.RAS.atlas.', measure, sprintf('.smooth%2.1fmm.nii.gz', smooth)))
  if (eddy == TRUE && smooth == 0)
    return(paste0('%s_dwi.dwi.RAS.correct.atlas.', measure, '.nii.gz'))
  if (eddy == FALSE && smooth == 0)
    return(paste0('%s_dwi.dwi.RAS.atlas.', measure, '.nii.gz'))
}

#' Stops analysis if desired type of analysis is not a valid analysis type
#' @param analysis_type string denoting desired type of analysis to be performed
#'
get_bs_file_list <- function(analysis_type) {
  valid_analysis_types <- unlist(analysis_type_list, use.names = FALSE)
  if (!(analysis_type %in% valid_analysis_types)) {
    stop(sprintf('Valid brain analyses are %s', paste(unlist(analysis_type_list), collapse = ', ')),
         call. = FALSE)
  }


}
#' Returns a list of all ROI files for all subjects
#' @param bstr_data object of type `BstrData`
#'
get_roi_file_list <- function(bstr_data) {

  bids_flag_and_filelist <- check_bids_compatibility_and_get_filelist(bstr_data, type="roi", hemi="", smooth="", measure="")
  return(bids_flag_and_filelist$filelist)

}
#' Returns a list of the cortical surface files for all subjects
#' @param bstr_data object of type `BstrData`
#' @param hemi designates which hemisphere is of interest
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
get_sba_file_list <- function(bstr_data, hemi, smooth = 0) {

  bids_flag_and_filelist <- check_bids_compatibility_and_get_filelist(bstr_data, type="sba", hemi, smooth)
  return(bids_flag_and_filelist$filelist)
}
#' Returns a list of the tensor-based files for all subjects
#' @param bstr_data object of type `BstrData`
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
get_tbm_file_list <- function(bstr_data, smooth = 0) {

  bids_flag_and_filelist <- check_bids_compatibility_and_get_filelist(bstr_data, type="tbm", hemi="left", smooth)
  return(bids_flag_and_filelist$filelist)
}

#' Returns a list of the diffusion files for all subjects
#' @param bstr_data object of type `BstrData`
#' @param measure numeric value denoting the measures used to create the output
#' @param smooth numeric value designating the smoothing used (default is 0)
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#'
get_dba_file_list <- function(bstr_data, measure, smooth = 0, eddy = TRUE) {

  bids_flag_and_filelist <- check_bids_compatibility_and_get_filelist(bstr_data, type="dba", hemi="left", smooth, measure)
  return(bids_flag_and_filelist$filelist)
}
#' Read the BrainSuite atlas prefix from the atlas
#' @param atlas filepath for atlas
#'
get_brainsute_custom_volume_atlas_prefix <- function(atlas) {

  # Check if the parent directory exists for the atlas
  check_file_exists(dirname(atlas), raise_error = TRUE)

  # If atlas points to a nifti image, return the prefix of the atlas (everything until the bfc.nii.gz)
  if (substr(atlas, nchar(atlas) - 9, nchar(atlas)) == bs_atlas_files_suffix$atlas_custom_suffix_tbm)
    return(substr(atlas, 1, nchar(atlas) - 11))

  # If atlas points to a valid prefix, return atlas
  if (!identical(Sys.glob(file.path(paste(atlas, "*", bs_atlas_files_suffix$atlas_custom_suffix_tbm, sep = ""))), character(0)))
    return(atlas)

  # Otherwise raise an exception
  stop(sprintf('Invalid custom atlas prefix/path: %s', atlas), call. = FALSE)

}
#' Read the BrainSuite atlas path from the svreg log file.
#' @param logfile path to the svreg.log file present in an individual subject directory
#'
get_brainsuite_atlas_path_from_logfile <- function(logfile) {
  fid = file(logfile, "rt")
  log_lines <- readLines(fid, n=2)
  bs_atlas_path <- unlist(strsplit(log_lines[2], ' ', fixed = TRUE))[3]
  close(fid)
  return(bs_atlas_path)
}

#' Read the BrainSuite atlas identifier from the svreg log file.
#' @param logfile path to the svreg.log file present in an individual subject directory
#' @details Valid atlases are BrainSuiteAtlas1 or BCI-DNI_brain_atlas
#'
get_brainsuite_atlas_id_from_logfile <- function(logfile) {
  fid = file(logfile, "rt")
  log_lines <- readLines(fid, n=2)
  close(fid)
  if (grepl("BrainSuiteAtlas1", log_lines[2], fixed = TRUE))
    return("BrainSuiteAtlas1")
  else if (grepl("DNI_brain_atlas", log_lines[2], fixed = TRUE))
    return("BCI-DNI_brain_atlas")
  else
    stop(paste("Could not determine the BrainSuite atlas used for registration.\n",
               "Please check the log file ", logfile, ", and check if the subject directory is valid. If using a custom atlas, supply the path in the atlas= argument.", sep = ""), call. = FALSE)
}

#' Gets the BraisSuite svreg.log file's filepath
#' @param subjdir individual subject directory that the svreg.log file exists in
#' @param csv csv file for the svreg.log file
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @export
get_brainsuite_logfilename <- function(subjdir, csv, exclude_col="") {
  # Open the svreg.log file and get the atlas file name
  if ( identical(tools::file_ext(csv), 'tsv') | identical(tools::file_ext(csv), 'csv')  ) {
    demog <- read_demographics(csv, exclude_col)
  }
  # The first column has to contain subject IDs which are same as subject directories
  first_subjid <- demog[[1]][1]
  # Get atlas names from log files.
  svreg_log_file <- Sys.glob(file.path(subjdir, first_subjid, "*", '*.svreg.log')) #Extra * due to BIDS compatibility (extra dir anat or dwi)

  if (length(svreg_log_file) == 0){
    svreg_log_file <- Sys.glob(file.path(subjdir, first_subjid, '*.svreg.log'))  #If subject data is not BIDS compatible
  }
  # svreg_log_file <- file.path(subjdir, first_subjid, sprintf('%s.svreg.log', first_subjid))
  if (check_file_exists(svreg_log_file, raise_error = TRUE,
                    errmesg = sprintf('Could not find svreg.log in the subject directory %s/%s. Please check if the subject directory is valid.', subjdir, first_subjid)))
    return(tools::file_path_as_absolute(svreg_log_file))
}

#' Get the svreg log file for each subject
#' @param subjdir individual subject directory that the svreg.log file exists in
#' @param csv csv file for the svreg.log file
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @export
get_brainsuite_logfilename_for_all_subjects <- function(subjdir, csv, exclude_col="") {
  # Open the svreg.log file and get the atlas file name
  if ( identical(tools::file_ext(csv), 'tsv') | identical(tools::file_ext(csv), 'csv')  ) {
    demog <- read_demographics(csv, exclude_col)
  }
  # The first column has to contain subject IDs which are same as subject directories
  first_subjid <- demog[[1]][1]
  # Get atlas names from log files.
  svreg_log_files <- Sys.glob(file.path(subjdir, demog[[1]], "*", '*.svreg.log'))

  if (length(svreg_log_files) == 0){
    svreg_log_files <- Sys.glob(file.path(subjdir, demog[[1]], '*.svreg.log'))
  }
  # svreg_log_file <- file.path(subjdir, first_subjid, sprintf('%s.svreg.log', first_subjid))
  if (check_multiple_files_exists(svreg_log_files, errmesg = 'Could not find svreg.log in a few subject directories.')) {
    return(svreg_log_files)
  }
}

#' Get the desired cortical surface atlas
#' @param brainsuite_atlas_id individual subject directory that the svreg.log file exists in
#' @param hemi designates which hemisphere of the brain
#'
get_sba_atlas <- function(brainsuite_atlas_id, hemi) {

  if (! brainsuite_atlas_id %in% c("BrainSuiteAtlas1", "BCI-DNI_brain_atlas"))
    stop('Valid values for hemi are BrainSuiteAtlas1 or BCI-DNI_brain_atlas.', call. = FALSE)

  if (! hemi %in% c("left", "right"))
    stop('Valid values for hemi are left or right.', call. = FALSE)

  brainsuite_install_path <- get_brainsuite_install_path()
  if (brainsuite_atlas_id == "BrainSuiteAtlas1") {
    if (hemi == "left") {
      lh_surf_atlas <- file.path(brainsuite_install_path, bs_atlas_files$lh_atlas_BS1_sba)
      check_file_exists(lh_surf_atlas, raise_error = TRUE)
      return(lh_surf_atlas)
    }
    else if (hemi == "right") {
      rh_surf_atlas <- file.path(brainsuite_install_path, bs_atlas_files$rh_atlas_BS1_sba)
      check_file_exists(rh_surf_atlas, raise_error = TRUE)
      return(rh_surf_atlas)
    }
  }

  if (brainsuite_atlas_id == "BCI-DNI_brain_atlas") {
    if (hemi == "left") {
      lh_surf_atlas <- file.path(brainsuite_install_path, bs_atlas_files$lh_atlas_BCIDNI_sba)
      check_file_exists(lh_surf_atlas, raise_error = TRUE)
      return(lh_surf_atlas)
    }
    else if (hemi == "right") {
      rh_surf_atlas <- file.path(brainsuite_install_path, bs_atlas_files$rh_atlas_BCIDNI_sba)
      check_file_exists(rh_surf_atlas, raise_error = TRUE)
      return(rh_surf_atlas)
    }
  }
}
#' Get the tensor based morphometry atlas and mask
#' @param brainsuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_tbm_atlas_and_mask <- function(brainsuite_atlas_id) {

  if (! brainsuite_atlas_id %in% c("BrainSuiteAtlas1", "BCI-DNI_brain_atlas"))
    stop('Valid values for atlas are BrainSuiteAtlas1 or BCI-DNI_brain_atlas.', call. = FALSE)
  brainsuite_install_path <- get_brainsuite_install_path()
  if (brainsuite_atlas_id == "BrainSuiteAtlas1") {
    nii_atlas <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BS1_tbm)
    nii_atlas_mask <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BS1_mask_tbm)
  }
  if (brainsuite_atlas_id == "BCI-DNI_brain_atlas") {
    nii_atlas <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BCIDNI_tbm)
    nii_atlas_mask <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BCIDNI_mask_tbm)
  }
  check_file_exists(nii_atlas, raise_error = TRUE)
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(list("nii_atlas" = nii_atlas, "nii_atlas_mask" = nii_atlas_mask))
}

#' Check and return the maskfile and the atlas file
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param maskfile filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
check_tbm_atlas_and_mask <- function(subjdir="", csv="", atlas="", maskfile="", exclude_col) {

  if (maskfile == ""  && atlas == "") {
    brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv, exclude_col))
    tbm_atlas_and_mask <- get_tbm_atlas_and_mask(brainsuite_atlas_id)
  }
  else {
    if (maskfile != "" && atlas == "") {
      brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv, exclude_col))
      check_file_exists(maskfile, raise_error = TRUE)
      tbm_atlas_and_mask$nii_atlas_mask <- maskfile
      tbm_atlas_and_mask$nii_atlas <- get_tbm_atlas(brainsuite_atlas_id)
    }
    else if (maskfile== "" && atlas != "") {
      brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv, exclude_col))
      check_file_exists(atlas, raise_error = TRUE)
      tbm_atlas_and_mask$nii_atlas <- atlas
      tbm_atlas_and_mask$nii_atlas_mask <- get_tbm_mask(brainsuite_atlas_id)
    }
    else if (maskfile!= "" && atlas != "") {
      check_file_exists(maskfile, raise_error = TRUE)
      check_file_exists(atlas, raise_error = TRUE)
      tbm_atlas_and_mask <- list("nii_atlas" = atlas, "nii_atlas_mask" = maskfile)
    }
  }

  return(tbm_atlas_and_mask)

}



#' Get the BrainSuite tensor based morphometry atlas
#' @param brainsuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_tbm_atlas <- function(brainsuite_atlas_id) {

  if (! brainsuite_atlas_id %in% c("BrainSuiteAtlas1", "BCI-DNI_brain_atlas"))
    stop('Valid values for atlas are BrainSuiteAtlas1 or BCI-DNI_brain_atlas.', call. = FALSE)
  brainsuite_install_path <- get_brainsuite_install_path()
  if (brainsuite_atlas_id == "BrainSuiteAtlas1") {
    nii_atlas <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BS1_tbm)
  }
  if (brainsuite_atlas_id == "BCI-DNI_brain_atlas") {
    nii_atlas <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BCIDNI_tbm)
  }
  check_file_exists(nii_atlas, raise_error = TRUE)
  return(nii_atlas)
}

#' Get the BrainSuite tensor based morphometry mask
#' @param brainsuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_tbm_mask <- function(brainsuite_atlas_id) {

  if (! brainsuite_atlas_id %in% c("BrainSuiteAtlas1", "BCI-DNI_brain_atlas"))
    stop('Valid values for atlas are BrainSuiteAtlas1 or BCI-DNI_brain_atlas.', call. = FALSE)
  brainsuite_install_path <- get_brainsuite_install_path()
  if (brainsuite_atlas_id == "BrainSuiteAtlas1") {
    nii_atlas_mask <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BS1_mask_tbm)
  }
  if (brainsuite_atlas_id == "BCI-DNI_brain_atlas") {
    nii_atlas_mask <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BCIDNI_mask_tbm)
  }
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(nii_atlas_mask)
}


#' Check that cortical surface atlas exists
#' @param atlas filepath for sba atlas
#' @param maskfile filepath for the atlas mask file
#'
get_custom_sba_atlas_and_mask <- function(atlas, maskfile="") {
  #TODO Only the atlas file is implemented
  check_file_exists(atlas, raise_error = TRUE)
  return(atlas)
}
#' Get the custom tensor based morphometry atlas and mask
#' @param brainsuite_custom_atlas_prefix file prefix for atlas
#'
get_custom_tbm_atlas_and_mask <- function(brainsuite_custom_atlas_prefix) {

  brainsuite_custom_atlas_prefix <- get_brainsute_custom_volume_atlas_prefix(brainsuite_custom_atlas_prefix)
  nii_atlas <- paste(brainsuite_custom_atlas_prefix, ".", bs_atlas_files_suffix$atlas_custom_suffix_tbm, sep="")
  nii_atlas_mask <- paste(brainsuite_custom_atlas_prefix, ".", bs_atlas_files_suffix$atlas_custom_mask_suffix_tbm, sep="")
  check_file_exists(nii_atlas, raise_error = TRUE)
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(list("nii_atlas" = nii_atlas, "nii_atlas_mask" = nii_atlas_mask))
}
#' Get the diffusion atlas and mask
#' @param brainsuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_dba_atlas_and_mask <- function(brainsuite_atlas_id) {

  if (! brainsuite_atlas_id %in% c("BrainSuiteAtlas1", "BCI-DNI_brain_atlas"))
    stop('Valid values for atlas are BrainSuiteAtlas1 or BCI-DNI_brain_atlas.', call. = FALSE)
  brainsuite_install_path <- get_brainsuite_install_path()
  if (brainsuite_atlas_id == "BrainSuiteAtlas1") {
    nii_atlas <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BS1_tbm)
    nii_atlas_mask <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BS1_mask_dba)
  }
  if (brainsuite_atlas_id == "BCI-DNI_brain_atlas") {
    nii_atlas <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BCIDNI_tbm)
    nii_atlas_mask <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BCIDNI_mask_dba)
  }
  check_file_exists(nii_atlas, raise_error = TRUE)
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(list("nii_atlas" = nii_atlas, "nii_atlas_mask" = nii_atlas_mask))
}

#' Get the BrainSuite tensor based morphometry atlas
#' @param brainsuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_dba_atlas <- function(brainsuite_atlas_id) {

  if (! brainsuite_atlas_id %in% c("BrainSuiteAtlas1", "BCI-DNI_brain_atlas"))
    stop('Valid values for atlas are BrainSuiteAtlas1 or BCI-DNI_brain_atlas.', call. = FALSE)
  brainsuite_install_path <- get_brainsuite_install_path()
  if (brainsuite_atlas_id == "BrainSuiteAtlas1") {
    nii_atlas <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BS1_tbm)
  }
  if (brainsuite_atlas_id == "BCI-DNI_brain_atlas") {
    nii_atlas <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BCIDNI_tbm)
  }
  check_file_exists(nii_atlas, raise_error = TRUE)
  return(nii_atlas)
}

#' Get the BrainSuite tensor based morphometry mask
#' @param brainsuite_atlas_id individual subject directory that the svreg.log file exists in
#'
get_dba_mask <- function(brainsuite_atlas_id) {

  if (! brainsuite_atlas_id %in% c("BrainSuiteAtlas1", "BCI-DNI_brain_atlas"))
    stop('Valid values for atlas are BrainSuiteAtlas1 or BCI-DNI_brain_atlas.', call. = FALSE)
  brainsuite_install_path <- get_brainsuite_install_path()
  if (brainsuite_atlas_id == "BrainSuiteAtlas1") {
    nii_atlas_mask <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BS1_mask_dba)
  }
  if (brainsuite_atlas_id == "BCI-DNI_brain_atlas") {
    nii_atlas_mask <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BCIDNI_mask_dba)
  }
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(nii_atlas_mask)
}



#' Get the diffusion custom atlas and mask
#' @param brainsuite_custom_atlas_prefix file prefix for atlas
#'
get_custom_dba_atlas_and_mask <- function(brainsuite_custom_atlas_prefix) {

  brainsuite_custom_atlas_prefix <- get_brainsute_custom_volume_atlas_prefix(brainsuite_custom_atlas_prefix)
  nii_atlas <- paste(brainsuite_custom_atlas_prefix, ".", bs_atlas_files_suffix$atlas_custom_suffix_dba, sep = "")
  nii_atlas_mask <- paste(brainsuite_custom_atlas_prefix, ".", bs_atlas_files_suffix$atlas_custom_mask_suffix_dba, sep = "")
  check_file_exists(nii_atlas, raise_error = TRUE)
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(list("nii_atlas" = nii_atlas, "nii_atlas_mask" = nii_atlas_mask))
}
#' Reads the demographics from an inputted csv file
#' @param csvfile csv file containing the demographics
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
read_demographics <- function(csvfile, exclude_col="") {

  switch(tools::file_ext(csvfile),
         "tsv" = {demog <- read.table(file = csvfile, sep = "\t", header = T)},
         "csv" = {demog <- read.csv(csvfile)})
  # demo <- read.csv(csvfile)
  colnames(demog)[1] <- "subjID"
  if (exclude_col != "") {
    if (! exclude_col %in% colnames(demog))
      stop(sprintf("Exclude column specified as %s does not exist in %s.", exclude_col, csvfile), call. = FALSE)
    demog <- demog[demog[, exclude_col] == 0,]
  }
  return(demog)
}

check_bids_compatibility_and_get_filelist <- function(bstr_data, type="sba", hemi, smooth = 0, measure="FA", eddy = TRUE, roimeas = "gmthickness") {

  valid_types <- c("sba", "tbm", "roi", "dba")
  if (! type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  if (type == "sba") {
    filelist <- file.path(bstr_data@subjdir, bstr_data@demographics$subjID, bs_surface_file_string(hemi, smooth))
    bids_filelist <- file.path(bstr_data@subjdir, bstr_data@demographics$subjID, 'anat', bs_surface_file_string(hemi, smooth))
  }
  else if (type == "tbm") {
    filelist <- file.path(bstr_data@subjdir, bstr_data@demographics$subjID, sprintf(bs_volume_jacobian_file_string(smooth), bstr_data@demographics$subjID))
    bids_filelist <- file.path(bstr_data@subjdir, bstr_data@demographics$subjID, 'anat', sprintf(bs_BIDS_volume_jacobian_file_string(smooth), bstr_data@demographics$subjID))
  }
  else if (type == "dba") {
    valid_dba_measures <- c('FA', 'MD', 'axial', 'radial', 'mADC', 'FRT_GFA')
    if (! measure %in% valid_dba_measures) {
      stop(sprintf("Valid dba measures are %s.", paste(valid_dba_measures, collapse = ', ')), call. = FALSE)
    }
    filelist <- file.path(bstr_data@subjdir, bstr_data@demographics$subjID, sprintf(bs_diffusion_file_string(measure, smooth, eddy), bstr_data@demographics$subjID))
    bids_filelist <- file.path(bstr_data@subjdir, bstr_data@demographics$subjID, 'dwi', sprintf(bs_BIDS_diffusion_file_string(measure, smooth, eddy), bstr_data@demographics$subjID))

  }
  else if (type == "roi") {
    valid_roi_measures <- c(roi_types$svreg_roi_types, roi_types$swm_roi_types, roi_types$swm_roi_types)
    if (roimeas %in% roi_types$svreg_roi_types || roimeas %in% roi_types$swm_roi_types) {
      bids_filelist <- file.path(bstr_data@subjdir, bstr_data@demographics$subjID, 'anat',
                                     sprintf('%s%s%s', bstr_data@demographics$subjID, '_T1w', bs_file_formats$roi_txt))
      filelist <- file.path(bstr_data@subjdir, bstr_data@demographics$subjID,
                                 sprintf('%s%s%s', bstr_data@demographics$subjID, '_T1w', bs_file_formats$roi_txt))

    }
    else if (roimeas %in% roi_types$wm_roi_types) {
      bids_filelist <- file.path(bstr_data@subjdir, bstr_data@demographics$subjID, 'dwi',
                                     sprintf('%s%s%s', bstr_data@demographics$subjID, '_T1w', bs_file_formats$wm_roi_txt))
      filelist <- file.path(bstr_data@subjdir, bstr_data@demographics$subjID,
                                 sprintf('%s%s%s', bstr_data@demographics$subjID, '_T1w', bs_file_formats$wm_roi_txt))

    }
  }

  # Check BIDS compatibility
  if (any(file.exists(filelist))) {
    bids_compatible = FALSE
  }
  else if (any(file.exists(bids_filelist))) {
    bids_compatible = TRUE
    filelist <- bids_filelist
  }
  else {
    message('Following subjects have missing files')
    print(filelist[which(!file.exists(filelist))], row.names = FALSE)
    stop('\nCheck if SVREG was run succesfully on all the subjects. Optionally, check the smoothing level (smooth= under [subject]).\nIt is possible that surface or volume files at the specified smoothing level do not exist.',
         call. = FALSE)

#    stop('\nCould not understand the directory hierarchy. Check if svreg was run succesfully on all the subjects. Also check the smoothing level (smooth= under [subject]).\nIt is possible that surface files at the specified smoothing level do not exist.',
#         call. = FALSE)
  }

  # Check if all files exist
  if ( !all(file.exists(filelist))) {
    message('Following subjects have missing files')
    print(filelist[which(!file.exists(filelist))], row.names = FALSE)
    stop('\nCheck if SVREG was run succesfully on all the subjects. Optionally, check the smoothing level (smooth= under [subject]).\nIt is possible that surface or volume files at the specified smoothing level do not exist.',
         call. = FALSE)
  }
  return (list("bids_compatible" = bids_compatible, "filelist" = filelist))

}
