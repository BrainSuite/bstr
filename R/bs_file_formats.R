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

#' List of file formats used in BrainSuite
#' @export
## TODO: use closures for this in the future
bs_file_formats <- list(
  jacdet = '*.svreg.inv.map.jacdet.*.nii.gz',
  roi_txt = '.roiwise.stats.txt',
  svreg_log = '*.svreg.log',
  surf_atlas_left = 'mri.left.mid.cortex.dfs',
  surf_atlas_right = 'mri.right.mid.cortex.dfs',
  nii_atlas = 'mri.bfc.nii.gz',
  nii_maskfile = 'mri.cerebrum.mask.nii.gz'
)

#' List of atlas files used in BrainSuite
#' @export
bs_atlas_files <- list(
  atlas_BS1_tbm = 'svreg/BrainSuiteAtlas1/mri.bfc.nii.gz',
  atlas_BS1_mask_tbm = 'svreg/BrainSuiteAtlas1/mri.mask.nii.gz',
  atlas_BS1_mask_dbm = 'svreg/BrainSuiteAtlas1/mri.cortex.dewisp.mask.nii.gz',
  lh_atlas_BS1_cbm = 'svreg/BrainSuiteAtlas1/mri.left.mid.cortex.dfs',
  rh_atlas_BS1_cbm = 'svreg/BrainSuiteAtlas1/mri.right.mid.cortex.dfs',

  atlas_BCIDNI_tbm = 'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.bfc.nii.gz',
  atlas_BCIDNI_mask_tbm = 'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.mask.nii.gz',
  atlas_BCIDNI_mask_dbm = 'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.cortex.dewisp.mask.nii.gz',
  lh_atlas_BCIDNI_cbm = 'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.left.mid.cortex.dfs',
  rh_atlas_BCIDNI_cbm = 'svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.right.mid.cortex.dfs'
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
  atlas_custom_suffix_dbm = 'bfc.nii.gz',
  atlas_custom_mask_suffix_dbm = 'wm.mask.nii.gz'
)

analysis_type_list <- list(
  cbm = 'cbm',
  tbm = 'tbm',
  roi = 'roi',
  dbm = 'dbm',
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
#' @param bss_data object of type \code{BssData}
#'
get_roi_file_list <- function(bss_data) {

  roi_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID,
                            sprintf('%s%s', bss_data@demographics$subjID, bs_file_formats$roi_txt))
  roi_bids_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, 'anat',
                            sprintf('%s%s', bss_data@demographics$subjID, bs_file_formats$roi_txt))
  # Check if all subjects have roi text files
  if ( !all(file.exists(roi_filelist)) & !all(file.exists(roi_bids_filelist))) {
    message('Following subjects have missing roi text files')
    print(roi_filelist[which(!file.exists(roi_filelist))], row.names = FALSE)
    stop('\nCheck if svreg was run succesfully and if roiwise.txt are present in all the subjects.', call. = FALSE)
  }
  if(all(file.exists(roi_filelist)))
    return(roi_filelist)
  else
    return(roi_bids_filelist)
}
#' Returns a list of the cortical surface files for all subjects
#' @param bss_data object of type \code{BssData}
#' @param hemi designates which hemisphere is of interest
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
get_cbm_file_list <- function(bss_data, hemi, smooth = 0) {

  cbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, bs_surface_file_string(hemi, smooth))
  cbm_bids_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, 'anat', bs_surface_file_string(hemi, smooth))
  # Check if all subjects have dfs files
  if ( !all(file.exists(cbm_filelist)) & !all(file.exists(cbm_bids_filelist))) {
    message('Following subjects have missing dfs files')
    print(cbm_filelist[which(!file.exists(cbm_filelist))], row.names = FALSE)
    stop('\nCheck if svreg was run succesfully on all the subjects. Also check the smoothing level (smooth= under [subject]).\nIt is possible that surface files at the specified smoothing level do not exist.',
         call. = FALSE)
  }
  if(all(file.exists(cbm_filelist)))
    return(cbm_filelist)
  else
    return(cbm_bids_filelist)
}
#' Returns a list of the tensor-based files for all subjects
#' @param bss_data object of type \code{BssData}
#' @param smooth numeric value designating the smoothing used (default is 0)
#'
get_tbm_file_list <- function(bss_data, smooth = 0) {

  tbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, sprintf(bs_volume_jacobian_file_string(smooth), bss_data@demographics$subjID))
  tbm_bids_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, 'anat', sprintf(bs_BIDS_volume_jacobian_file_string(smooth), bss_data@demographics$subjID))
  # Check if all subjects have nii.gz files
  if ( !all(file.exists(tbm_filelist)) & !all(file.exists(tbm_bids_filelist))) {
    message('Following subjects have missing nii.gz files')
    print(tbm_filelist[which(!file.exists(tbm_filelist))], row.names = FALSE)
    stop('\nCheck if svreg was run succesfully and if jacobian* files exist for all the subjects.\nAlso check if smoothing was performed.', call. = FALSE)
  }
  if(all(file.exists(tbm_filelist)))
    return(tbm_filelist)
  else
    return(tbm_bids_filelist)
}

#' Returns a list of the diffusion files for all subjects
#' @param bss_data object of type \code{BssData}
#' @param measure numeric value denoting the measures used to create the output
#' @param smooth numeric value designating the smoothing used (default is 0)
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#'
get_dbm_file_list <- function(bss_data, measure, smooth = 0, eddy = TRUE) {

  valid_dbm_measures <- c('FA', 'MD', 'axial', 'radial', 'mADC', 'FRT_GFA')
  if (! measure %in% valid_dbm_measures) {
    stop(sprintf("Valid dbm measures are %s.", paste(valid_dbm_measures, collapse = ', ')), call. = FALSE)
  }
  dbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, sprintf(bs_diffusion_file_string(measure, smooth, eddy), bss_data@demographics$subjID))
  dbm_bids_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, 'dwi', sprintf(bs_BIDS_diffusion_file_string(measure, smooth, eddy), bss_data@demographics$subjID))
  # Check if all subjects have nii.gz files
  if ( !all(file.exists(dbm_filelist)) & !all(file.exists(dbm_bids_filelist))) {
    message('Following subjects have missing nii.gz files')
    print(dbm_filelist[which(!file.exists(dbm_filelist))], row.names = FALSE)
    stop('\nCheck if svreg_apply_map was run succesfully and if the *.dwi.*.atlas.*.nii.gz files exist for all the subjects.\nAlso check if smoothing was performed.', call. = FALSE)
  }
  if(all(file.exists(dbm_filelist)))
    return(dbm_filelist)
  else
    return(dbm_bids_filelist)
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
  if (grepl("/BrainSuiteAtlas1/", log_lines[2], fixed = TRUE))
    return("BrainSuiteAtlas1")
  else if (grepl("/BCI-DNI_brain_atlas/", log_lines[2], fixed = TRUE))
    return("BCI-DNI_brain_atlas")
  else
    stop(paste("Could not determine the BrainSuite atlas used for registration.\n",
               "Please check the log file ", logfile, ", and check if the subject directory is valid. If using a custom atlas, supply the path in the atlas= argument.", sep = ""), call. = FALSE)
}

#' Gets the BraisSuite svreg.log file's filepath
#' @param subjdir individual subject directory that the svreg.log file exists in
#' @param csv csv file for the svreg.log file
#' @export
get_brainsuite_logfilename <- function(subjdir, csv) {
  # Open the svreg.log file and get the atlas file name
  if ( identical(tools::file_ext(csv), 'tsv') | identical(tools::file_ext(csv), 'csv')  ) {
    demo <- read_demographics(csv)
  }
  # The first column has to contain subject IDs which are same as subject directories
  first_subjid <- demo[[1]][1]
  # Get atlas names from log files.
  svreg_log_file <- Sys.glob(file.path(subjdir, first_subjid, "*", '*.svreg.log'))

  if (length(svreg_log_file) == 0){
    svreg_log_file <- Sys.glob(file.path(subjdir, first_subjid, '*.svreg.log'))
  }
  # svreg_log_file <- file.path(subjdir, first_subjid, sprintf('%s.svreg.log', first_subjid))
  if (check_file_exists(svreg_log_file, raise_error = TRUE,
                    errmesg = sprintf('Could not find svreg.log in the subject directory %s/%s. Please check if the subject directory is valid.', subjdir, first_subjid)))
    return(tools::file_path_as_absolute(svreg_log_file))
}

#' Get the svreg log file for each subject
#' @param subjdir individual subject directory that the svreg.log file exists in
#' @param csv csv file for the svreg.log file
#' @export
get_brainsuite_logfilename_for_all_subjects <- function(subjdir, csv) {
  # Open the svreg.log file and get the atlas file name
  if ( identical(tools::file_ext(csv), 'tsv') | identical(tools::file_ext(csv), 'csv')  ) {
    demo <- read_demographics(csv)
  }
  # The first column has to contain subject IDs which are same as subject directories
  first_subjid <- demo[[1]][1]
  # Get atlas names from log files.
  svreg_log_files <- Sys.glob(file.path(subjdir, demo[[1]], "*", '*.svreg.log'))

  if (length(svreg_log_files) == 0){
    svreg_log_files <- Sys.glob(file.path(subjdir, demo[[1]], '*.svreg.log'))
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
get_cbm_atlas <- function(brainsuite_atlas_id, hemi) {

  if (! brainsuite_atlas_id %in% c("BrainSuiteAtlas1", "BCI-DNI_brain_atlas"))
    stop('Valid values for hemi are BrainSuiteAtlas1 or BCI-DNI_brain_atlas.', call. = FALSE)

  if (! hemi %in% c("left", "right"))
    stop('Valid values for hemi are left or right.', call. = FALSE)

  brainsuite_install_path <- get_brainsuite_install_path()
  if (brainsuite_atlas_id == "BrainSuiteAtlas1") {
    if (hemi == "left") {
      lh_surf_atlas <- file.path(brainsuite_install_path, bs_atlas_files$lh_atlas_BS1_cbm)
      check_file_exists(lh_surf_atlas, raise_error = TRUE)
      return(lh_surf_atlas)
    }
    else if (hemi == "right") {
      rh_surf_atlas <- file.path(brainsuite_install_path, bs_atlas_files$rh_atlas_BS1_cbm)
      check_file_exists(rh_surf_atlas, raise_error = TRUE)
      return(rh_surf_atlas)
    }
  }

  if (brainsuite_atlas_id == "BCI-DNI_brain_atlas") {
    if (hemi == "left") {
      lh_surf_atlas <- file.path(brainsuite_install_path, bs_atlas_files$lh_atlas_BCIDNI_cbm)
      check_file_exists(lh_surf_atlas, raise_error = TRUE)
      return(lh_surf_atlas)
    }
    else if (hemi == "right") {
      rh_surf_atlas <- file.path(brainsuite_install_path, bs_atlas_files$rh_atlas_BCIDNI_cbm)
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

#' Check that cortical surface atlas exists
#' @param atlas filepath for cbm atlas
#' @param maskfile filepath for the atlas mask file
#'
get_custom_cbm_atlas_and_mask <- function(atlas, maskfile="") {
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
get_dbm_atlas_and_mask <- function(brainsuite_atlas_id) {

  if (! brainsuite_atlas_id %in% c("BrainSuiteAtlas1", "BCI-DNI_brain_atlas"))
    stop('Valid values for atlas are BrainSuiteAtlas1 or BCI-DNI_brain_atlas.', call. = FALSE)
  brainsuite_install_path <- get_brainsuite_install_path()
  if (brainsuite_atlas_id == "BrainSuiteAtlas1") {
    nii_atlas <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BS1_tbm)
    nii_atlas_mask <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BS1_mask_dbm)
  }
  if (brainsuite_atlas_id == "BCI-DNI_brain_atlas") {
    nii_atlas <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BCIDNI_tbm)
    nii_atlas_mask <- file.path(brainsuite_install_path, bs_atlas_files$atlas_BCIDNI_mask_dbm)
  }
  check_file_exists(nii_atlas, raise_error = TRUE)
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(list("nii_atlas" = nii_atlas, "nii_atlas_mask" = nii_atlas_mask))
}
#' Get the diffusion custom atlas and mask
#' @param brainsuite_custom_atlas_prefix file prefix for atlas
#'
get_custom_dbm_atlas_and_mask <- function(brainsuite_custom_atlas_prefix) {

  brainsuite_custom_atlas_prefix <- get_brainsute_custom_volume_atlas_prefix(brainsuite_custom_atlas_prefix)
  nii_atlas <- paste(brainsuite_custom_atlas_prefix, ".", bs_atlas_files_suffix$atlas_custom_suffix_dbm, sep = "")
  nii_atlas_mask <- paste(brainsuite_custom_atlas_prefix, ".", bs_atlas_files_suffix$atlas_custom_mask_suffix_dbm, sep = "")
  check_file_exists(nii_atlas, raise_error = TRUE)
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(list("nii_atlas" = nii_atlas, "nii_atlas_mask" = nii_atlas_mask))
}
#' Reads the demographics from an inputted csv file
#' @param csvfile csv file containing the demographics
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
read_demographics <- function(csvfile, exclude_col="") {

  switch(tools::file_ext(csvfile),
         "tsv" = {demo <- read.table(file = csvfile, sep = "\t", header = T)},
         "csv" = {demo <- read.csv(csvfile)})
  # demo <- read.csv(csvfile)
  colnames(demo)[1] <- "subjID"
  if (exclude_col != "") {
    if (! exclude_col %in% colnames(demo))
      stop(sprintf("Exclude column specified as %s does not exist in %s.", exclude_col, csvfile), call. = FALSE)
    demo <- subset(demo, demo[[exclude_col]] == 1)
  }
  return(demo)
}
