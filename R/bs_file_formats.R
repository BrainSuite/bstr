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
  corr_values_adjusted = "corr_values_adjusted"
)

get_render_image_filename <- function(outdir,voxelcoord, overlay_name, brain_sector_index, voxelcoord_index) {
  view_order <- c("sag","cor","ax")
  return(paste0(outdir, "/png_images_crosshairs/", view_order[brain_sector_index], voxelcoord[[voxelcoord_index]][brain_sector_index],"_",overlay_name,"_cluster",voxelcoord_index,".png"))
}

bs_surface_file_string <- function(hemi="left", smooth = 0) {

  if (smooth != 0)
    return(paste('atlas.pvc-thickness_0-6mm.', sprintf('smooth%2.1fmm.', smooth), hemi, '.mid.cortex.dfs', sep = ''))
  else
    return(paste('atlas.pvc-thickness_0-6mm.', hemi, '.mid.cortex.dfs', sep = ''))
}

bs_volume_jacobian_file_string <- function(smooth = 0) {

  if (smooth != 0)
    return(paste('%s.svreg.inv.jacobian.', sprintf('smooth%2.1fmm.nii.gz', smooth), sep = ''))
  else
    return('%s.svreg.inv.jacobian.nii.gz')
}

bs_diffusion_file_string <- function(measure = "FA", smooth = 0, eddy = TRUE) {

  valid_diffusion_measures <- c('FA', 'MD', 'AD', 'RD', 'ADC', 'GFA')
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

get_bs_file_list <- function(analysis_type) {
  valid_analysis_types <- unlist(analysis_type_list, use.names = FALSE)
  if (!(analysis_type %in% valid_analysis_types)) {
    stop(sprintf('Valid brain analyses are %s', paste(unlist(analysis_type_list), collapse = ', ')),
         call. = FALSE)
  }


}

get_roi_file_list <- function(bss_data) {

  roi_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID,
                            sprintf('%s%s', bss_data@demographics$subjID, bs_file_formats$roi_txt))

  # Check if all subjects have roi text files
  if ( !all(file.exists(roi_filelist)) ) {
    message('Following subjects have missing roi text files')
    print(roi_filelist[which(!file.exists(roi_filelist))], row.names = FALSE)
    stop('\nCheck if svreg was run succesfully and if roiwise.txt are present in all the subjects.', call. = FALSE)
  }
  return(roi_filelist)
}

get_cbm_file_list <- function(bss_data, hemi, smooth = 0) {

  cbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, bs_surface_file_string(hemi, smooth))
  # Check if all subjects have dfs files
  if ( !all(file.exists(cbm_filelist)) ) {
    message('Following subjects have missing dfs files')
    print(cbm_filelist[which(!file.exists(cbm_filelist))], row.names = FALSE)
    stop('\nCheck if svreg was run succesfully on all the subjects. Also check the smoothing level (smooth= under [subject]).\nIt is possible that surface files at the specified smoothing level do not exist.',
         call. = FALSE)
  }
  return(cbm_filelist)
}

get_tbm_file_list <- function(bss_data, smooth = 0) {

  tbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, sprintf(bs_volume_jacobian_file_string(smooth), bss_data@demographics$subjID))
  # Check if all subjects have nii.gz files
  if ( !all(file.exists(tbm_filelist)) ) {
    message('Following subjects have missing nii.gz files')
    print(tbm_filelist[which(!file.exists(tbm_filelist))], row.names = FALSE)
    stop('\nCheck if svreg was run succesfully and if jacobian* files exist for all the subjects.\nAlso check if smoothing was performed.', call. = FALSE)
  }
  return(tbm_filelist)
}

get_dbm_file_list <- function(bss_data, measure, smooth = 0, eddy = TRUE) {

  dbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, sprintf(bs_diffusion_file_string(measure, smooth, eddy), bss_data@demographics$subjID))
  # Check if all subjects have nii.gz files
  if ( !all(file.exists(dbm_filelist)) ) {
    message('Following subjects have missing nii.gz files')
    print(dbm_filelist[which(!file.exists(dbm_filelist))], row.names = FALSE)
    stop('\nCheck if svreg_apply_map was run succesfully and if the *.dwi.*.atlas.*.nii.gz files exist for all the subjects.\nAlso check if smoothing was performed.', call. = FALSE)
  }
  return(dbm_filelist)
}

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
#' @export
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

get_custom_cbm_atlas_and_mask <- function(atlas, maskfile="") {
  #TODO Only the atlas file is implemented
  check_file_exists(atlas, raise_error = TRUE)
  return(atlas)
}

get_custom_tbm_atlas_and_mask <- function(brainsuite_custom_atlas_prefix) {

  brainsuite_custom_atlas_prefix <- get_brainsute_custom_volume_atlas_prefix(brainsuite_custom_atlas_prefix)
  nii_atlas <- paste(brainsuite_custom_atlas_prefix, ".", bs_atlas_files_suffix$atlas_custom_suffix_tbm, sep="")
  nii_atlas_mask <- paste(brainsuite_custom_atlas_prefix, ".", bs_atlas_files_suffix$atlas_custom_mask_suffix_tbm, sep="")
  check_file_exists(nii_atlas, raise_error = TRUE)
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(list("nii_atlas" = nii_atlas, "nii_atlas_mask" = nii_atlas_mask))
}

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

get_custom_dbm_atlas_and_mask <- function(brainsuite_custom_atlas_prefix) {

  brainsuite_custom_atlas_prefix <- get_brainsute_custom_volume_atlas_prefix(brainsuite_custom_atlas_prefix)
  nii_atlas <- paste(brainsuite_custom_atlas_prefix, ".", bs_atlas_files_suffix$atlas_custom_suffix_dbm, sep = "")
  nii_atlas_mask <- paste(brainsuite_custom_atlas_prefix, ".", bs_atlas_files_suffix$atlas_custom_mask_suffix_dbm, sep = "")
  check_file_exists(nii_atlas, raise_error = TRUE)
  check_file_exists(nii_atlas_mask, raise_error = TRUE)
  return(list("nii_atlas" = nii_atlas, "nii_atlas_mask" = nii_atlas_mask))
}

read_demographics <- function(csvfile) {

  switch(tools::file_ext(csvfile),
         "tsv" = {demo <- read.table(file = csvfile, sep = "\t", header = T)},
         "csv" = {demo <- read.csv(csvfile)})
  # demo <- read.csv(csvfile)
  colnames(demo)[1] <- "subjID"
  return(demo)
}
