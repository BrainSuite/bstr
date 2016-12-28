#' List of file formats used in BrainSuite
#' @export

## TODO: use closures for this in the future
bs_file_formats <- list(
  surf_left = 'atlas.pvc-thickness_0-6mm.left.mid.cortex.dfs',
  surf_right = 'atlas.pvc-thickness_0-6mm.right.mid.cortex.dfs',
  surf_left_smooth = 'atlas.pvc-thickness_0-6mm.*mm.left.mid.cortex.dfs',
  surf_right_smooth = 'atlas.pvc-thickness_0-6mm.*mm.right.mid.cortex.dfs',
  jacdet = '*.svreg.inv.map.jacdet.*.nii.gz',
  roi_txt = '.roiwise.stats.txt',
  svreg_log = '*.svreg.log',
  surf_atlas_left = 'mri.left.mid.cortex.dfs',
  surf_atlas_right = 'mri.right.mid.cortex.dfs',
  nii_atlas = 'mri.bfc.nii.gz',
  nii_file = '%s.svreg.inv.map.jacdet.nii.gz',
  nii_maskfile = 'mri.cerebrum.mask.nii.gz',
  nii_file_smooth = '%s.svreg.inv.map.jacdet.smooth%2.1fmm.nii.gz'
)

analysis_type_list <- list(
  cbm = 'cbm',
  tbm = 'tbm',
  roi = 'roi',
  dbm = 'dbm',
  nca = 'nca'
)

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

get_cbm_file_list <- function(bss_data, hemi, smooth = NULL) {

  if (!is.null(smooth)) {
    if (identical(hemi, 'left')) {
      # Replace '*' in surf_left_smooth by smoothing value
      surf_left_smooth <- gsub("\\*", sprintf('smooth%2.1f', smooth), bs_file_formats$surf_left_smooth )
      cbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, surf_left_smooth)
    }
    else {
      surf_right_smooth <- gsub("\\*", sprintf('smooth%2.1f', smooth), bs_file_formats$surf_right_smooth )
      cbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, surf_right_smooth)
    }
  }
  else {
    if (identical(hemi, 'left'))
      cbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, bs_file_formats$surf_left)
    else
      cbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, bs_file_formats$surf_right)
  }
  # Check if all subjects have dfs files
  if ( !all(file.exists(cbm_filelist)) ) {
    message('Following subjects have missing dfs files')
    print(cbm_filelist[which(!file.exists(cbm_filelist))], row.names = FALSE)
    stop('\nCheck if svreg was run succesfully on all the subjects. Also check the smoothing level (smooth= under [subject]).\nIt is possible that surface files at the specified smoothing level do not exist.',
         call. = FALSE)
  }
  return(cbm_filelist)
}

get_tbm_file_list <- function(bss_data, smooth = NULL) {

  if (is.null(smooth)) {
    tbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, bs_file_formats$nii_file)
  }
  else {
    nii_smooth <- sprintf(bs_file_formats$nii_file_smooth, bss_data@demographics$subjID, smooth)
    tbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, nii_smooth)
    # tbm_filelist <- file.path(bss_data@subjdir, bss_data@demographics$subjID, bs_file_formats$nii_file_smooth)
  }
  # Check if all subjects have nii.gz files
  if ( !all(file.exists(tbm_filelist)) ) {
    message('Following subjects have missing nii.gz files')
    print(tbm_filelist[which(!file.exists(tbm_filelist))], row.names = FALSE)
    stop('\nCheck if svreg was run succesfully and if jacdet* files exist for all the subjects.\nAlso check if smoothing was performed.', call. = FALSE)
  }
  return(tbm_filelist)
}

get_brainsuite_atlas_path_from_logfile <- function(logfile) {
  fid = file(logfile, "rt")
  log_lines <- readLines(fid, n=2)
  bs_atlas_path <- unlist(strsplit(log_lines[2], ' ', fixed = TRUE))[3]
  close(fid)
  return(bs_atlas_path)

}

