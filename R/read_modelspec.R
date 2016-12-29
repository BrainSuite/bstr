read_modelspec <- function(modelspec) {

  mspec <- ini::read.ini(modelspec)

  if (is.null(mspec$subject)) {
    stop(sprintf('The modelspec file %s does not contain a [subject] section.', modelspec), call. = FALSE)
  }
  if (is.null(mspec$stats)) {
    stop(sprintf('The modelspec file %s does not contain a [stats] section.', modelspec), call. = FALSE)
  }

  # Check for existence of [subject] fields
  if (is.null(mspec$subject$subjdir))
    stop('Section subjdir under [subject] not found. It should point to the top level directory that contains individual subjects.', call. = FALSE)
  else
    mspec$subjdir <- mspec$subject$subjdir

  if (is.null(mspec$subject$demographics))
    stop('Section demographics under [subject] not found. It should point to the csv/xls demographics file.', call. = FALSE)
  else
    mspec$csv <- mspec$subject$demographics

  mspec$smooth <- mspec$subject$smooth
  if (!is.null(mspec$smooth)) {
    mspec$smooth <- as.numeric(mspec$smooth)
    message(sprintf('Using smoothing level %2.1f.', mspec$smooth))
  }
  else
    message('Using no smoothing.')


  # Check for existence of [stats] fields
  if (is.null(mspec$stats$type))
    stop('Section type under [stats] not found. It should be either cbm, tbm or croi.', call. = FALSE)
  else {
    # Check if type is cbm, tbm or roi
    if (identical(mspec$stats$type, 'tbm') || identical(mspec$stats$type, 'cbm') || identical(mspec$stats$type, 'croi')) {
      mspec$type <- mspec$stats$type
    } else
      stop('Section type under [stats] should be either cbm, tbm or croi.', call. = FALSE)
  }

  mspec$main_effect <- mspec$stats$main_effect
  mspec$covariates <- mspec$stats$covariates
  mspec$corr_var <- mspec$stats$corr_var

  # Either main_effect and covariates or corr_var must be present
  if ( is.null(mspec$main_effect) && is.null(mspec$covariates)) {
    if (is.null(mspec$corr_var)) {
      stop('Either the main_effect and covariates OR corr_var must be present.')
    }
  }
  if (is.null(mspec$main_effect)) {
    if (is.null(mspec$covariates)) {
      stop('If main_effect is specified, covariates must also be specified.')
    }
  }
  if (is.null(mspec$covariates)) {
    if (is.null(mspec$main_effect)) {
      stop('If covariates is specified, main_effect must also be specified.')
    }
  }

  # TODO: Check if main_effect, covariates, corr_var are present in the demographics files

  # Open the svreg.log file and get the atlas file name
  if ( identical(tools::file_ext(mspec$csv), 'csv') ) {
    demo <- read.csv(mspec$csv)
  }
  # The first column has to contain subject IDs which are same as subject directories
  first_subjid <- demo[[1]][1]
  # Get atlas names from log files.
  svreg_log_file <- file.path(mspec$subjdir, first_subjid, sprintf('%s.svreg.log', first_subjid))
  if ( !file.exists(svreg_log_file) ) {
    stop(sprintf('Subject %s does not contain the svreg log file %s.\nPlease check if this is a valid subject directory and if svreg was run on all subjects.', first_subjid, svreg_log_file), call. = FALSE)
  }

  bs_atlas_path <- get_brainsuite_atlas_path_from_logfile(svreg_log_file)
  if (identical(mspec$stats$type, 'cbm')) {
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
        stop('One or more atlas files are not specified in the modelspec file.', call. = FALSE)
      else {
        if (file.exists(mspec$subject$lh_atlas) && file.exists(mspec$subject$rh_atlas)) {
          message(sprintf('Using atlas files %s and %s for left and right hemispheres.', mspec$subject$lh_atlas, mspec$subject$rh_atlas))
          mspec$lh_surf_atlas <- mspec$subject$lh_atlas
          mspec$rh_surf_atlas <- mspec$subject$rh_atlas
        }
        else
          stop(sprintf('Atlas file %s or %s do not exist', mspec$subject$lh_atlas, mspec$subject$rh_atlas), call. = FALSE)
      }
    }

    # if ( is.null(mspec$subject$atlas) ){
    #   mspec$lh_surf_atlas <- file.path(dirname(bs_atlas_path), bs_file_formats$surf_atlas_left)
    #   mspec$rh_surf_atlas <- file.path(dirname(bs_atlas_path), bs_file_formats$surf_atlas_right)
    # } else {
    #   # Check if mspec$subject$atlas points to a valid path
    #   if (file.exists(mspec$subject$atlas)){
    #     mspec$lh_surf_atlas <- file.path(dirname(mspec$subject$atlas), bs_file_formats$surf_atlas_left)
    #     mspec$rh_surf_atlas <- file.path(dirname(mspec$subject$atlas), bs_file_formats$surf_atlas_right)
    #   } else {
    #     warning(sprintf('The specified atlas file atlas=%s does not exist. Will use the atlas from the svreg.log file', mspec$subject$atlas), call. = FALSE)
    #     mspec$lh_surf_atlas <- file.path(dirname(bs_atlas_path), bs_file_formats$surf_atlas_left)
    #     mspec$rh_surf_atlas <- file.path(dirname(bs_atlas_path), bs_file_formats$surf_atlas_right)
    #   }
    # }
    return(mspec)
  }

  if (identical(mspec$stats$type, 'tbm')) {
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
  if (identical(mspec$stats$type, 'croi')) {
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
}

check_files_exist <- function(filename, raise=FALSE) {
  if ( file.exists(filename) )
    return(TRUE)
  else {
    if (raise)
      stop(sprintf('File %s does not exist', filename))
    else
      return(FALSE)
  }
}

