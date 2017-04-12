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
bss_load_roi_data <- function(subjects_dir, csv, roiid, roimeas, outdir=NULL) {

  if (!dir.exists(subjects_dir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", subjects_dir), call. = FALSE)
  }
  if (is.null(outdir)) {
    cat(sprintf('Output directory is not specified. Using %s to save outputs.\n', subjects_dir))
    outdir <- subjects_dir
  }
  else {
    dir.create(file.path(outdir), showWarnings = FALSE)
  }

  demographics <- read.csv(csv)
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

  roi_data <- read_roi_data_for_all_subjects(demographics$File_roi, roiid, roimeas)
  roi_columnname <- paste("ROI_", as.character(roiid), sep = "")
  demographics[[roi_columnname]] <- roi_data
  out_csv <- file.path(outdir, basename(csv))
  write.csv(demographics, out_csv, row.names = FALSE)

  bss_data <- list(df=demographics, outdir=outdir, roiid=roiid, roimeas=roimeas, csv=out_csv)

  # Finally also include the command to load the data
  bss_data$load_data_command <- sprintf("bss_data <- bss_load_roi_data('%s', '%s', %d, '%s', outdir = '%s') ", subjects_dir, csv, roiid, roimeas, outdir)

  return(bss_data)
}



