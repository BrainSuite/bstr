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

#' S4 class for storing data for statistical analysis
#' @slot data_array matrix containing data of dimensions (N x T), where N = number of subjects and T = number of vertices/voxels.
#' @slot data_array_lh matrix containing data for left hemisphere.
#' @slot data_array_rh atrix containing data for right hemisphere.
#' @slot analysis_type character string denoting the type of analysis. Valid types are "cbm", "tbm" or "roi".
#' @slot data_type character string denoting the type of data. Valid types are "surface" or "nifti_image".
#' @slot demographics data.frame containing the demographic information. Usually loaded from a csv file.
#' @slot subjdir character string for subject directory.
#' @slot csv filename of a comma separated (csv) file containing the subject demographic information.
#' @slot filelist list of files belonging to N subjects.
#'
#' @export
BssData <- setClass(
  "BssData",
  slots = list(
    data_array = "matrix",
    data_array_lh = "matrix",
    data_array_rh = "matrix",
    analysis_type = "character",
    data_type = "character",
    demographics = "data.frame",
    subjdir = "character",
    csv = "character",
    filelist = "character"
  )
)

BssROIData <- setClass(
  "BssROIData",
  slots = list(roiid = "numeric",
               roimeas = "character"),
  contains = "BssData"
)

BssCBMData <- setClass(
  "BssCBMData",
  slots = list(atlas_filename = "character",
               atlas_surface = 'list'),
  contains = "BssData"
)

setOldClass("niftiImage") #Declare niftiImage (RNifti) so the @slot atlas_image can be defined

BssTBMData <- setClass(
  "BssTBMData",
  slots = list(atlas_filename = "character",
               atlas_image = 'niftiImage',
               maskfile = 'character',
               mask_idx = 'vector'),
  contains = "BssData"
)

BssDBMData <- setClass(
  "BssDBMData",
  slots = list(atlas_filename = "character",
               atlas_image = 'niftiImage',
               maskfile = 'character',
               mask_idx = 'vector'),
  contains = "BssData"
)

setMethod("initialize", valueClass = "BssData", signature = "BssData", function(.Object, subjdir, csv) {

  check_file_exists(subjdir, raise_error = TRUE)
  check_file_exists(csv, raise_error = TRUE)
  .Object@subjdir = subjdir
  .Object@csv <- csv
  .Object@demographics <- read.csv(csv)
  return(.Object)
})

#' A generic function to load data for statistical analysis.
#' @param bss_data object of type \code{\link{BssData}}
#' @param atlas_filename path name to the atlas
#' @param maskfile path name to the mask file
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param smooth numeric value denoting the smoothing level.
#' @param roiid numeric label identifier for the region of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. Should either be "gmthickness", "gmvolume", or "wmvolume".
#' @details
#' For the most part, the user will never have to call this function directly.
#' Instead the user should call \code{\link{load_bss_data}}.
#' @seealso \code{\link{load_bss_data}}
#'
#' @export
setGeneric("load_data", valueClass = "BssData", function(bss_data, atlas_filename = NULL, maskfile = NULL, hemi = "left", smooth = 0.0, roiid = NULL, roimeas = NULL) {
  standardGeneric("load_data")
})

setGeneric("load_demographics", valueClass = "BssData", function(object) {
  standardGeneric("load_demographics")
})

#' @rdname load_data
setMethod("load_data", signature = "BssData", function(bss_data, roiid = NULL, roimeas = NULL) {
  return(bss_data)
})

#' @rdname load_data
setMethod("load_data", signature = "BssCBMData", function(bss_data, atlas_filename, hemi, smooth) {

  bss_data@atlas_filename <- atlas_filename
  bss_data@atlas_surface <- readdfs(atlas_filename)
  cbm_filelist <- get_cbm_file_list(bss_data, hemi, smooth)
  attrib_siz <- bss_data@atlas_surface$hdr$nVertices
  bss_data@data_array <- read_dfs_attributes_for_all_subjects(cbm_filelist, attrib_siz)
  bss_data@filelist <- cbm_filelist
  bss_data@analysis_type <- "cbm"
  bss_data@data_type <- bs_data_types$surface
  return(bss_data)
})

#' @rdname load_data
setMethod("load_data", signature = "BssTBMData", function(bss_data, atlas_filename, maskfile = NULL, smooth) {

  bss_data@atlas_filename <- atlas_filename
  bss_data@atlas_image <- RNifti::readNifti(atlas_filename)
  bss_data@filelist <- get_tbm_file_list(bss_data, smooth)
  attrib_siz <- length(bss_data@atlas_image)
  if ( !is.null(maskfile) ) {
    bss_data@maskfile <- maskfile
    mask_image <- as.vector(RNifti::readNifti(maskfile))
    if ( length(mask_image) != attrib_siz) {
      stop(sprintf('Dimensions of atlas file %s and maskfile %s do not match', atlas_filename, maskfile), call. = FALSE)
    }
    bss_data@mask_idx <- which(mask_image > 0)
  }
  else
    bss_data@mask_idx = 1:attrib_siz

  bss_data@data_array <- read_nii_images_for_all_subjects(bss_data@filelist, attrib_siz, bss_data@mask_idx)
  bss_data@analysis_type <- "tbm"
  bss_data@data_type <- bs_data_types$nifti_image
  return(bss_data)
})

#' @rdname load_data
setMethod("load_data", signature = "BssDBMData", function(bss_data, atlas_filename, maskfile = NULL, smooth) {
  
  bss_data@atlas_filename <- atlas_filename
  bss_data@atlas_image <- RNifti::readNifti(atlas_filename)
  bss_data@filelist <- get_dbm_file_list(bss_data, smooth)
  attrib_siz <- length(bss_data@atlas_image)
  if ( !is.null(maskfile) ) {
    bss_data@maskfile <- maskfile
    mask_image <- as.vector(RNifti::readNifti(maskfile))
    if ( length(mask_image) != attrib_siz) {
      stop(sprintf('Dimensions of atlas file %s and maskfile %s do not match', atlas_filename, maskfile), call. = FALSE)
    }
    bss_data@mask_idx <- which(mask_image > 0)
  }
  else
    bss_data@mask_idx = 1:attrib_siz
  
  bss_data@data_array <- read_nii_images_for_all_subjects(bss_data@filelist, attrib_siz, bss_data@mask_idx)
  bss_data@analysis_type <- "dbm"
  bss_data@data_type <- bs_data_types$nifti_image
  return(bss_data)
})

#' @rdname load_data
setMethod("load_data", signature = "BssROIData", function(bss_data, roiid = NULL, roimeas = NULL) {

  # if (is.null(outdir)) {
  #   cat(sprintf('Output directory is not specified. Using %s to save outputs.\n', subjects_dir))
  #   outdir <- subjects_dir
  # }
  # else {
  #   dir.create(file.path(outdir), showWarnings = FALSE)
  # }
  demographics <- bss_data@demographics
  if("File_roi" %in% colnames(demographics)){
    warning(sprintf("The file %s already contains a File_roi column.\nWill overwrite this column.\n", bss_data@csv), call.=FALSE)
  }

  demographics$subjID <- as.character(demographics$subjID)
  roiwise_file_list <- get_roi_file_list(bss_data)
  demographics$File_roi <- unlist(roiwise_file_list)

  roi_data <- read_roi_data_for_all_subjects(demographics$File_roi, roiid, roimeas)
  roi_columnname <- paste("ROI_", as.character(roiid), sep = "")
  demographics[[roi_columnname]] <- roi_data
  # out_csv <- file.path(outdir, basename(csv))
  # write.csv(demographics, out_csv, row.names = FALSE)
  # object@
  bss_data@analysis_type <- "roi"
  bss_data@roiid <- roiid
  bss_data@roimeas <- roimeas
  bss_data@demographics <- demographics

  # Finally also include the command to load the data
  # bss_data$load_data_command <- sprintf("bss_data <- bss_load_roi_data('%s', '%s', %d, '%s', outdir = '%s') ", subjects_dir, csv, roiid, roimeas, outdir)

  return(bss_data)

})


setMethod ("load_demographics", "BssData", function(object) {
  object@demographics <- read.csv(object@csv)
  return(object)
})



# signature(c(object = "BssData", csv = "character"))

# setMethod("load_data", signature(object = "BssROIData", subjects_dir = "character", csv = "character", roiid = "numeric",
#                                  roimeas = "character", outdir = "character"),
#           function(object, subjects_dir, csv, roiid, roimeas, outdir=NULL) {
#             print("hi")
#             }
#           )

#' Load data for statistical analysis.
#'
#' Loading data is usually the first step before running any statistical analysis.
#' Prior to using this function, BrainSuite and svreg should be run on all subjects.
#' If required, smoothing should be performed on cortical surface or volumetric image based measures.
#' A csv file containing subject demographic information should exist. The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#'
#' @param type character string denoting type of analysis. Should be cbm, tbm, or roi.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param smooth numeric value denoting the smoothing level.
#' @param roiid numeric label identifier for the region of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. Should either be "gmthickness", "gmvolume", or "wmvolume".
#' @examples
#' \dontrun{
#' my_cbm_data <- load_bss_data(type="cbm", subjdir = "/path/to/my/subjectdirectory",
#' csv = "/path/to/my/demographics.csv", hemi = "left", smooth = 2.5)
#'
#' my_roi_data <- load_bss_data(type="roi", subjdir = "/path/to/my/subjectdirectory",
#' csv="/path/to/my/demographics.csv", roiid=501, roimeas="gmthickness")
#' }
#'
#' @export
load_bss_data <- function(type="cbm", subjdir="", csv="", hemi="left", smooth=0.0, roiid=0, roimeas="gmthickness") {

  valid_types <- c("cbm", "tbm", "roi","dbm","nca")
  if (! type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  switch(type,
         cbm = { bss_data <- load_cbm_data(subjdir=subjdir, csv=csv, hemi=hemi, smooth = smooth) },
         tbm = { bss_data <- load_tbm_data(subjdir=subjdir, csv=csv, smooth=smooth) },
         dbm = { bss_data <- load_dbm_data(subjdir=subjdir, csv=csv, smooth=smooth) },
         roi = { bss_data <- load_roi_data(subjdir, csv, roiid, roimeas) }
  )
  return(bss_data)
}


load_cbm_data <- function(subjdir="", csv="", hemi="left", smooth=0.0) {

  bss_cbm_data <- new("BssCBMData", subjdir, csv)
  brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv))
  cbm_surf_atlas <- get_cbm_atlas(brainsuite_atlas_id, hemi)
  bss_cbm_data <- load_data(bss_cbm_data, atlas_filename = cbm_surf_atlas, hemi = hemi, smooth=smooth)
  bss_cbm_data@data_type <- bs_data_types$surface
  return(bss_cbm_data)
}

load_tbm_data <- function(subjdir="", csv="", smooth=0.0) {

  bss_tbm_data <- new("BssTBMData", subjdir, csv)
  brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv))
  tbm_atlas_and_mask <- get_tbm_atlas_and_mask(brainsuite_atlas_id)
  bss_tbm_data <- load_data(bss_tbm_data, atlas_filename = tbm_atlas_and_mask$nii_atlas, maskfile = tbm_atlas_and_mask$nii_atlas_mask, smooth=smooth)
  bss_tbm_data@data_type <- bs_data_types$nifti_image
  return(bss_tbm_data)
}

load_dbm_data <- function(subjdir="", csv="", smooth=0.0) {
  
  bss_dbm_data <- new("BssDBMData", subjdir, csv)
  brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv))
  tbm_atlas_and_mask <- get_tbm_atlas_and_mask(brainsuite_atlas_id)
  bss_dbm_data <- load_data(bss_dbm_data, atlas_filename = tbm_atlas_and_mask$nii_atlas, maskfile = tbm_atlas_and_mask$nii_atlas_mask, smooth=smooth)
  bss_dbm_data@data_type <- bs_data_types$nifti_image
  return(bss_dbm_data)
}


load_roi_data <- function(subjdir="", csv="", roiid="", roimeas="") {
  bss_roi_data <- new("BssROIData", subjdir, csv)
  bss_roi_data <- load_data(bss_roi_data, roiid = roiid, roimeas = roimeas)
  return(bss_roi_data)
}

check_files <- function(object){
  if (!dir.exists(object@subjdir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", object@subjdir), call. = FALSE)
  }

  if (!file.exists(object@csv)) {
    stop(sprintf("Demographics csv file %s does not exist.\n", object@csv), call. = FALSE)
  }
}
