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

#' S4 class for storing data for statistical analysis
#' @slot data_array matrix containing data of dimensions (N x T), where N = number of subjects and T = number of vertices/voxels.
#' @slot data_array_lh matrix containing data for left hemisphere.
#' @slot data_array_rh matrix containing data for right hemisphere.
#' @slot analysis_type character string denoting the type of analysis. Valid types are "sba", "tbm" or "roi".
#' @slot data_type character string denoting the type of data. Valid types are "surface" or "nifti_image".
#' @slot demographics data.frame containing the demographic information. Usually loaded from a csv file.
#' @slot subjdir character string for subject directory.
#' @slot csv filename of a comma separated (csv) file containing the subject demographic information.
#' @slot smooth numeric value used to smooth the data.
#' @slot measure character string denoting the type of measure used.
#' @slot filelist list of files belonging to N subjects.
#' @slot load_data_command character string for the command used to load the data
#' @slot hemi character string to specifiy the hemisphere. Valid values are "left" or "right" or "both"
#' @slot nvertices_lh numeric value containing the number of vertices for the left cortical surface atlas
#' @slot nvertices_rh numeric value containing the number of vertices for the right cortical surface atlas
#'
#' @export
BstrData <- setClass(
  "BstrData",
  slots = list(
    data_array = "matrix",
    data_array_lh = "matrix",
    data_array_rh = "matrix",
    analysis_type = "character",
    data_type = "character",
    demographics = "data.frame",
    subjdir = "character",
    csv = "character",
    smooth = "numeric",
    measure = "character",
    filelist = "character",
    load_data_command = "character",
    hemi = "character",
    nvertices_lh = "numeric",
    nvertices_rh = "numeric"
  )
)

BstrROIData <- setClass(
  "BstrROIData",
  slots = list(roiids = "numeric",
               roimeas = "character"),
  contains = "BstrData"
)

BstrSBAData <- setClass(
  "BstrSBAData",
  slots = list(atlas_filename = "character",
               atlas_filename_lh = "character",
               atlas_filename_rh = "character",
               atlas_surface = 'list',
               atlas_surface_lh = 'list',
               atlas_surface_rh = 'list'),
  contains = "BstrData"
)

setOldClass("niftiImage") #Declare niftiImage (RNifti) so the @slot atlas_image can be defined

BstrTBMData <- setClass(
  "BstrTBMData",
  slots = list(atlas_filename = "character",
               atlas_image = 'niftiImage',
               maskfile = 'character',
               mask_idx = 'vector'),
  contains = "BstrData"
)

BstrDBAData <- setClass(
  "BstrDBAData",
  slots = list(atlas_filename = "character",
               atlas_image = 'niftiImage',
               maskfile = 'character',
               mask_idx = 'vector'),
  contains = "BstrData"
)

setMethod("initialize", valueClass = "BstrData", signature = "BstrData", function(.Object, subjdir, csv, exclude_col) {

  check_file_exists(subjdir, raise_error = TRUE)
  check_file_exists(csv, raise_error = TRUE)
  .Object@subjdir = subjdir
  .Object@csv <- csv
  .Object@demographics <- read_demographics(csv, exclude_col = exclude_col)
  return(.Object)
})

#' A generic function to load data for statistical analysis.
#' @param bstr_data object of type [BstrData()]
#' @param atlas_filename path name to the atlas
#' @param maskfile path name to the mask file
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right" or "both".
#' @param measure character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".
#' @param smooth numeric value denoting the smoothing level.
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#' @param roiids numeric label identifier for the region of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. Should either be "gmthickness", "gmvolume", or "wmvolume".
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @details
#' For the most part, the user will never have to call this function directly.
#' Instead the user should call [load_bstr_data()].
#' @seealso [load_bstr_data()]
#'
#' @export
setGeneric("load_data", valueClass = "BstrData", function(bstr_data, atlas_filename = NULL, maskfile = NULL, hemi = "left", measure = "", smooth = 0.0, eddy = TRUE, roiids = NULL, roimeas = NULL, exclude_col) {
  standardGeneric("load_data")
})

setGeneric("load_demographics", valueClass = "BstrData", function(object) {
  standardGeneric("load_demographics")
})

#' @rdname load_data
setMethod("load_data", signature = "BstrData", function(bstr_data, roiids = NULL, roimeas = NULL, exclude_col) {
  return(bstr_data)
})

#' @rdname load_data
setMethod("load_data", signature = "BstrSBAData", function(bstr_data, atlas_filename, hemi, smooth) {

  bstr_data@atlas_filename <- atlas_filename
  bstr_data@atlas_surface <- readdfs(atlas_filename)
  sba_filelist <- get_sba_file_list(bstr_data, hemi, smooth)
  attrib_siz <- bstr_data@atlas_surface$hdr$nVertices
  bstr_data@data_array <- read_dfs_attributes_for_all_subjects(sba_filelist, attrib_siz)
  bstr_data@filelist <- sba_filelist
  bstr_data@analysis_type <- "sba"
  bstr_data@smooth <- smooth
  bstr_data@data_type <- bs_data_types$surface
  return(bstr_data)
})

#' @rdname load_data
setMethod("load_data", signature = "BstrTBMData", function(bstr_data, atlas_filename, maskfile = NULL, smooth) {

  bstr_data@atlas_filename <- atlas_filename
  bstr_data@atlas_image <- RNifti::readNifti(atlas_filename)
  bstr_data@filelist <- get_tbm_file_list(bstr_data, smooth)
  bstr_data@smooth <- smooth
  bstr_data@measure <- ""
  attrib_siz <- length(bstr_data@atlas_image)
  if ( !is.null(maskfile) ) {
    bstr_data@maskfile <- maskfile
    mask_image <- as.vector(RNifti::readNifti(maskfile))
    if ( length(mask_image) != attrib_siz) {
      stop(sprintf('Dimensions of atlas file %s and maskfile %s do not match', atlas_filename, maskfile), call. = FALSE)
    }
    bstr_data@mask_idx <- which(mask_image > 0)
  }
  else
    bstr_data@mask_idx = 1:attrib_siz

  first_subject_file <- as.vector(RNifti::readNifti(bstr_data@filelist[1]))
  # Check if the dimensions of first subject file and atlas match (the dimensions of atlas and mask are already checked above)
  if ( length(first_subject_file) != attrib_siz) {
    stop(sprintf('Dimensions of the atlas file %s and subject %s do not match. Check if you are using the correct atlas', atlas_filename, bstr_data@filelist[1]), call. = FALSE)
  }

  bstr_data@data_array <- read_nii_images_for_all_subjects(bstr_data@filelist, attrib_siz, bstr_data@mask_idx)
  bstr_data@analysis_type <- "tbm"
  bstr_data@data_type <- bs_data_types$nifti_image
  bstr_data@hemi <- "NA"
  return(bstr_data)
})

#' @rdname load_data
setMethod("load_data", signature = "BstrDBAData", function(bstr_data, atlas_filename, maskfile = NULL, measure, smooth, eddy) {

  bstr_data@atlas_filename <- atlas_filename
  bstr_data@atlas_image <- RNifti::readNifti(atlas_filename)
  bstr_data@filelist <- get_dba_file_list(bstr_data, measure, smooth, eddy)
  bstr_data@smooth <- smooth
  bstr_data@measure <- measure
  attrib_siz <- length(bstr_data@atlas_image)
  if ( !is.null(maskfile) ) {
    bstr_data@maskfile <- maskfile
    mask_image <- as.vector(RNifti::readNifti(maskfile))
    if ( length(mask_image) != attrib_siz) {
      stop(sprintf('Dimensions of atlas file %s and maskfile %s do not match', atlas_filename, maskfile), call. = FALSE)
    }
    bstr_data@mask_idx <- which(mask_image > 0)
  }
  else
    bstr_data@mask_idx = 1:attrib_siz

  first_subject_file <- as.vector(RNifti::readNifti(bstr_data@filelist[1]))
  # Check if the dimensions of first subject file and atlas match (the dimensions of atlas and mask are already checked above)
  if ( length(first_subject_file) != attrib_siz) {
    stop(sprintf('Dimensions of the atlas file %s and subject %s do not match', atlas_filename, bstr_data@filelist[1]), call. = FALSE)
  }

  bstr_data@data_array <- read_nii_images_for_all_subjects(bstr_data@filelist, attrib_siz, bstr_data@mask_idx)
  bstr_data@analysis_type <- "dba"
  bstr_data@data_type <- bs_data_types$nifti_image
  bstr_data@hemi <- "NA"
  return(bstr_data)
})

#' @rdname load_data
setMethod("load_data", signature = "BstrROIData", function(bstr_data, roiids = NULL, roimeas = NULL, exclude_col) {

  bstr_data@analysis_type <- "roi"
  bstr_data@roiids <- roiids
  bstr_data@roimeas <- roimeas
  all_subjects <- bstr_load_roi_data(subjects_dir = bstr_data@subjdir,
                    csv = bstr_data@csv,
                    roiids = bstr_data@roiids,
                    roimeas = bstr_data@roimeas,
                    exclude_col)

  bstr_data@demographics <- as.data.frame(all_subjects[[1]])
  bstr_data@data_array <- matrix(nrow=nrow(bstr_data@demographics),ncol = length(bstr_data@roiids))
  for (col in 1:length(bstr_data@roiids)){
    current_col <- which(colnames(bstr_data@demographics) == paste0(bstr:::get_roi_tag(label_desc_df = bstr:::read_label_desc(),roiid=bstr_data@roiids[col])[[1]],"(",bstr_data@roiids[col],")"))
    bstr_data@data_array[,col] <- bstr_data@demographics[,current_col]
  }
  bstr_data@load_data_command <- sprintf("bstr_data <- load_bstr_data(type= 'roi',subjdir = '%s',csv= '%s',roiids= c( %s), roimeas= '%s', exclude_col='%s')",
                                        bstr_data@subjdir, bstr_data@csv, paste(bstr_data@roiids,collapse = ", "), bstr_data@roimeas, exclude_col)


  return(bstr_data)

})


setMethod ("load_demographics", "BstrData", function(object) {
  object@demographics <- read_demographics(object@csv)
  return(object)
})



# signature(c(object = "BstrData", csv = "character"))

# setMethod("load_data", signature(object = "BstrROIData", subjects_dir = "character", csv = "character", roiids = "numeric",
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
#' should have the subject identifiers. Subject identifiers can be alphanumeric
#' and should be exactly equal to the individual subject directory names.
#'
#' @param type character string denoting type of analysis. Should be sba, tbm, or roi.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param smooth numeric value denoting the smoothing level.
#' @param roiids numeric label identifiers for the regions of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. Should either be "gmthickness", "gmvolume", or "wmvolume".
#' @param measure character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#' @param maskfile filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @examples
#' \dontrun{
#' my_sba_data <- load_bstr_data(type="sba", subjdir = "/path/to/my/subjectdirectory",
#' csv = "/path/to/my/demographics.csv", hemi = "left", smooth = 2.5)
#'
#' my_roi_data <- load_bstr_data(type="roi", subjdir = "/path/to/my/subjectdirectory",
#' csv="/path/to/my/demographics.csv", roiids=501, roimeas="gmthickness")
#' }
#'
#' @export
load_bstr_data <- function(type="sba", subjdir="", csv="", hemi="left",
                          smooth=0.0, roiids=0, roimeas="gmthickness", measure="", atlas="", maskfile = "", eddy=TRUE, exclude_col = "") {

  atlas <- path.expand(atlas)
  maskfile <- path.expand(maskfile)
  subjdir <- path.expand(subjdir)

  valid_types <- c("sba", "tbm", "roi","dba","nca")
  if (! type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  switch(type,
         sba = { bstr_data <- load_sba_data_both_hemi(subjdir=subjdir, csv=csv, hemi=hemi, smooth = smooth, atlas=atlas, exclude_col=exclude_col) },
         tbm = { bstr_data <- load_tbm_data(subjdir=subjdir, csv=csv, smooth=smooth, atlas=atlas, maskfile=maskfile, exclude_col=exclude_col) },
         dba = { bstr_data <- load_dba_data(subjdir=subjdir, csv=csv, measure=measure, smooth=smooth, atlas=atlas, maskfile=maskfile, eddy=eddy, exclude_col=exclude_col) },
         roi = { bstr_data <- load_roi_data(subjdir, csv, roiids, roimeas, exclude_col=exclude_col) }
  )
  return(bstr_data)
}

#' Load data for statistical analysis from a filelist
#'
#' Loading data is usually the first step before running any statistical analysis.
#' Prior to using this function, BrainSuite and svreg should be run on all subjects.
#' If required, smoothing should be performed on cortical surface or volumetric image based measures.
#' Unlike [load_bstr_data()], this function loads data from a csv that contains a column for filelist
#' A csv file containing subject demographic information should exist. The first column of this csv file
#' should have the subject identifiers. Subject identifiers can be alphanumeric
#' and should be exactly equal to the individual subject directory names.
#' This csv file should also contain a column that contains a full path to the data file to be loaded.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param type character string denoting type of analysis. Should be sba, tbm, or roi.
#' @param file_col character string for the full file path.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param maskfile optional filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @examples
#' \dontrun{
#' my_data <- load_bstr_data_from_filelist(csv = "/path/to/my/demographics.csv",
#' type="sba", file_col = "COL_NAME", atlast = "/path/to/atlas",
#' maskfile = "/path/to/maskfile")
#' }
#'
#' @export
load_bstr_data_from_filelist <- function(csv="", subjdir="", hemi = "left", type="sba", file_col="", atlas="", maskfile = "") {

#  bstr_sba_data <- new("BstrSBAData", subjdir=subjdir, csv=csv, exclude_col="")

  switch(type,
         tbm = { bstr_data <- load_tbm_data_from_filelist(subjdir = subjdir, csv = csv, file_col = file_col, atlas = atlas, maskfile = maskfile) },
         sba = { bstr_data <- load_sba_data_from_filelist(subjdir=subjdir, csv=csv, hemi = hemi, file_col = file_col, atlas=atlas) }
  )
  return(bstr_data)
}

#' Load cortical surface data for statistical analysis.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param smooth numeric value denoting the smoothing level.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
load_sba_data <- function(subjdir="", csv="", hemi="left", smooth=0.0, atlas="", exclude_col) {

  bstr_sba_data <- new("BstrSBAData", subjdir, csv, exclude_col)
  if (atlas == "") {
    brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv, exclude_col))
    sba_surf_atlas <- get_sba_atlas(brainsuite_atlas_id, hemi)
  }
  else
    sba_surf_atlas <- get_custom_sba_atlas_and_mask(atlas)

  bstr_sba_data <- load_data(bstr_sba_data, atlas_filename = sba_surf_atlas, hemi = hemi, smooth=smooth)
  bstr_sba_data@data_type <- bs_data_types$surface
  bstr_sba_data@hemi <- hemi
  return(bstr_sba_data)
}

#' Load cortical surface data for statistical analysis.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param file_col character string for the full file path.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @export
load_sba_data_from_filelist <- function(subjdir="", csv="", file_col="", atlas="", hemi="left") {

  bstr_data <- new("BstrSBAData", subjdir, csv, exclude_col="")

  bstr_data@atlas_filename <- atlas
  bstr_data@atlas_surface <- readdfs(atlas)
  sba_filelist <- bstr_data@demographics[, file_col]
  attrib_siz <- bstr_data@atlas_surface$hdr$nVertices
  bstr_data@data_array <- read_dfs_attributes_for_all_subjects(sba_filelist, attrib_siz)
  bstr_data@filelist <- sba_filelist
  bstr_data@analysis_type <- "sba"
  bstr_data@smooth <- 0.0
  bstr_data@data_type <- bs_data_types$surface
  bstr_data@hemi <- hemi
  return(bstr_data)
}

load_sba_data_both_hemi <- function(subjdir="", csv="", hemi="left", smooth=0.0, atlas="", exclude_col) {
    if (hemi == "both") {
      bstr_data_lh <- load_sba_data(subjdir=subjdir, csv=csv, hemi="left", smooth = smooth, atlas=atlas, exclude_col=exclude_col)
      bstr_data_rh <- load_sba_data(subjdir=subjdir, csv=csv, hemi="right", smooth = smooth, atlas=atlas, exclude_col=exclude_col)
      bstr_sba_data_both_hemi <- bstr_data_lh
      bstr_sba_data_both_hemi@data_array <- cbind(bstr_data_lh@data_array, bstr_data_rh@data_array)
      bstr_sba_data_both_hemi@nvertices_lh <- dim(bstr_data_lh@data_array)[2]
      bstr_sba_data_both_hemi@nvertices_rh <- dim(bstr_data_rh@data_array)[2]
      bstr_sba_data_both_hemi@atlas_surface_lh <- bstr_data_lh@atlas_surface
      bstr_sba_data_both_hemi@atlas_surface_rh <- bstr_data_rh@atlas_surface
      bstr_sba_data_both_hemi@hemi <- "both"
      bstr_sba_data_both_hemi@atlas_filename_lh <- bstr_data_lh@atlas_filename
      bstr_sba_data_both_hemi@atlas_filename_rh <- bstr_data_rh@atlas_filename
      return(bstr_sba_data_both_hemi)
    }
    else {
      bstr_data <- load_sba_data(subjdir=subjdir, csv=csv, hemi=hemi, smooth = smooth, atlas=atlas, exclude_col=exclude_col)
      return(bstr_data)
    }
}


#' Load tensor-based morphometry data for statistical analysis.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param smooth numeric value denoting the smoothing level.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param maskfile filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
load_tbm_data <- function(subjdir="", csv="", smooth=0.0, atlas="", maskfile="", exclude_col) {

  bstr_tbm_data <- new("BstrTBMData", subjdir, csv, exclude_col)
  tbm_atlas_and_mask = list()

  tbm_atlas_and_mask <- check_tbm_atlas_and_mask(subjdir, csv, atlas, maskfile, exclude_col)

  bstr_tbm_data <- load_data(bstr_tbm_data, atlas_filename = tbm_atlas_and_mask$nii_atlas, maskfile = tbm_atlas_and_mask$nii_atlas_mask, smooth=smooth)
  bstr_tbm_data@data_type <- bs_data_types$nifti_image
  return(bstr_tbm_data)
}

#' Load TBM data from file list for statistical analysis.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param file_col character string for the full file path.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' @param maskfile filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#' @export
#'
load_tbm_data_from_filelist <- function(subjdir="", csv="", file_col="", atlas="", maskfile="", exclude_col = "") {

  bstr_data <- new("BstrTBMData", subjdir, csv, exclude_col="")
  tbm_atlas_and_mask <- check_tbm_atlas_and_mask(subjdir, csv, atlas, maskfile, exclude_col)

  bstr_data@atlas_filename <- tbm_atlas_and_mask$nii_atlas
  bstr_data@atlas_image <- RNifti::readNifti(tbm_atlas_and_mask$nii_atlas)
  bstr_data@filelist <- bstr_data@demographics[, file_col]
  bstr_data@smooth <- 0.0
  bstr_data@measure <- ""

  attrib_siz <- length(bstr_data@atlas_image)
  if ( !is.null(maskfile) ) {
    bstr_data@maskfile <- maskfile
    mask_image <- as.vector(RNifti::readNifti(tbm_atlas_and_mask$nii_atlas))
    if ( length(mask_image) != attrib_siz) {
      stop(sprintf('Dimensions of atlas file %s and maskfile %s do not match', bstr_data@atlas_filename, maskfile), call. = FALSE)
    }
    bstr_data@mask_idx <- which(mask_image > 0)
  }
  else
    bstr_data@mask_idx = 1:attrib_siz
  first_subject_file <- as.vector(RNifti::readNifti(bstr_data@filelist[1]))
  # Check if the dimensions of first subject file and atlas match (the dimensions of atlas and mask are already checked above)
  if ( length(first_subject_file) != attrib_siz) {
    stop(sprintf('Dimensions of the atlas file %s and subject %s do not match. Check if you are using the correct atlas', bstr_data@atlas_filename, bstr_data@filelist[1]), call. = FALSE)
  }
  bstr_data@data_array <- read_nii_images_for_all_subjects(bstr_data@filelist, attrib_siz, bstr_data@mask_idx)
  bstr_data@analysis_type <- "tbm"
  bstr_data@data_type <- bs_data_types$nifti_image
  bstr_data@hemi <- "NA"

  return(bstr_data)
}

#' Load diffusion data for statistical analysis.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param measure character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".
#' @param smooth numeric value denoting the smoothing level.
#' @param atlas character specifying the file path prefix (all characters in the file name upto the first ".") for the custom atlas. If empty, the atlas will be read from the svreg.log file in the subject directory.
#' Otherwise, for example, if the atlas for tensor based morphometry is located at /path/to/atlas/myatlas.mri.bfc.nii.gz, then specify atlas="/path/to/atlas/myatlas".
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#' @param maskfile filename of the mask for tbm or diffusion parameter analysis. The mask has to be in the atlas space.
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'

load_dba_data <- function(subjdir="", csv="", measure="", smooth=0.0, atlas="", eddy=TRUE, maskfile="", exclude_col) {

  bstr_dba_data <- new("BstrDBAData", subjdir, csv, exclude_col)
  dba_atlas_and_mask = list()

  dba_atlas_and_mask = list()
  if (maskfile == ""  && atlas == "") {
    brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv, exclude_col))
    dba_atlas_and_mask <- get_dba_atlas_and_mask(brainsuite_atlas_id)
  }
  else {
    if (maskfile != "" && atlas == "") {
      brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv, exclude_col))
      check_file_exists(maskfile, raise_error = TRUE)
      dba_atlas_and_mask$nii_atlas_mask <- maskfile
      dba_atlas_and_mask$nii_atlas <- get_dba_atlas(brainsuite_atlas_id)
    }
    else if (maskfile== "" && atlas != "") {
      brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv, exclude_col))
      check_file_exists(atlas, raise_error = TRUE)
      dba_atlas_and_mask$nii_atlas <- atlas
      dba_atlas_and_mask$nii_atlas_mask <- get_dba_mask(brainsuite_atlas_id)
    }
    else if (maskfile!= "" && atlas != "") {
      check_file_exists(maskfile, raise_error = TRUE)
      check_file_exists(atlas, raise_error = TRUE)
      dba_atlas_and_mask <- list("nii_atlas" = atlas, "nii_atlas_mask" = maskfile)
    }
  }


  bstr_dba_data <- load_data(bstr_dba_data, atlas_filename = dba_atlas_and_mask$nii_atlas, maskfile = dba_atlas_and_mask$nii_atlas_mask, measure=measure, smooth=smooth, eddy=eddy)
  bstr_dba_data@data_type <- bs_data_types$nifti_image
  return(bstr_dba_data)
}

#' Load ROI data for statistical analysis.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param roiids numeric label identifiers for the regions of interest (ROI) type analysis.
#' @param roimeas character string for the ROI measure. Should either be "gmthickness", "gmvolume", or "wmvolume".
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
load_roi_data <- function(subjdir="", csv="", roiids="", roimeas="", exclude_col) {
  bstr_roi_data <- new("BstrROIData", subjdir, csv, exclude_col)
  bstr_roi_data <- load_data(bstr_roi_data, roiids = roiids, roimeas = roimeas, exclude_col=exclude_col)
  return(bstr_roi_data)
}

# #' Check that subject directory and demographics csv files exist
# #' @param object object of type \code{BstrData}
# #'
# check_files <- function(object){
#   if (!dir.exists(object@subjdir)) {
#     stop(sprintf("Subjects directory %s does not exist.\n", object@subjdir), call. = FALSE)
#   }
#
#   if (!file.exists(object@csv)) {
#     stop(sprintf("Demographics csv file %s does not exist.\n", object@csv), call. = FALSE)
#   }
# }

#' Package data for reproducible statistical analysis.
#'
#' Takes same parameters as [load_bstr_data()] and copies the data to a new directory specified by outdir. You can repeatedly call this function to copy data of different types (tbm -- nii.gz, sba -- .dfs files etc.) to the same output directory.
#' Prior to using this function, BrainSuite and svreg should be run on all subjects.
#' If required, smoothing should be performed on cortical surface or volumetric image based measures.
#' A csv file containing subject demographic information should exist. The first column of this csv file
#' should have the subject identifiers. Subject identifiers can be alphanumeric
#' and should be exactly equal to the individual subject directory names.
#'
#' @param type character string denoting type of analysis. Should be sba, tbm, or roi.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param sbasmooth numeric value denoting the smoothing level for sba.
#' @param tbmsmooth numeric value denoting the smoothing level for tbm.
#' @param dbasmooth numeric value denoting the smoothing level for dba.
#' @param measure character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".
#' @param atlas path name to the atlas
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#' @param outdir output directory that will contain the copied data
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
#' @export
package_data <- function(type="sba", subjdir=NULL, csv="", hemi="left",
                         sbasmooth=0.0, tbmsmooth=0.0,
                         dbasmooth=0.0, measure="FA", atlas="", eddy=TRUE, outdir=NULL, exclude_col="") {

  valid_types <- c("sba", "tbm", "roi","dba","nca", "all")
  if (! type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  if (is.null(subjdir))
    stop(sprintf("Subject directory is not specified."), call. = FALSE)

  if (is.null(outdir))
    stop(sprintf("Output directory is not specified."), call. = FALSE)

  if (outdir == "") {
    stop("Output directory is an empty string.", call.=FALSE)
  }
  if (dir.exists(outdir))
    stop("Output directory exists. Please specify a new output directory that does not exist.", call.=FALSE)

  switch(type,
         sba = { copy_sba_data(subjdir=subjdir, csv=csv, hemi=hemi, smooth = sbasmooth, outdir=outdir, exclude_col=exclude_col) },
         tbm = { copy_tbm_data(subjdir=subjdir, csv=csv, smooth=tbmsmooth, atlas=atlas, outdir=outdir, exclude_col=exclude_col) },
         dba = { copy_dba_data(subjdir=subjdir, csv=csv, measure=measure,
                                           smooth=dbasmooth, atlas=atlas, eddy=eddy, outdir=outdir, exclude_col=exclude_col) },
         roi = { copy_roi_data(subjdir, csv, outdir=outdir, exclude_col=exclude_col) },
         all = {
           copy_sba_data(subjdir=subjdir, csv=csv, hemi="left", smooth = sbasmooth, outdir=outdir, exclude_col=exclude_col)
           copy_sba_data(subjdir=subjdir, csv=csv, hemi="right", smooth = sbasmooth, outdir=outdir, exclude_col=exclude_col)
           copy_tbm_data(subjdir=subjdir, csv=csv, smooth=tbmsmooth, atlas=atlas, outdir=outdir, exclude_col=exclude_col)
           copy_dba_data(subjdir=subjdir, csv=csv, measure=measure,
                         smooth=dbasmooth, atlas=atlas, eddy=eddy, outdir=outdir, exclude_col=exclude_col)
           copy_roi_data(subjdir, csv, outdir=outdir, exclude_col=exclude_col)
         }
  )
  # Copy the spreadsheet
  file.copy(csv, file.path(outdir, basename(csv)))
}

#' Package cortical surface data for reproducible statistical analysis.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param hemi chaaracter string denoting the brain hemisphere. Should either be "left" or "right".
#' @param smooth numeric value denoting the smoothing level.
#' @param outdir output directory that will contain the copied data
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
copy_sba_data <- function(subjdir="", csv="", hemi="left", smooth=0.0, outdir, exclude_col="") {

  bstr_data <- new("BstrSBAData", subjdir, csv, exclude_col )
  logfilenames <- get_brainsuite_logfilename_for_all_subjects(subjdir, csv, exclude_col)
  brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv, exclude_col))
  sba_atlas_filename <- get_sba_atlas(brainsuite_atlas_id, hemi)
  sba_filelist <- get_sba_file_list(bstr_data, hemi, smooth)
  src_filelist <- c(sba_filelist, logfilenames, sba_atlas_filename)

  # Create subdirectories for subject IDs in outdir
  dir.create(file.path(outdir), showWarnings = FALSE)
  Vectorize(dir.create)(file.path(outdir, bstr_data@demographics$subjID), showWarnings = FALSE)
  dest_filelist <- file.path(outdir, bstr_data@demographics$subjID, basename(sba_filelist))
  dest_filelist <- c(dest_filelist, file.path(outdir, bstr_data@demographics$subjID, basename(logfilenames)),
                     file.path(outdir, basename(sba_atlas_filename)))
  # Copy files
  file_copy(src_filelist, dest_filelist, messg = "Copying sba data")
}

#' Package tensor-based morphometry data for reproducible statistical analysis.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param smooth numeric value denoting the smoothing level.
#' @param atlas path name to the atlas
#' @param outdir output directory that will contain the copied data
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
copy_tbm_data <- function(subjdir="", csv="", smooth=0.0, atlas, outdir, exclude_col="") {

  bstr_data <- new("BstrTBMData", subjdir, csv, exclude_col)
  logfilenames <- get_brainsuite_logfilename_for_all_subjects(subjdir, csv, exclude_col)
  brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv, exclude_col))
  tbm_atlas_mask_filename <- get_tbm_atlas_and_mask(brainsuite_atlas_id)
  tbm_filelist <- get_tbm_file_list(bstr_data, smooth)
  src_filelist <- c(tbm_filelist, logfilenames, tbm_atlas_mask_filename$nii_atlas, tbm_atlas_mask_filename$nii_atlas_mask)

  # Create subdirectories for subject IDs in outdir
  dir.create(file.path(outdir), showWarnings = FALSE)
  Vectorize(dir.create)(file.path(outdir, bstr_data@demographics$subjID), showWarnings = FALSE)
  dest_filelist <- file.path(outdir, bstr_data@demographics$subjID, basename(tbm_filelist))
  dest_filelist <- c(dest_filelist, file.path(outdir, bstr_data@demographics$subjID, basename(logfilenames)),
                     file.path(outdir, basename(tbm_atlas_mask_filename$nii_atlas)),
                     file.path(outdir, basename(tbm_atlas_mask_filename$nii_atlas_mask))
                     )
  # Copy files
  message("Copying tbm data", appendLF = FALSE)
  file_copy(src_filelist, dest_filelist)
}

#' Package diffusion data for reproducible statistical analysis.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param measure character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".
#' @param atlas path name to the atlas
#' @param eddy boolean for specifying if the diffusion images were eddy-current corrected or not.
#' @param smooth numeric value denoting the smoothing level.
#' @param outdir output directory that will contain the copied data
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
copy_dba_data <- function(subjdir="", csv="", measure="FA", atlas="", eddy=TRUE, smooth=0.0, outdir, exclude_col="") {

  bstr_data <- new("BstrDBAData", subjdir, csv, exclude_col)
  logfilenames <- get_brainsuite_logfilename_for_all_subjects(subjdir, csv, exclude_col)
  brainsuite_atlas_id <- get_brainsuite_atlas_id_from_logfile(get_brainsuite_logfilename(subjdir, csv, exclude_col))
  dba_atlas_mask_filename <- get_dba_atlas_and_mask(brainsuite_atlas_id)

  message("Copying dba data", appendLF = FALSE)
  valid_diffusion_measures <- c('FA', 'MD', 'axial', 'radial', 'mADC', 'FRT_GFA')
  for (jj in valid_diffusion_measures) {
    dba_filelist <- get_dba_file_list(bstr_data, measure=jj, smooth = smooth, eddy = TRUE)
    src_filelist <- c(dba_filelist, logfilenames, dba_atlas_mask_filename$nii_atlas, dba_atlas_mask_filename$nii_atlas_mask)

    # Create subdirectories for subject IDs in outdir
    dir.create(file.path(outdir), showWarnings = FALSE)
    Vectorize(dir.create)(file.path(outdir, bstr_data@demographics$subjID), showWarnings = FALSE)
    dest_filelist <- file.path(outdir, bstr_data@demographics$subjID, basename(dba_filelist))
    dest_filelist <- c(dest_filelist, file.path(outdir, bstr_data@demographics$subjID, basename(logfilenames)),
                       file.path(outdir, basename(dba_atlas_mask_filename$nii_atlas)),
                       file.path(outdir, basename(dba_atlas_mask_filename$nii_atlas_mask))
    )
    # Copy files
    file_copy(src_filelist, dest_filelist)
  }
}

#' Package ROI data for reproducible statistical analysis.
#' @param subjdir subject directory containing BrainSuite processed data.
#' @param csv filename of a comma separated (csv) file containing the subject demographic information.
#' The first column of this csv file
#' should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
#' and should be exactly equal to the individual subject directory name.
#' @param outdir output directory that will contain the copied data
#' @param exclude_col character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.
#'
copy_roi_data <- function(subjdir="", csv="", outdir, exclude_col="") {

  bstr_data <- new("BstrROIData", subjdir, csv, exclude_col)
  roiwise_file_list <- get_roi_file_list(bstr_data)
  dest_filelist <- file.path(outdir, bstr_data@demographics$subjID, basename(roiwise_file_list))

  # Create subdirectories for subject IDs in outdir
  dir.create(file.path(outdir), showWarnings = FALSE)
  Vectorize(dir.create)(file.path(outdir, bstr_data@demographics$subjID), showWarnings = FALSE)

  # Copy files
  message("Copying ROI data", appendLF = FALSE)
  file_copy(roiwise_file_list, dest_filelist)
}

#' Copy files from the inputted source to the inputted destination
#' @param src_filelist list of source files
#' @param dest_filelist list of destination files
#' @param messg character string of displayed message
#' @param progress logical flag set TRUE to display progress bar in terminal
#'
file_copy <- function(src_filelist, dest_filelist, messg="Copying ", progress = TRUE) {

  if (length(src_filelist) != length(dest_filelist))
    stop(sprintf('The lengths of the source and the destination files do not match.'), call. = FALSE)

  message(messg, appendLF = FALSE)
  # If progress == TRUE, then display progressbar
  if (progress == TRUE) {
    pb <- txtProgressBar(max = length(src_filelist), style = 3)
    for (i in 1:length(src_filelist)) {
      file.copy(src_filelist[i], dest_filelist[i])
      setTxtProgressBar(pb, pb$getVal()+1)
    }
    close(pb)
  }
  else
    file.copy(src_filelist, dest_filelist)
}

