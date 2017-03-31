#' Load BrainSuite processed data for statistical analysis
#' Defines an S4 class for storing/loading data
#' @export

check_files <- function(object){
  if (!dir.exists(object@subjdir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", object@subjdir), call. = FALSE)
  }

  if (!file.exists(object@csv)) {
    stop(sprintf("Demographics csv file %s does not exist.\n", object@csv), call. = FALSE)
  }
}

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
  ),
  validity = check_files
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

setMethod("initialize", valueClass = "BssData", signature = "BssData", function(.Object, subjdir, csv) {
  if (dir.exists(subjdir)) {
    .Object@subjdir = subjdir
  }

  if (file.exists(csv)) {
    .Object@csv <- csv
    .Object@demographics <- read.csv(csv)
  }
  return(.Object)
})

#' @export
setGeneric("load_data", valueClass = "BssData", function(bss_data, atlas_filename = NULL, maskfile = NULL, hemi = "left", smooth = 0.0, roiid = NULL, roimeas = NULL) {
  standardGeneric("load_data")
})

setGeneric("load_demographics", valueClass = "BssData", function(object) {
  standardGeneric("load_demographics")
})

setMethod("load_data", signature = "BssData", function(bss_data, roiid = NULL, roimeas = NULL) {
  return(bss_data)
})

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

#' @export
load_bss_data <- function(type="cbm", subjdir="", csv="", hemi="left", smooth=0.0, roiid=0, roimeas="gmthickness") {

  valid_types <- c("cbm", "tbm", "roi")
  if (! type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  switch(type,
         cbm = { bss_data <- load_cbm_data(subjdir=subjdir, csv=csv, hemi=hemi, smooth = smooth) },
         tbm = { bss_data <- load_tbm_data(subjdir=subjdir, csv=csv, smooth=smooth) },
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

load_roi_data <- function(subjdir="", csv="", roiid="", roimeas="") {
  bss_roi_data <- new("BssROIData", subjdir, csv)
  bss_roi_data <- load_data(bss_roi_data, roiid = roiid, roimeas = roimeas)
  return(bss_roi_data)
}
