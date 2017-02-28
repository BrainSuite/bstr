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

BssData <- setClass(
  "BssData",
  slots = list(
    data_array = "matrix",
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
setGeneric("load_data", valueClass = "BssData", function(bss_data, atlas_filename = NULL, maskfile = NULL, hemi = NULL, smooth = NULL, roiid = NULL, roimeas = NULL) {
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
  bss_data@data_type <- "roi"
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
