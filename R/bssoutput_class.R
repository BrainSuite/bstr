#' Bss output class
#' Defines an S4 class for statistical outputs
#' @export

check_files <- function(object){
  if (!dir.exists(object@subjdir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", object@subjdir), call. = FALSE)
  }

  if (!file.exists(object@csv)) {
    stop(sprintf("Demographics csv file %s does not exist.\n", object@csv), call. = FALSE)
  }
}

BssOutput <- setClass(
  "BssOutput",
  slots = list(
    outdir = "character"
  ),
  validity = check_files
)

setMethod("initialize", valueClass = "BssOutput", signature = "BssOutput", function(.Object, outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir)
    .Object@outdir <- outdir
  }
  else {
    .Object@outdir <- outdir
    message(sprintf("The output directory %s already exists. Will overwrite it's contents.", outdir))
  }
  return(.Object)
})

#' @export
setGeneric("save_out", valueClass = "BssOutput", function(bss_out, bss_data, bss_model) {
  standardGeneric("save_out")
})

#' @export
setMethod("save_out", valueClass = "BssOutput", signature = "BssOutput", function(bss_out, bss_data, bss_model) {
  print("in save")
})

BssCBMOutput <- setClass(
  "BssCBMOutput",
  contains = "BssOutput"
)

BssTBMOutput <- setClass(
  "BssTBMOutput",
  contains = "BssOutput"
)

BssROIOutput <- setClass(
  "BssROIOutput",
  contains = "BssOutput"
)

setMethod("save_out", valueClass = "BssCBMOutput", signature = "BssCBMOutput", function(bss_out, bss_data, bss_model) {

  s1 <- bss_data@atlas_surface
  log_pvalues <- log10_transform(bss_model@pvalues)
  s1$attributes <- log_pvalues
  bss_cmap <- new("BssColormap", "log_pvalues", "RdYlBu", log_pvalues)
  s1$vColor <- bss_cmap@rgbcolors
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  outprefix <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.dfs', sep = '')
  writedfs(file.path(bss_out@outdir, outprefix), s1)
  cbar_filename <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  save_colorbar(file.path(bss_out@outdir,cbar_filename), bss_cmap@lut, bss_cmap@vmin, bss_cmap@vmax, "p-values")

  log_pvalues_adjusted <- log10_transform(sign(bss_model@pvalues) * p.adjust(abs(bss_model@pvalues),
                                                                             method = 'BY'))
  s1$attributes <- log_pvalues_adjusted
  bss_cmap <- new("BssColormap", "log_pvalues_adjusted", "RdYlBu", log_pvalues_adjusted)
  s1$vColor <- bss_cmap@rgbcolors
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  outprefix <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.dfs', sep = '')
  writedfs(file.path(bss_out@outdir, outprefix), s1)
  cbar_filename <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  save_colorbar(file.path(bss_out@outdir,cbar_filename), bss_cmap@lut, bss_cmap@vmin, bss_cmap@vmax, "p-values")

  bss_model@tvalues[abs(log_pvalues) <= -1*log10(0.05)] <- 0
  bss_cmap <- new("BssColormap", "tvalues", "RdYlBu", bss_model@tvalues)
  s1$attributes <- bss_model@tvalues
  s1$vColor <- bss_cmap@rgbcolors
  #s1$vColor <- get_colors(bss_cmap)
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  outprefix <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.dfs', sep = '')
  writedfs(file.path(bss_out@outdir, outprefix), s1)
  cbar_filename <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  save_colorbar(file.path(bss_out@outdir,cbar_filename), bss_cmap@lut, bss_cmap@vmin, bss_cmap@vmax, "t-values")

  # Copy modelspec file to the output directory
  file.copy(bss_model@mspec_file, bss_out@outdir)
  return(bss_out)
  }
)

setMethod("save_out", valueClass = "BssTBMOutput", signature = "BssTBMOutput", function(bss_out, bss_data, bss_model) {

  log_pvalues <- rep(1, length(bss_data@atlas_image))
  log_pvalues[bss_data@mask_idx] <- log10_transform(bss_model@pvalues)
  dim(log_pvalues) <- dim(bss_data@atlas_image)
  bss_cmap <- new("BssColormap", "log_pvalues", "RdYlBu", as.numeric(log_pvalues))

  outprefix <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.nii.gz', sep = '')
  RNifti::writeNifti(log_pvalues, file.path(bss_out@outdir, outprefix), template = bss_data@atlas_image)
  cbar_filename <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  save_colorbar(file.path(bss_out@outdir,cbar_filename), bss_cmap@lut, bss_cmap@vmin, bss_cmap@vmax, "p-values")
  # save_BrainSuiteLUT
  lut_fileprefix <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.lut', sep = '')
  save_BrainSuiteLUT(file.path(bss_out@outdir, lut_fileprefix), bss_cmap@lut)

  log_pvalues_adjusted <- log10_transform(sign(bss_model@pvalues) * p.adjust(abs(bss_model@pvalues), method = 'BY'))
  bss_cmap <- new("BssColormap", "log_pvalues_adjusted", "RdYlBu", as.numeric(log_pvalues_adjusted))
  outprefix <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.nii.gz', sep = '')
  RNifti::writeNifti(as.matrix(log_pvalues_adjusted), file.path(bss_out@outdir, outprefix), template = bss_data@atlas_image)
  cbar_filename <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  save_colorbar(file.path(bss_out@outdir,cbar_filename), bss_cmap@lut, bss_cmap@vmin, bss_cmap@vmax, "p-values")
  # save_BrainSuiteLUT
  lut_fileprefix <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.lut', sep = '')
  save_BrainSuiteLUT(file.path(bss_out@outdir, lut_fileprefix), bss_cmap@lut)

  # TODO: Save t-values
  # browser()
  # bss_model@tvalues[abs(as.numeric(log_pvalues)) <= -1*log10(0.05)] <- 0
  # bss_cmap <- new("BssColormap", "tvalues", "RdYlBu", as.numeric(bss_model@tvalues))
  # outprefix <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
  #   basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.nii.gz', sep = '')
  # RNifti::writeNifti(as.matrix(log_pvalues_adjusted), file.path(bss_out@outdir, outprefix), template = bss_data@atlas_image)
  # cbar_filename <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
  #   basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  # save_colorbar(file.path(bss_out@outdir,cbar_filename), bss_cmap@lut, bss_cmap@vmin, bss_cmap@vmax, "t-values")
  # # save_BrainSuiteLUT
  # lut_fileprefix <- paste(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(
  #   basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.lut', sep = '')
  # save_BrainSuiteLUT(file.path(bss_out@outdir, lut_fileprefix), bss_cmap@lut)



    # Copy modelspec file to the output directory
  file.copy(bss_model@mspec_file, bss_out@outdir)
  return(bss_out)
  }
)

setMethod("save_out", valueClass = "BssROIOutput", signature = "BssROIOutput", function(bss_out, bss_data, bss_model) {

  # Get the absolute path of outdir
  outdir <- tools::file_path_as_absolute(bss_out@outdir)

  # Copy demographics csv file to output directory
  csvfilename <- file.path(bss_out@outdir, basename(bss_data@csv))
  write.csv(bss_data@demographics, csvfilename)

  nb_header <- "### BrainSuite ROI statistical analysis report"
#  nb_libraries <-"```{r librar_cmds, echo=FALSE}
#  library('bssr')
#  ```"

  nb_load_data <- sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, load_data}\n")
  # nb_load_data <- paste(nb_load_data, bss_data@load_data_command, sep = "")
  nb_load_data <- paste(nb_load_data, "\nDT::datatable(bss_data@demographics)\n", sep = "")
  nb_load_data <- paste(nb_load_data, "\n```\n", sep = "")


  nb_commands <- "```{r warning=FALSE, run_command}\n"

  for (i in bss_model@stats_commands) {
    nb_commands <- paste(nb_commands, i, "\n", sep = "")
  }
  nb_commands <- paste(nb_commands, "```\n", sep = "")
  nb_commands <- paste(nb_commands, sprintf("\n#### Main effect of %d %s on %s controlling for %s
                                            ", bss_data@roiid, bss_data@roimeas, bss_model@main_effect,
                                            bss_model@covariates ), sep = "")

  rmdfileconn<-file(file.path(outdir, "report.Rmd"))
  writeLines(c(nb_header, nb_load_data, nb_commands), rmdfileconn)
  close(rmdfileconn)

  # Render the markdown
  rmarkdown::render(file.path(outdir, "report.Rmd"), output_file=file.path(outdir, "report.html"), quiet = TRUE)

  # Copy modelspec file to the output directory
  file.copy(bss_model@mspec_file, bss_out@outdir)
  return(bss_out)
  }
)

#' @export
save_bss_out <- function(bss_data, bss_model, outdir="") {

  valid_types <- c("cbm", "tbm", "roi")
  if (! bss_data@data_type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  switch(bss_data@data_type,
         cbm = { bss_out <- new("BssCBMOutput", outdir)},
         tbm = { bss_out <- new("BssTBMOutput", outdir) },
         roi = { bss_out <- new("BssROIOutput", outdir) }
  )
  bss_out <- save_out(bss_out, bss_data, bss_model)
  invisible(bss_out)
}
