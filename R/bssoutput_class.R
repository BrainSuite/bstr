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

  log_pvalues <- log10_transform(bss_model@pvalues)
  outdir <- bss_out@outdir
  log_pvalues_adjusted <- log10_transform(sign(bss_model@pvalues) * p.adjust(abs(bss_model@pvalues),
                                                                             method = 'BY'))
  bss_model@tvalues[abs(log_pvalues) <= -1*log10(0.05)] <- 0

  switch(bss_model@model_type,
         bss_anova = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@main_effect, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@main_effect, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues_adjusted, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(bss_model@tvalues, bss_model@main_effect, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(bss_model@tvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
         },
         bss_lm = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@main_effect, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@main_effect, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues_adjusted, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(bss_model@tvalues, bss_model@main_effect, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(bss_model@tvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           },
         bss_corr = {
           bss_model@corr_values[abs(log_pvalues) <= -1*log10(0.05)] <- 0
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@corr_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@corr_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues_adjusted, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(bss_model@corr_values, bss_model@corr_var, "corr_values", bss_data, bss_model, outdir)
           save_bss_out_surface(bss_model@corr_values, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           },
         pairedttest = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@group_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@group_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues_adjusted, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(bss_model@tvalues, bss_model@group_var, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(bss_model@tvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
         },
         unpairedttest = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@group_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@group_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues_adjusted, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(bss_model@tvalues, bss_model@group_var, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(bss_model@tvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
         }
  )

  # Copy modelspec file to the output directory
  file.copy(bss_model@mspec_file, bss_out@outdir)
  return(bss_out)
  }
)

setMethod("save_out", valueClass = "BssTBMOutput", signature = "BssTBMOutput", function(bss_out, bss_data, bss_model) {

  log_pvalues <- rep(1, length(bss_data@atlas_image))
  log_pvalues[bss_data@mask_idx] <- log10_transform(bss_model@pvalues)
  dim(log_pvalues) <- dim(bss_data@atlas_image)

  log_pvalues_adjusted <- rep(1, length(bss_data@atlas_image))
  log_pvalues_adjusted[bss_data@mask_idx] <- log10_transform(sign(bss_model@pvalues) * p.adjust(abs(bss_model@pvalues), method = 'BY'))
  dim(log_pvalues_adjusted) <- dim(bss_data@atlas_image)
  outdir <- bss_out@outdir

  tvalues <- rep(0, length(bss_data@atlas_image))
  tvalues[bss_data@mask_idx] <- bss_model@tvalues
  dim(tvalues) <- dim(bss_data@atlas_image)

  if (bss_model@model_type == "bss_corr") {
    corr_values <- rep(0, length(bss_data@atlas_image))
    corr_values[bss_data@mask_idx] <- bss_model@corr_values
    dim(corr_values) <- dim(bss_data@atlas_image)
  }

  switch(bss_model@model_type,
         bss_anova = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@main_effect, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@main_effect, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues_adjusted, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(tvalues, bss_model@main_effect, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(tvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
         },

         bss_lm = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@main_effect, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@main_effect, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues_adjusted, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(tvalues, bss_model@main_effect, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(tvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
         },
         bss_corr = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@corr_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@corr_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues_adjusted, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(corr_values, bss_model@corr_var, "corr_values", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(corr_values, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
         },
         pairedttest = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@group_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@group_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues_adjusted, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(tvalues, bss_model@group_var, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(tvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
         },
         unpairedttest = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@group_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@group_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues_adjusted, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(tvalues, bss_model@group_var, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(tvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
         }
  )

    # Copy modelspec file to the output directory
  file.copy(bss_model@mspec_file, bss_out@outdir)
  invisible(bss_out)
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
  if (! bss_data@analysis_type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  switch(bss_data@analysis_type,
         cbm = { bss_out <- new("BssCBMOutput", outdir)},
         tbm = { bss_out <- new("BssTBMOutput", outdir) },
         roi = { bss_out <- new("BssROIOutput", outdir) }
  )
  bss_out <- save_out(bss_out, bss_data, bss_model)
  invisible(bss_out)
}

save_bss_color_files <- function(measure, var_name, cmap_title, bss_data, bss_model, outdir) {

  measure <- as.numeric(measure)
  bss_cmap <- new("BssColormap", cmap_title, "RdYlBu", measure)
  cbar_filename <- paste(paste(bss_model@model_type, var_name, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  save_colorbar(file.path(outdir,cbar_filename), bss_cmap@lut, bss_cmap@vmin, bss_cmap@vmax, cmap_title)

  # save the color LUT
  lut_fileprefix <- paste(paste(bss_model@model_type, var_name, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.lut', sep = '')
  save_BrainSuiteLUT(file.path(outdir, lut_fileprefix), bss_cmap@lut)

  return(bss_cmap)
}

save_bss_out_surface <- function(measure, var_name, bss_cmap, bss_data, bss_model, outdir) {

  s1 <- bss_data@atlas_surface
  s1$attributes <- measure
  s1$vColor <- bss_cmap@rgbcolors
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  outprefix <- paste(paste(bss_model@model_type, var_name, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), bss_data@data_type, sep = '')
  writedfs(file.path(outdir, outprefix), s1)
}

save_bss_out_nifti_image <- function(measure, var_name, bss_cmap, bss_data, bss_model, outdir) {

  outprefix <- paste(paste(bss_model@model_type, var_name, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), bss_data@data_type, sep = '')
  RNifti::writeNifti(measure, file.path(outdir, outprefix), template = bss_data@atlas_image)
}
