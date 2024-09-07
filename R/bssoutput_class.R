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

#' Check that subject directory and demographics files exist
#' @param object object of type `BstrOutput`
#'
check_files <- function(object){
  if (!dir.exists(object@subjdir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", object@subjdir), call. = FALSE)
  }

  if (!file.exists(object@csv)) {
    stop(sprintf("Demographics csv file %s does not exist.\n", object@csv), call. = FALSE)
  }
}

#' S4 class for saving results of statistical analysis
#' @slot outdir output directory to save the results
#'
BstrOutput <- setClass(
  "BstrOutput",
  slots = list(
    outdir = "character"
  ),
  validity = check_files
)

setMethod("initialize", valueClass = "BstrOutput", signature = "BstrOutput", function(.Object, outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir)
    .Object@outdir <- outdir
  }
  else {
    .Object@outdir <- outdir
    message(sprintf("The output directory %s already exists.", outdir))
  }
  return(.Object)
})

#' Generic save function for `BstrOutput`
#' @param bstr_out object of type `BstrOutput`
#' @param bstr_data object of type `BstrData`
#' @param bstr_model object of type `BstrModel`
#' @param overwrite logical parameter denoting if existing output directory should be overwritten or not (default is FALSE)
#' @param ... Extra named arguments passed to save_out
#' @details
#' For the most part, the user will never have to call this function directly.
#' Instead the user should call [save_bstr_out()].
#' @seealso [save_bstr_out()]
#' @export
setGeneric("save_out", valueClass = "BstrOutput", function(bstr_out, bstr_data, bstr_model, overwrite = FALSE, ...) {
  standardGeneric("save_out")
})

BstrSBAOutput <- setClass(
  "BstrSBAOutput",
  contains = "BstrOutput"
)

BstrTBMOutput <- setClass(
  "BstrTBMOutput",
  contains = "BstrOutput"
)

BstrDBAOutput <- setClass(
  "BstrDBAOutput",
  contains = "BstrOutput"
)

BstrROIOutput <- setClass(
  "BstrROIOutput",
  contains = "BstrOutput"
)


#' @rdname save_out
#' @inheritParams save_out
setMethod("save_out", valueClass = "BstrSBAOutput", signature = "BstrSBAOutput", function(bstr_out, bstr_data, bstr_model, overwrite = F) {

  # If output directory is not empty, then empty if overwrite is true or stop if overwrite is false
  if (overwrite == TRUE){
    delete_and_recreate_dir(bstr_out@outdir)
  }
  # else {
  #   if (length(list.files(bstr_out@outdir, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) != 0){
  #     stop(sprintf("Output directory %s is not empty.\n", bstr_out@outdir), call. = FALSE)
  #   }
  # }

  log_pvalues <- log10_transform(bstr_model@pvalues)
  outdir <- bstr_out@outdir
  log_pvalues_adjusted <- log10_transform(bstr_model@pvalues_adjusted)
  bstr_model@tvalues[abs(log_pvalues) <= -1*log10(0.05)] <- 0

  switch(bstr_model@model_type,
         bstr_anova = {
           bstr_cmap <- save_bstr_color_files(log_pvalues, bstr_model@main_effect, "log_pvalues", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(log_pvalues, bstr_model@main_effect, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(log_pvalues_adjusted, bstr_model@main_effect, "log_pvalues_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(log_pvalues_adjusted, bstr_model@main_effect, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@tvalues, bstr_model@main_effect, "tvalues", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@tvalues, bstr_model@main_effect, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@tvalues_adjusted, bstr_model@main_effect, "tvalues_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@tvalues_adjusted, bstr_model@main_effect, bstr_cmap, bstr_data, bstr_model, outdir)
           save_bstr_rds(bstr_model@pvalues, bstr_model@main_effect, "pvalues", bstr_data, bstr_model, outdir) # Save pvalues as a rds file
         },
         bstr_lm = {
           bstr_cmap <- save_bstr_color_files(log_pvalues, bstr_model@main_effect, "log_pvalues", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(log_pvalues, bstr_model@main_effect, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(log_pvalues_adjusted, bstr_model@main_effect, "log_pvalues_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(log_pvalues_adjusted, bstr_model@main_effect, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@tvalues, bstr_model@main_effect, "tvalues", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@tvalues, bstr_model@main_effect, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@tvalues_adjusted, bstr_model@main_effect, "tvalues_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@tvalues_adjusted, bstr_model@main_effect, bstr_cmap, bstr_data, bstr_model, outdir)
           save_bstr_rds(bstr_model@pvalues, bstr_model@main_effect, "pvalues", bstr_data, bstr_model, outdir) # Save pvalues as a rds file
           },
         bstr_corr = {
           bstr_model@corr_values[abs(log_pvalues) <= -1*log10(0.05)] <- 0
           bstr_cmap <- save_bstr_color_files(log_pvalues, bstr_model@corr_var, "log_pvalues", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(log_pvalues, bstr_model@corr_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(log_pvalues_adjusted, bstr_model@corr_var, "log_pvalues_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(log_pvalues_adjusted, bstr_model@corr_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@tvalues, bstr_model@corr_var, "tvalues", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@tvalues, bstr_model@corr_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@tvalues_adjusted, bstr_model@corr_var, "tvalues_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@tvalues_adjusted, bstr_model@corr_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@corr_values, bstr_model@corr_var, "corr_values", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@corr_values, bstr_model@corr_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@corr_values_masked_adjusted, bstr_model@corr_var, "corr_values_masked_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@corr_values_masked_adjusted, bstr_model@corr_var, bstr_cmap, bstr_data, bstr_model, outdir)
           save_bstr_rds(bstr_model@pvalues, bstr_model@corr_var, "pvalues", bstr_data, bstr_model, outdir) # Save pvalues as a rds file
           },
         pairedttest = {
           bstr_cmap <- save_bstr_color_files(log_pvalues, bstr_model@group_var, "log_pvalues", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(log_pvalues, bstr_model@group_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(log_pvalues_adjusted, bstr_model@group_var, "log_pvalues_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(log_pvalues_adjusted, bstr_model@group_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@tvalues, bstr_model@group_var, "tvalues", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@tvalues, bstr_model@group_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@tvalues_adjusted, bstr_model@group_var, "tvalues_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@tvalues_adjusted, bstr_model@group_var, bstr_cmap, bstr_data, bstr_model, outdir)
           save_bstr_rds(bstr_model@pvalues, bstr_model@group_var, "pvalues", bstr_data, bstr_model, outdir) # Save pvalues as a rds file
         },
         unpairedttest = {
           bstr_cmap <- save_bstr_color_files(log_pvalues, bstr_model@group_var, "log_pvalues", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(log_pvalues, bstr_model@group_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(log_pvalues_adjusted, bstr_model@group_var, "log_pvalues_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(log_pvalues_adjusted, bstr_model@group_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@tvalues, bstr_model@group_var, "tvalues", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@tvalues, bstr_model@group_var, bstr_cmap, bstr_data, bstr_model, outdir)
           bstr_cmap <- save_bstr_color_files(bstr_model@tvalues_adjusted, bstr_model@group_var, "tvalues_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_surface_both_hemi(bstr_model@tvalues_adjusted, bstr_model@group_var, bstr_cmap, bstr_data, bstr_model, outdir)
           save_bstr_rds(bstr_model@pvalues, bstr_model@group_var, "pvalues", bstr_data, bstr_model, outdir) # Save pvalues as a rds file
         }
  )

  # Copy modelspec file to the output directory
  file.copy(bstr_model@mspec_file, bstr_out@outdir)
  return(bstr_out)
  }
)


#' @rdname save_out
#' @inheritParams save_out
#' @param nclusters numeric parameter denoting number of clusters (default is 10)
setMethod("save_out", valueClass = "BstrTBMOutput", signature = "BstrTBMOutput", function(bstr_out, bstr_data, bstr_model, overwrite = F, nclusters = 10) {

  # If output directory is not empty, then empty if overwrite is true or stop if overwrite is false
  if (overwrite == TRUE){
    delete_and_recreate_dir(bstr_out@outdir)
  } else {
    if (length(list.files(bstr_out@outdir, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) != 0){
      stop(sprintf("Output directory %s is not empty.\n", bstr_out@outdir), call. = FALSE)
    }
  }

  log_pvalues <- rep(1, length(bstr_data@atlas_image))
  log_pvalues[bstr_data@mask_idx] <- log10_transform(bstr_model@pvalues)
  dim(log_pvalues) <- dim(bstr_data@atlas_image)

  log_pvalues_adjusted <- rep(1, length(bstr_data@atlas_image))
  log_pvalues_adjusted[bstr_data@mask_idx] <- log10_transform(bstr_model@pvalues_adjusted)
  dim(log_pvalues_adjusted) <- dim(bstr_data@atlas_image)
  outdir <- bstr_out@outdir

  tvalues <- rep(0, length(bstr_data@atlas_image))
  tvalues[bstr_data@mask_idx] <- bstr_model@tvalues
  dim(tvalues) <- dim(bstr_data@atlas_image)

  tvalues_adjusted <- rep(0, length(bstr_data@atlas_image))
  tvalues_adjusted[bstr_data@mask_idx] <- bstr_model@tvalues_adjusted
  dim(tvalues_adjusted) <- dim(bstr_data@atlas_image)
  measure <- NULL
  switch(bstr_model@model_type,
         bstr_anova = {
           var_name = bstr_model@main_effect
         },
         bstr_lm = {
           var_name = bstr_model@main_effect
         },
         bstr_corr = {
           corr_values <- rep(0, length(bstr_data@atlas_image))
           corr_values[bstr_data@mask_idx] <- bstr_model@corr_values
           dim(corr_values) <- dim(bstr_data@atlas_image)

           corr_values_masked_adjusted <- rep(0, length(bstr_data@atlas_image))
           corr_values_masked_adjusted[bstr_data@mask_idx] <- bstr_model@corr_values_masked_adjusted
           dim(corr_values_masked_adjusted) <- dim(bstr_data@atlas_image)

           var_name = bstr_model@corr_var
         },
         pairedttest = {
           var_name = bstr_model@group_var
         },
         unpairedttest = {
           var_name = bstr_model@group_var
         }
  )

  if (bstr_model@model_type == "bstr_corr") {
    save_vol_stats_out(log_pvalues, log_pvalues_adjusted, tvalues, tvalues_adjusted,
                       var_name, bstr_data, bstr_model, outdir, corr_values, corr_values_masked_adjusted)
  } else {
    save_vol_stats_out(log_pvalues, log_pvalues_adjusted, tvalues, tvalues_adjusted, var_name, bstr_data, bstr_model, outdir)
  }

  voxelcoord <- get_voxelcoord(bstr_out, bstr_data, bstr_model, outdir, nclusters)

  # Check if voxelcoord is empty
  if (length(voxelcoord) == 0 | voxelcoord[[1]][1] == -1) {
    sink(file.path(outdir,  sprintf("report_%s_%s.Rmd", bstr_model@model_type, var_name)), type = "output")
    cat("---\n")
    bstr_report_title_str <- create_rmd_report_title_str(bstr_data, bstr_model)
    cat(sprintf("title: Bstr Report -- %s \n", bstr_report_title_str))
    cat("output: html_document\n")
    cat("---\n\n\n")
    cat("No detected clusters above significance threshold.")
    sink()
    rmarkdown::render(file.path(outdir,  sprintf("report_%s_%s.Rmd", bstr_model@model_type, var_name)))
    stop("No detected clusters above significance threshold.")
  }

  # add create an R6 class function from here
  bstrmd_volout <- BstrRmdVolumeOutput$new()
  bstrmd_volout$save_out(bstr_data, bstr_model, voxelcoord = voxelcoord, outdir)


    # Copy modelspec file to the output directory
  file.copy(bstr_model@mspec_file, bstr_out@outdir)
  invisible(bstr_out)
  }
)


#' @rdname save_out
#' @inheritParams save_out
setMethod("save_out", valueClass = "BstrDBAOutput", signature = "BstrDBAOutput", function(bstr_out, bstr_data, bstr_model, overwrite = F, nclusters = 10) {

  # If output directory is not empty, then empty if overwrite is true or stop if overwrite is false
  if (overwrite == TRUE){
    delete_and_recreate_dir(bstr_out@outdir)
  } else {
    if (length(list.files(bstr_out@outdir, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) != 0){
      stop(sprintf("Output directory %s is not empty.\n", bstr_out@outdir), call. = FALSE)
    }
  }

  log_pvalues <- rep(1, length(bstr_data@atlas_image))
  log_pvalues[bstr_data@mask_idx] <- log10_transform(bstr_model@pvalues)
  dim(log_pvalues) <- dim(bstr_data@atlas_image)

  log_pvalues_adjusted <- rep(1, length(bstr_data@atlas_image))
  log_pvalues_adjusted[bstr_data@mask_idx] <- log10_transform(bstr_model@pvalues_adjusted)
  dim(log_pvalues_adjusted) <- dim(bstr_data@atlas_image)
  outdir <- bstr_out@outdir

  tvalues <- rep(0, length(bstr_data@atlas_image))
  tvalues[bstr_data@mask_idx] <- bstr_model@tvalues
  dim(tvalues) <- dim(bstr_data@atlas_image)

  tvalues_adjusted <- rep(0, length(bstr_data@atlas_image))
  tvalues_adjusted[bstr_data@mask_idx] <- bstr_model@tvalues_adjusted
  dim(tvalues_adjusted) <- dim(bstr_data@atlas_image)
  measure <- NULL
  switch(bstr_model@model_type,
         bstr_anova = {
           var_name = bstr_model@main_effect
         },
         bstr_lm = {
           var_name = bstr_model@main_effect
         },
         bstr_corr = {
           corr_values <- rep(0, length(bstr_data@atlas_image))
           corr_values[bstr_data@mask_idx] <- bstr_model@corr_values
           dim(corr_values) <- dim(bstr_data@atlas_image)

           corr_values_masked_adjusted <- rep(0, length(bstr_data@atlas_image))
           corr_values_masked_adjusted[bstr_data@mask_idx] <- bstr_model@corr_values_masked_adjusted
           dim(corr_values_masked_adjusted) <- dim(bstr_data@atlas_image)

           var_name = bstr_model@corr_var
         },
         pairedttest = {
           var_name = bstr_model@group_var
         },
         unpairedttest = {
           var_name = bstr_model@group_var
         }
  )

  if (bstr_model@model_type == "bstr_corr") {
    save_vol_stats_out(log_pvalues, log_pvalues_adjusted, tvalues, tvalues_adjusted,
                       var_name, bstr_data, bstr_model, outdir, corr_values, corr_values_masked_adjusted)
  } else {
    save_vol_stats_out(log_pvalues, log_pvalues_adjusted, tvalues, tvalues_adjusted, var_name, bstr_data, bstr_model, outdir)
  }

  voxelcoord <- get_voxelcoord(bstr_out, bstr_data, bstr_model, outdir, nclusters)

  # Check if voxelcoord is empty
  if (length(voxelcoord) == 0 | voxelcoord[[1]][1] == -1) {
    sink(file.path(outdir,  sprintf("report_%s_%s.Rmd", bstr_model@model_type, var_name)), type = "output")
    cat("---\n")
    bstr_report_title_str <- create_rmd_report_title_str(bstr_data, bstr_model)
    cat(sprintf("title: Bstr Report -- %s \n", bstr_report_title_str))
    cat("output: html_document\n")
    cat("---\n\n\n")
    cat("No detected clusters above significance threshold.")
    sink()
    rmarkdown::render(file.path(outdir,  sprintf("report_%s_%s.Rmd", bstr_model@model_type, var_name)))
    stop("No detected clusters above significance threshold.")
  }

  # Create a new R6 class object here
  bstrmd_volout <- BstrRmdVolumeOutput$new()
  bstrmd_volout$save_out(bstr_data, bstr_model, voxelcoord = voxelcoord, outdir)

  # Copy modelspec file to the output directory
  file.copy(bstr_model@mspec_file, bstr_out@outdir)
  invisible(bstr_out)
}
)

#' @rdname save_out
#' @inheritParams save_out
setMethod("save_out", valueClass = "BstrROIOutput", signature = "BstrROIOutput", function(bstr_out, bstr_data, bstr_model, overwrite = F) {

  # # Create the output directory
  # if (is.null(outdir)) {
  #   # Outputted directory (if not specified) is of the same format (csv or tsv) as the inputted file
  #   if (tools::file_ext(file.path(csv)) == "csv"){
  #     outdir <- paste0(tools::file_path_sans_ext(csv), "_roidata.csv")
  #     data_separator = ','
  #   } else if (tools::file_ext(file.path(csv)) == "tsv"){
  #     outdir <- paste0(tools::file_path_sans_ext(csv), "_roidata.tsv")
  #     data_separator = '\t'
  #   }
  #   cat(sprintf('Output directory is not specified. Using %s to save outputs.\n', outdir))
  # }
  # else {
  #   dir.create(file.path(outdir), showWarnings = FALSE)
  # }
  #
  # write.table(bstr_data@demographics, outdir, sep = data_separator,row.names = FALSE)

  # If output directory is not empty, then empty if overwrite is true or stop if overwrite is false
  if (overwrite == TRUE){
    delete_and_recreate_dir(bstr_out@outdir)
  } else {
    if (length(list.files(bstr_out@outdir, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) != 0){
      stop(sprintf("Output directory %s is not empty.\n", bstr_out@outdir), call. = FALSE)
    }
  }

  # Get the absolute path of outdir
  outdir <- tools::file_path_as_absolute(bstr_out@outdir)

  # Copy demographics csv file to output directory
  csvfilename <- paste0(file.path(bstr_out@outdir, tools::file_path_sans_ext(basename(bstr_data@csv))),".csv")
  write.csv(bstr_data@demographics, csvfilename)

  nb_header <- "---\ntitle: 'BrainSuite ROI statistical analysis report'\noutput: html_document\n---"

  nb_libraries <-"```{r librar_cmds, echo=FALSE}\n"
  nb_libraries <- paste(nb_libraries, "\nlibrary('bstr')\nlibrary('ggplot2')\n", sep="")
  nb_libraries <- paste(nb_libraries, "\n```\n", sep="")

  nb_data_header_one <- "The following command loads the data."

  # Add commands to load the data
  nb_data_command_one <- sprintf("\n```{r message=FALSE, warning=FALSE, results='hide', data_command_one}\n")
  nb_data_command_one <- paste(nb_data_command_one,"\n",
                            bstr_data@load_data_command,"\n",sep = "")
  nb_data_command_one <- paste(nb_data_command_one, "\n```\n", sep = "")


  nb_load_data <- sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, load_data}\n")
  nb_load_data <- paste0(nb_load_data, "\nDT::datatable(bstr_data@demographics[,-(ncol(bstr_data@demographics))], rownames = FALSE)\n")
  nb_load_data <- paste0(nb_load_data, "\n```\n")


  nb_data_header_two <- "The following command creates the model used to analyze the data."

  nb_data_command_two <- sprintf("\n```{r message=FALSE, warning=FALSE, results='hide', data_command_two}\n")
  nb_data_command_two <- paste(nb_data_command_two, "\n",bstr_model@load_data_command,"\n")
  nb_data_command_two <- paste0(nb_data_command_two, "\n```\n")

  nb_data_header_three <- "The final command (below) was used to render this document."

  command_save_bstr_out <- sprintf("save_bstr_out(bstr_data, bstr_model, outdir = '%s')", bstr_out@outdir)

  nb_data_command_three <-  sprintf("\n```{r eval=FALSE, message=FALSE, warning=FALSE, data_command_three}\n")
  nb_data_command_three <- paste(nb_data_command_three, "\n", command_save_bstr_out,"\n")
  nb_data_command_three <- paste0(nb_data_command_three, "\n```\n")


  nb_commands <- lapply(as.character(bstr_data@roiids),function(x){sprintf("```{r warning=FALSE, run_command_%s}\n",x)})
  nb_plots <- lapply(as.character(bstr_data@roiids),function(x){sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, plot_%s}\n",x)})
  if (bstr_model@model_type != 'bstr_corr'){
    nb_calculations <- lapply(as.character(bstr_data@roiids),function(x){sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, pval_%s}\n",x)})
  }

  selected_col <- rep(NA, length(bstr_data@roiids))
  for (i in 1:length(bstr_data@roiids)){
    selected_col[i] <- paste0(as.character(get_roi_tag(read_label_desc(),bstr_data@roiids[i])), "(",bstr_data@roiids[i],")")
    bstr_data@demographics[,selected_col[i]]
  }

  if (bstr_model@model_type != 'pairedttest' & bstr_model@model_type != 'unpairedttest') {
    if (bstr_model@model_type=="bstr_lmer"){ comparison_stat = "Chisq"}
    else {comparison_stat = "F"}

    for (m in 1:length(bstr_data@roiids)){
      if (bstr_model@model_type != 'bstr_corr'){
        for (i in 1:length(bstr_model@stats_commands[[m]])) {
          nb_commands[[m]] <- paste0(nb_commands[[m]], bstr_model@stats_commands[[m]][i], "\n")
        }
        nb_commands[[m]] <- paste0(nb_commands[[m]], "```\n\n")
        nb_calculations[[m]] <- paste0(nb_calculations[[m]],"anova_table <- ",
                                       substr(bstr_model@stats_commands[[m]][3],nchar("pander::pander(")+1,nchar(bstr_model@stats_commands[[m]][3])-1),
                                       "\np_val_",bstr_data@roiids[m],"<- round(anova_table$`Pr(>",comparison_stat,")`[2],digits=4)\n",
                                       "pval_string <- paste('pvalue:', as.character(p_val_",bstr_data@roiids[m],"))\n```\n\n")
      }
      if (inherits(bstr_data@demographics[,gsub("([[:alnum:]_]+).*", "\\1", bstr_model@fullmodel)], "integer") |
          inherits(bstr_data@demographics[,gsub("([[:alnum:]_]+).*", "\\1", bstr_model@fullmodel)], "double") |
          bstr_model@model_type == 'bstr_corr' | bstr_model@model_type == 'bstr_lm' | bstr_model@model_type == 'bstr_anova'){
        if (bstr_model@model_type == "bstr_corr"){
          x_var = bstr_model@corr_var
          if (bstr_model@group_var == "") {
            annotate_label = paste0("paste('corr val:',",round(bstr_model@corr_values[m],5),")")
            pdf_filename <- sprintf("%s_roi%d_%s_vs_%s.pdf", as.character(get_roi_tag(read_label_desc(),bstr_data@roiids[m])), bstr_data@roiids[m], bstr_data@roimeas, x_var)
            nb_plots[[m]] <- sprintf("ggplot2::ggplot(bstr_data@demographics, aes(x=%s, y=`%s`)) +
            ggplot2::geom_point(shape=21, size=1.5) + ggplot2::geom_smooth(method=lm, se=TRUE) +
            ggplot2::theme(aspect.ratio=1) + ggplot2::ggtitle('%s') + ggplot2::labs(y='%s') +
            ggplot2::theme(axis.title.x = element_text(size = rel(1.6))) + ggplot2::theme(axis.title.y = element_text(size = rel(1.6))) +
            ggplot2::theme(plot.title= ggplot2::element_text(size=20), hjust = 0.5) +
            ggplot2::theme(axis.text.x=ggplot2::element_text(size=14,face='bold')) +
            ggplot2::theme(axis.text.y=ggplot2::element_text(size=14,face='bold')) +
            ggplot2::scale_x_continuous(expand = c(0.2, 0.2)) +
            ggplot2::scale_y_continuous(expand = c(0.2, 0.2)) +
            paletteer::scale_color_paletteer_d('ggsci::category10_d3') +
        		paletteer::scale_fill_paletteer_d('ggsci::category10_d3')\n
            ggplot2::ggsave(filename='%s')\n```\n",
                                       x_var, selected_col[m],
                                       as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bstr_data@roiids[m])[[1]]),
                                       bstr_data@roimeas,
                                       pdf_filename)
          }
          else {
            group_var = bstr_model@group_var
            annotate_label = paste0("paste('corr val:',",round(bstr_model@corr_values[m],5),")")
            pdf_filename <- sprintf("%s_roi%d_%s_vs_%s.pdf", as.character(get_roi_tag(read_label_desc(),bstr_data@roiids[m])), bstr_data@roiids[m], bstr_data@roimeas, x_var)
            nb_plots[[m]] <- sprintf("ggplot2::ggplot(bstr_data@demographics, aes(x=%s, y=`%s`, color=%s)) +
            ggplot2::geom_point(shape=21, size=1.5, aes(fill=%s)) + ggplot2::geom_smooth(aes(fill=%s),  method=lm, se=TRUE) +
            ggplot2::theme(aspect.ratio=1) + ggplot2::ggtitle('%s') + ggplot2::labs(y='%s') + ggplot2::labs(x='%s') +
            ggplot2::theme(axis.title.x = element_text(size = rel(1.6))) + ggplot2::theme(axis.title.y = element_text(size = rel(1.6))) +
            ggplot2::theme(plot.title= ggplot2::element_text(size=20)) +
            ggplot2::theme(axis.text.x=ggplot2::element_text(size=14,face='bold')) +
            ggplot2::theme(axis.text.y=ggplot2::element_text(size=14,face='bold')) +
            ggplot2::scale_x_continuous(expand = c(0.2, 0.2)) +
            ggplot2::scale_y_continuous(expand = c(0.2, 0.2)) +
            paletteer::scale_color_paletteer_d('ggsci::category10_d3') +
        		paletteer::scale_fill_paletteer_d('ggsci::category10_d3')\n
            ggplot2::ggsave(filename='%s')\n```\n",
                                       x_var, selected_col[m], bstr_model@group_var,
                                       bstr_model@group_var, bstr_model@group_var,
                                       as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bstr_data@roiids[m])[[1]]), bstr_data@roimeas, x_var,
                                       pdf_filename)
          }

        } else if (bstr_model@model_type == "bstr_lm" || bstr_model@model_type == 'bstr_anova') {
          x_var = gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel)
          annotate_label = paste0("paste('pvalue:', as.character(p_val_",bstr_data@roiids[m],"))")
          pdf_filename <- sprintf("%s_roi%d_%s_vs_%s.pdf", as.character(get_roi_tag(read_label_desc(),bstr_data@roiids[m])), bstr_data@roiids[m], bstr_data@roimeas, x_var)
          if (is.factor(bstr_data@demographics[[bstr_model@main_effect]]) ){
            bstr_data@demographics[[bstr_model@main_effect]] <- as.factor(bstr_data@demographics[[bstr_model@main_effect]])
            nb_plots[[m]] <- sprintf("```{r echo=TRUE, warning=FALSE, message=FALSE}\n\n ggplot2::ggplot(bstr_data@demographics, aes(x=as.factor(%s), y=`%s`, fill=%s)) +
            ggplot2::geom_boxplot(width = 0.4, alpha=0.5, position=position_dodge(width = 0.4), outlier.colour='black')  +
      		  geom_violin(alpha = 0.5, trim=FALSE, position = position_dodge(width = 0.4)  ) +
      		  ggplot2::geom_point(shape=21, position = position_jitterdodge(seed = 1, dodge.width = 0.4), size=0.3) +
            ggplot2::theme(aspect.ratio=1) + ggplot2::ggtitle('%s') + ggplot2::labs(y='%s') + ggplot2::labs(x='%s') +
            ggplot2::theme(axis.title.x = element_text(size = rel(1.6))) + ggplot2::theme(axis.title.y = element_text(size = rel(1.6))) +
            ggplot2::theme(plot.title= ggplot2::element_text(size=20, hjust = 0.5)) +
            ggplot2::theme(axis.text.x=ggplot2::element_text(size=14,face='bold')) +
            ggplot2::theme(axis.text.y=ggplot2::element_text(size=14,face='bold')) +
            ggplot2::scale_x_discrete(expand = c(0.4, 0.4)) +
            ggplot2::scale_y_continuous(expand = c(0.25, 0.25)) +
            paletteer::scale_color_paletteer_d('ggsci::category10_d3') +
        	  paletteer::scale_fill_paletteer_d('ggsci::category10_d3')\n
            ggplot2::ggsave(filename='%s')\n```\n",
                                       x_var, selected_col[m], bstr_model@main_effect,
                                       as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bstr_data@roiids[m])[[1]]),
                                       bstr_data@roimeas, bstr_model@main_effect,
                                       pdf_filename)
          }
          else {
            # TODO Add provision to group by a covariate, which is a factor.
            # Will involve testing each covariate whether it can be converted into a factor.
            # Choose that covariate with the minimum number of levels
            nb_plots[[m]] <- sprintf("```{r echo=TRUE, warning=FALSE, message=FALSE}\n\n ggplot2::ggplot(bstr_data@demographics, aes(x=%s, y=`%s`)) +
            ggplot2::geom_point(shape=21, size=1.5) + ggplot2::geom_smooth(method=lm, se=TRUE) +
            ggplot2::theme(aspect.ratio=1) + ggplot2::ggtitle('%s') + ggplot2::labs(y='%s') +
            ggplot2::theme(axis.title.x = element_text(size = rel(1.6))) + ggplot2::theme(axis.title.y = element_text(size = rel(1.6))) +
            ggplot2::theme(plot.title= ggplot2::element_text(size=20, hjust = 0.5)) +
            ggplot2::theme(axis.text.x=ggplot2::element_text(size=14,face='bold')) +
            ggplot2::theme(axis.text.y=ggplot2::element_text(size=14,face='bold')) +
            ggplot2::scale_x_continuous(expand = c(0.2, 0.2)) +
            ggplot2::scale_y_continuous(expand = c(0.2, 0.2)) +
            paletteer::scale_color_paletteer_d('ggsci::category10_d3') +
        		paletteer::scale_fill_paletteer_d('ggsci::category10_d3')\n
            ggplot2::ggsave(filename='%s')\n```\n",
                                       x_var, selected_col[m],
                                       as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bstr_data@roiids[m])[[1]]),
                                       bstr_data@roimeas,
                                       pdf_filename)
          }
        }

      } else if (class(bstr_data@demographics[,gsub("([[:alnum:]_]+).*", "\\1", bstr_model@fullmodel)])=="factor"){
        nb_plots[[m]]<-paste0(nb_plots[[m]],
                              "mean_lengths <- rep(NA, length(levels(bstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel),")))
                              std_devs <- rep(NA, length(levels(bstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel),")))
                              for (i in 1:length(levels(bstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel),"))){
                              mean_lengths[i] <-  mean(bstr_data@demographics[bstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel),
                              " == levels(bstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel),")[i],]$`",
                              selected_col[m],"`)
                              std_devs[i] <- sd (bstr_data@demographics[bstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel),
                              " == levels(bstr_data@demographics$",gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel),")[i],]$`",
                              selected_col[m],"`)
                              } \n",
"modified_df <- data.frame(len = mean_lengths, sd = std_devs,",gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel)," = levels(bstr_data@demographics$",
gsub('([A-Z_a-z]+).*', '\\1', bstr_model@fullmodel),"))\n",
"ggplot2::ggplot(data=modified_df, ggplot2::aes(x = ",gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel),
", y = len, color = ",gsub('([[:alnum:]_]+).*', '\\1', bstr_model@fullmodel),
")) +
ggplot2::geom_bar(stat = 'identity') + ggplot2::ggtitle('",
as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bstr_data@roiids[m])[[1]]),
" ", bstr_data@roimeas,
" vs ", bstr_model@main_effect,"') +
ggplot2::labs(y='",as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bstr_data@roiids[m])[[1]]),
"') +
ggplot2::geom_errorbar(mapping=aes(ymin=len-sd, ymax=len+sd)) +
ggplot2::theme(axis.title=ggplot2::element_text(size=16,face='bold')) +
ggplot2::theme(plot.title=ggplot2::element_text(size=18,face='bold')) +
ggplot2::scale_color_discrete(name= pval_string)\n
ggplot2::ggsave(filename='",
paste0(as.character(get_roi_tag(read_label_desc(),bstr_data@roiids[m]))), "_roi",bstr_data@roiids[m],
" ", bstr_data@roimeas,
" vs ", bstr_model@main_effect, ".pdf',device='pdf')\n```\n")
      }
      if (bstr_model@model_type == 'bstr_corr'){
        nb_commands[[m]] <- paste0(sprintf("\n#### Main effect of %s (%d) %s on %s \n",
                                           as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bstr_data@roiids[m])[[1]]),
                                           bstr_data@roiids[m],
                                           bstr_data@roimeas, bstr_model@corr_var),
                                   nb_commands[[m]])
      } else{
        nb_commands[[m]] <- paste0(sprintf("\n#### Main effect of %s (%d) %s on %s controlling for %s \n",
                                           as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bstr_data@roiids[m])[[1]]),
                                           bstr_data@roiids[m],
                                           bstr_data@roimeas, bstr_model@main_effect,bstr_model@covariates ),
                                   nb_commands[[m]])
      }

      }
  } else if (bstr_model@model_type == 'pairedttest' ||  bstr_model@model_type == 'unpairedttest' ) {
    for (i in 1:length(bstr_data@roiids)){
      current_t_test_table<-paste0("data.frame('P-Values' = ",round(bstr_model@pvalues[i],6),",'T-Values' = ",round(bstr_model@tvalues[i],6),
                                       ",'Adjusted P-Values' = ",round(bstr_model@pvalues_adjusted[i],6), ",'Adjusted T-Values' = ",round(bstr_model@tvalues_adjusted[i],6),
                                   ",row.names = paste0('roiid ',",bstr_data@roiids[i],"))")
      nb_calculations[[i]] <- paste0(sprintf("\n#### T-test output for differences between means of brain imaging phenotypes for %s for roiid %d \n",
                                             bstr_model@group_var,bstr_data@roiids[i]),
                                     nb_calculations[[i]],"\nDT::formatStyle(DT::datatable(",current_t_test_table,
                                     "),column = 'P.Values',color = ifelse(",abs(round(bstr_model@pvalues[i],6)),
                                     "<=0.05,'red','black'))\n```\n")
      nb_plots[[i]]<-paste0(nb_plots[[i]],"ggplot2::ggplot(data=bstr_data@demographics, ggplot2::aes(x=factor(",
                            bstr_model@group_var,"), y = `",selected_col[i],"`, color = factor(",bstr_model@group_var,"))) + ggplot2::geom_boxplot() + ggplot2::ggtitle('",
                            as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bstr_data@roiids[i])[[1]]),
                            " ", bstr_data@roimeas,
                            " vs ", bstr_model@group_var, "') + ggplot2::labs(x = '",bstr_model@group_var,"', y='",
                            as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bstr_data@roiids[i])[[1]]),
                            "', colour = '",bstr_model@group_var,"') +
                            ggplot2::theme(axis.title=ggplot2::element_text(size=16,face='bold')) +
                            ggplot2::theme(plot.title=ggplot2::element_text(size=18,face='bold')) \n
                            ggplot2::ggsave(filename='",
                            paste0(as.character(get_roi_tag(read_label_desc(),bstr_data@roiids[i]))), "_roi",bstr_data@roiids[i],
                            "_", bstr_data@roimeas,
                            "_vs_", bstr_model@group_var, ".pdf',device='pdf')\n```\n")
    }
  }



  rmdfileconn<-file(file.path(outdir, sprintf("report_%s_%s.Rmd", bstr_model@model_type, switch(bstr_model@model_type,
                                                                                               bstr_corr = {bstr_model@corr_var},
                                                                                               unpairedttest = {bstr_model@group_var},
                                                                                               pairedttest = {bstr_model@group_var},
                                                                                               bstr_anova = {bstr_model@main_effect},
                                                                                               bstr_lm = {bstr_model@main_effect},
                                                                                               bstr_lmer = {bstr_model@main_effect}))))
  if (bstr_model@model_type == 'bstr_corr'){
    writeLines(c(nb_header, nb_libraries, nb_data_header_one, nb_data_command_one,
                 nb_load_data, nb_data_header_two, nb_data_command_two,
                 nb_data_header_three, nb_data_command_three,
                 unlist(mapply(paste0, nb_commands, nb_plots))), rmdfileconn)
  } else {
    writeLines(c(nb_header, nb_libraries, nb_data_header_one, nb_data_command_one,
                 nb_load_data, nb_data_header_two, nb_data_command_two,
                 nb_data_header_three, nb_data_command_three,
                 unlist(mapply(paste0, nb_commands, nb_calculations, nb_plots))), rmdfileconn)
  }
  close(rmdfileconn)

  # Render the markdown
  rmarkdown::render(file.path(outdir, sprintf("report_%s_%s.Rmd", bstr_model@model_type, switch(bstr_model@model_type,
                                                                                               bstr_corr = {bstr_model@corr_var},
                                                                                               unpairedttest = {bstr_model@group_var},
                                                                                               pairedttest = {bstr_model@group_var},
                                                                                               bstr_anova = {bstr_model@main_effect},
                                                                                               bstr_lm = {bstr_model@main_effect},
                                                                                               bstr_lmer = {bstr_model@main_effect}))),
                    output_file=file.path(outdir, sprintf("report_%s_%s.html", bstr_model@model_type, switch(bstr_model@model_type,
                                                                                                           bstr_corr = {bstr_model@corr_var},
                                                                                                           unpairedttest = {bstr_model@group_var},
                                                                                                           pairedttest = {bstr_model@group_var},
                                                                                                           bstr_anova = {bstr_model@main_effect},
                                                                                                           bstr_lm = {bstr_model@main_effect},
                                                                                                           bstr_lmer = {bstr_model@main_effect}))), quiet = TRUE)

  # Copy modelspec file to the output directory
  file.copy(bstr_model@mspec_file, bstr_out@outdir)
  return(bstr_out)
  }
)

#' Save the statistical analysis output
#' @param bstr_data object of type `BstrData`
#' @param bstr_model object of type `BstrModel`
#' @param outdir output directory to save the results
#' @param overwrite logical parameter denoting if existing output directory should be overwritten or not (default is FALSE)
#' @param nclusters numeric value denoting number of clusters (default is 10)
#' @export
save_bstr_out <- function(bstr_data, bstr_model, outdir="", overwrite = F, nclusters = 10) {

  outdir <- path.expand(outdir)
  valid_types <- c("sba", "tbm", "roi", "dba", "nca")
  if (! bstr_data@analysis_type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  switch(bstr_data@analysis_type,
         sba = { bstr_out <- new("BstrSBAOutput", outdir)},
         tbm = { bstr_out <- new("BstrTBMOutput", outdir) },
         dba = { bstr_out <- new("BstrDBAOutput", outdir) },
         roi = { bstr_out <- new("BstrROIOutput", outdir) }
  )
  if (bstr_data@analysis_type == "tbm" | bstr_data@analysis_type == "dba"){
    bstr_out <- save_out(bstr_out, bstr_data, bstr_model, overwrite = overwrite, nclusters = nclusters)
    invisible(bstr_out)
  # } else if (bstr_data@analysis_type == "sba") {
  #   bstr_out <- save_bstr_out_sba_both_hemi(bstr_out, bstr_data, bstr_model, overwrite = overwrite)
  #   invisible(bstr_out)
  }
  else {
    bstr_out <- save_out(bstr_out, bstr_data, bstr_model, overwrite = overwrite)
    invisible(bstr_out)
  }
}


save_bstr_out_sba_both_hemi <- function(bstr_out, bstr_data, bstr_model, outdir="", overwrite = F) {

  if (bstr_data@hemi == "both") {
    # Split the data array, model and all variables into left and right hemispheres
    save_out(bstr_out, bstr_data, bstr_model, overwrite = overwrite)
    # bstr_data_lh <- bstr_data
    # bstr_data_lh@atlas_filename = bstr_data@atlas_filename_lh
    # bstr_data_lh@atlas_surface <- bstr_data_both@atlas_surface_lh
    # save_out(bstr_out, bstr_data_lh, bstr_model, overwrite = T)
    #
    # bstr_data_rh <- bstr_data
    # bstr_data_rh@atlas_filename = bstr_data@atlas_filename_rh
    # bstr_data_rh@atlas_surface <- bstr_data_both@atlas_surface_rh
    # save_out(bstr_out, bstr_data_rh, bstr_model, overwrite = F)

  }
  else {
    save_out(bstr_out, bstr_data, bstr_model, overwrite = overwrite)
  }
}

#' Save the color LUT and ini colormap images
#' @param measure numeric value denoting the measures used to create the color file
#' @param var_name string denoting name of variable that the color file is being created for
#' @param cmap_title string denoting the type of color map
#' @param bstr_data object of type `BstrData`
#' @param bstr_model object of type `BstrModel`
#' @param outdir string specifying output directory to save the results in
#' @export

save_bstr_color_files <- function(measure, var_name, cmap_title, bstr_data, bstr_model, outdir) {

  measure <- as.numeric(measure)
  bstr_cmap <- new("BstrColormap", cmap_title, "RdYlBu", measure)
  colorbar_label_text <- bs_stat_overlays_mapping_to_label[[cmap_title]]

  if (bstr_data@hemi == "both")
    cbar_filename <- paste(paste(bstr_model@model_type, var_name, 'both_hemi', bstr_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  else
      cbar_filename <- paste(paste(bstr_model@model_type, var_name, tools::file_path_sans_ext(
      basename(bstr_data@atlas_filename)), bstr_cmap@cmap_type, sep = '_'), '_cbar.pdf', sep = '')
  save_colorbar(file.path(outdir,cbar_filename), bstr_cmap@lut, bstr_cmap@vmin, bstr_cmap@vmax, colorbar_label_text)

  if (bstr_data@hemi == "both")
    cbar_filename <- paste(paste(bstr_model@model_type, var_name, 'both_hemi', bstr_cmap@cmap_type, sep = '_'), '_cbar.png', sep = '')
  else
    cbar_filename <- paste(paste(bstr_model@model_type, var_name, tools::file_path_sans_ext(
      basename(bstr_data@atlas_filename)), bstr_cmap@cmap_type, sep = '_'), '_cbar.png', sep = '')

  save_colorbar(file.path(outdir,cbar_filename), bstr_cmap@lut, bstr_cmap@vmin, bstr_cmap@vmax, colorbar_label_text)

  # save the color LUT
  if (bstr_data@hemi == "both")
    lut_fileprefix <- paste(paste(bstr_model@model_type, var_name, 'both_hemi', bstr_cmap@cmap_type, sep = '_'), '.lut', sep = '')
  else
    lut_fileprefix <- paste(paste(bstr_model@model_type, var_name, tools::file_path_sans_ext(
      basename(bstr_data@atlas_filename)), bstr_cmap@cmap_type, sep = '_'), '.lut', sep = '')

  save_BrainSuiteLUT(file.path(outdir, lut_fileprefix), bstr_cmap@lut)

  # save the cbar json file with colormap for BrainSuite statmap
  #cbarlist <- list("colorbar", bstr_cmap@vmin, bstr_cmap@vmax, t(col2rgb(bstr_cmap@lut)/255))



  cbarlist <- list("colorbar", bstr_cmap@cnegmax, bstr_cmap@cnegmin, bstr_cmap@cposmin, bstr_cmap@cposmax, t(col2rgb(bstr_cmap@lut)/255))
  #names(cbarlist) <- c("jsonid", "cbarmin", "cbarmax", "colormap")
  names(cbarlist) <- c("jsonid", "cbarmin", "cbarlowerthresh", "cbarupperthresh", "cbarmax", "colormap")
  cbar_json <- jsonlite::toJSON(cbarlist, pretty = TRUE, auto_unbox=TRUE)
  if (bstr_data@hemi == "both")
    cbar_json_filename <- paste(paste(bstr_model@model_type, var_name, 'both_hemi', bstr_cmap@cmap_type, sep = '_'), '.cbar', sep = '')
  else
    cbar_json_filename <- paste(paste(bstr_model@model_type, var_name, tools::file_path_sans_ext(
      basename(bstr_data@atlas_filename)), bstr_cmap@cmap_type, sep = '_'), '.cbar', sep = '')

  write(cbar_json, file.path(outdir,cbar_json_filename))

  # save ini colormap with ranges
  if (bstr_data@hemi == "both")
    ini_fileprefix <- paste(paste(bstr_model@model_type, var_name, 'both_hemi', bstr_cmap@cmap_type, sep = '_'), '.ini', sep = '')
  else
    ini_fileprefix <- paste(paste(bstr_model@model_type, var_name, tools::file_path_sans_ext(
      basename(bstr_data@atlas_filename)), bstr_cmap@cmap_type, sep = '_'), '.ini', sep = '')

  save_colormap_to_ini(file.path(outdir, ini_fileprefix), bstr_cmap)
  return(bstr_cmap)
}

#' Save the surface output to the given ouput directory
#' @param measure numeric value denoting the measures used to create the output
#' @param var_name string denoting name of variable used by the function
#' @param bstr_data object of type `BstrData`
#' @param bstr_model object of type `BstrModel`
#' @param bstr_cmap object of type `BstrColormap`
#' @param outdir string specifying output directory to save the results in
#' @export

save_bstr_out_surface <- function(measure, var_name, bstr_cmap, bstr_data, bstr_model, outdir) {

  s1 <- bstr_data@atlas_surface
  s1$attributes <- measure
  s1$vColor <- bstr_cmap@rgbcolors
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  outprefix <- paste(paste(bstr_model@model_type, var_name, tools::file_path_sans_ext(
    basename(bstr_data@atlas_filename)), bstr_cmap@cmap_type, sep = '_'), bstr_data@data_type, sep = '')
  writedfs(file.path(outdir, outprefix), s1)
}

save_bstr_out_surface_both_hemi <- function(measure, var_name, bstr_cmap, bstr_data, bstr_model, outdir) {

  if (bstr_data@hemi == "both") {
    idx_lh <- 1:bstr_data@nvertices_lh
    idx_rh <- (bstr_data@nvertices_lh+1):(bstr_data@nvertices_lh+bstr_data@nvertices_rh)

    bstr_data_lh <- bstr_data
    bstr_data_lh@atlas_filename = bstr_data@atlas_filename_lh
    bstr_data_lh@atlas_surface <- bstr_data@atlas_surface_lh
    measure_lh <- measure[idx_lh]
    bstr_cmap_lh <- bstr_cmap
    bstr_cmap_lh@values <- bstr_cmap@values[idx_lh]
    bstr_cmap_lh@rgbcolors <- bstr_cmap@rgbcolors[idx_lh,1:3]
    save_bstr_out_surface(measure_lh, var_name, bstr_cmap_lh, bstr_data_lh, bstr_model, outdir)

    bstr_data_rh <- bstr_data
    bstr_data_rh@atlas_filename = bstr_data@atlas_filename_rh
    bstr_data_rh@atlas_surface <- bstr_data@atlas_surface_rh
    measure_rh <- measure[idx_rh]
    bstr_cmap_rh <- bstr_cmap
    bstr_cmap_rh@values <- bstr_cmap@values[idx_rh]
    bstr_cmap_rh@rgbcolors <- bstr_cmap@rgbcolors[idx_rh,1:3]
    save_bstr_out_surface(measure_rh, var_name, bstr_cmap_rh, bstr_data_rh, bstr_model, outdir)

  }
  else {
      save_bstr_out_surface(measure, var_name, bstr_cmap, bstr_data, bstr_model, outdir)
  }
}
#' Save the nifti image to the output file
#' @param measure denotes the measure used to create the output
#' @param var_name string denoting name of variable used by the function
#' @param bstr_data object of type `BstrData`
#' @param bstr_model object of type `BstrModel`
#' @param bstr_cmap object of type `BstrColormap`
#' @param outdir string specifying output directory to save the results in
#' @export

save_bstr_out_nifti_image <- function(measure, var_name, bstr_cmap, bstr_data, bstr_model, outdir) {

  outprefix <- paste0(paste(bstr_model@model_type, var_name, tools::file_path_sans_ext(
    basename(bstr_data@atlas_filename)), bstr_cmap@cmap_type, sep = '_'), bstr_data@data_type)
  RNifti::writeNifti(measure, file.path(outdir, outprefix), template = bstr_data@atlas_image)
}

#' Save the measure to the output directory
#' @param measure numeric value denoting the measures used to create the output
#' @param var_name string denoting name of variable used by the function
#' @param label string denoting the label for the object
#' @param bstr_data object of type `BstrData`
#' @param bstr_model object of type `BstrModel`
#' @param outdir string specifying output directory to save the results in
#' @export

save_bstr_rds <- function(measure, var_name, label, bstr_data, bstr_model, outdir) {

  outprefix <- paste(paste(bstr_model@model_type, var_name, tools::file_path_sans_ext(
    basename(bstr_data@atlas_filename)), label, sep = '_'), ".rds", sep = '')
  saveRDS(measure, file=file.path(outdir, outprefix))
}

#' Save the volume statistics
#' @param log_pvalues log transformed p-values
#' @param log_pvalues_adjusted log transformed adjusted p-values
#' @param tvalues t-values
#' @param tvalues_adjusted adjusted t-values
#' @param var_name string denoting name of variable used by the function
#' @param bstr_data object of type `BstrData`
#' @param bstr_model object of type `BstrModel`
#' @param outdir string specifying output directory to save the results in
#' @param corr_values correlation values
#' @param corr_values_masked_adjusted adjusted, masked correlation values
#' @export

save_vol_stats_out <- function(log_pvalues, log_pvalues_adjusted, tvalues, tvalues_adjusted,
                                var_name, bstr_data, bstr_model, outdir, corr_values = NULL, corr_values_masked_adjusted = NULL) {

  bstr_cmap <- save_bstr_color_files(log_pvalues, var_name, "log_pvalues", bstr_data, bstr_model, outdir)
  save_bstr_out_nifti_image(log_pvalues, var_name, bstr_cmap, bstr_data, bstr_model, outdir)

  bstr_cmap <- save_bstr_color_files(log_pvalues_adjusted, var_name, "log_pvalues_adjusted", bstr_data, bstr_model, outdir)
  save_bstr_out_nifti_image(log_pvalues_adjusted, var_name, bstr_cmap, bstr_data, bstr_model, outdir)

  bstr_cmap <- save_bstr_color_files(tvalues, var_name, "tvalues", bstr_data, bstr_model, outdir)
  save_bstr_out_nifti_image(tvalues, var_name, bstr_cmap, bstr_data, bstr_model, outdir)

  bstr_cmap <- save_bstr_color_files(tvalues_adjusted, var_name, "tvalues_adjusted", bstr_data, bstr_model, outdir)
  save_bstr_out_nifti_image(tvalues_adjusted, var_name, bstr_cmap, bstr_data, bstr_model, outdir)
  save_bstr_rds(bstr_model@pvalues, var_name, "pvalues", bstr_data, bstr_model, outdir) # Save pvalues as a rds file

  switch(bstr_model@model_type,
         bstr_corr = {
           bstr_cmap <- save_bstr_color_files(corr_values, var_name, "corr_values", bstr_data, bstr_model, outdir)
           save_bstr_out_nifti_image(corr_values, var_name, bstr_cmap, bstr_data, bstr_model, outdir)

           bstr_cmap <- save_bstr_color_files(corr_values_masked_adjusted, var_name, "corr_values_masked_adjusted", bstr_data, bstr_model, outdir)
           save_bstr_out_nifti_image(corr_values_masked_adjusted, var_name, bstr_cmap, bstr_data, bstr_model, outdir)
         }
  )
}

#' Get voxel coordinates of all significant clusters (up to number of clusters)
#' @param bstr_out object of type `BstrOut`
#' @param bstr_data object of type `BstrData`
#' @param bstr_model object of type `BstrModel`
#' @param outdir string specifying output directory to save the results in
#' @param nclusters numeric value specifying number of clusters
#' @export

get_voxelcoord <- function(bstr_out, bstr_data, bstr_model, outdir, nclusters){
  # Call cluster code from terminal
  switch(bstr_model@model_type,
         bstr_anova = {var_name = bstr_model@main_effect},
         bstr_lm = {var_name = bstr_model@main_effect},
         bstr_lmer = {var_name = bstr_model@main_effect},
         bstr_corr = {var_name = bstr_model@corr_var},
         pairedttest = {var_name = bstr_model@group_var},
         unpairedttest = {var_name = bstr_model@group_var}
  )

  if (bstr_model@model_type == "bstr_corr") {
    file_suffix_adjusted <- bs_stat_overlays$corr_values_masked_adjusted
    file_suffix <- bs_stat_overlays$corr_values
  } else {
    file_suffix_adjusted <- bs_stat_overlays$tvalues_adjusted
    file_suffix <- bs_stat_overlays$tvalues
  }

  voxelcoord_call <- paste0("\"",file.path(get_brainsuite_install_path(),bs_binary_files$clustermap),"\" -i \"", outdir,"/",bstr_model@model_type,
                            "_",var_name,"_", tools::file_path_sans_ext(basename(bstr_data@atlas_filename)), "_", file_suffix_adjusted, ".nii.gz\""," -m \"", bstr_data@maskfile, "\" -o \"", outdir, "/cluster.tsv\"", " -n ", nclusters)

  # voxelcoord_call <- paste0("\"",file.path(get_brainsuite_install_path(),bs_binary_files$clustermap),"\" -i \"", outdir,"/",bstr_model@model_type,
  #                           "_",var_name,"_", tools::file_path_sans_ext(basename(bstr_data@atlas_filename)),
  #                           "_tvalues_adjusted.nii.gz\""," -m \"", bstr_data@maskfile, "\" -o \"", outdir, "/cluster.tsv\"", " -n ", nclusters)
  system(voxelcoord_call,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
  if(file.info(paste0(outdir, "/cluster.tsv"))$size!=0){
    vox_table <- read.table(paste0(outdir,"/cluster.tsv"),header=F,sep="\t")
    voxelcoord <- vector("list",nrow(vox_table))
    for (individ_vox in 1:nrow(vox_table)){
      for (vox_component in 1:3){
        voxelcoord[[individ_vox]]<- c(voxelcoord[[individ_vox]],vox_table[individ_vox,3+vox_component])
      }
    }
  }
  # Use tvalues instead of adjusted tvalues if no clusters are found
  if (file.info(paste0(outdir, "/cluster.tsv"))$size==0){
    voxelcoord_call <- paste0("\"",file.path(get_brainsuite_install_path(),bs_binary_files$clustermap),"\" -i \"", outdir,"/",bstr_model@model_type,
                              "_",var_name,"_", tools::file_path_sans_ext(basename(bstr_data@atlas_filename)), "_", file_suffix_adjusted, ".nii.gz\""," -m \"", bstr_data@maskfile, "\" -o \"", outdir, "/cluster.tsv\"", " -n ", nclusters)

#    voxelcoord_call <- paste0("\"",file.path(get_brainsuite_install_path(),bs_binary_files$clustermap),"\" -i \"", outdir,"/",bstr_model@model_type,
#                              "_",var_name,"_", tools::file_path_sans_ext(basename(bstr_data@atlas_filename)), "_tvalues.nii.gz\"",
#                              " -m \"", bstr_data@maskfile, "\" -o \"", outdir, "/cluster.tsv\"", " -n ", nclusters)
    system(voxelcoord_call,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
    if(file.info(paste0(outdir, "/cluster.tsv"))$size!=0){
      vox_table <- read.table(paste0(outdir,"/cluster.tsv"),header=F,sep="\t")
      voxelcoord <- vector("list",nrow(vox_table))
      for (individ_vox in 1:nrow(vox_table)){
        for (vox_component in 1:3){
          voxelcoord[[individ_vox]]<- c(voxelcoord[[individ_vox]],vox_table[individ_vox,3+vox_component])
        }
      }
    } else {voxelcoord <- list(-1)}
  }
  return(voxelcoord)
}
