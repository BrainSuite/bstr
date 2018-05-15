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

#' Generic save function for \code{BssOutput}
#' @param bss_out object of type \code{BssOutput}
#' @param bss_data object of type \code{BssData}
#' @param bss_model object of type \code{BssModel}
#' @details
#' For the most part, the user will never have to call this function directly.
#' Instead the user should call \code{\link{save_bss_out}}.
#' @seealso \code{\link{save_bss_out}}
#'
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

BssDBMOutput <- setClass(
  "BssDBMOutput",
  contains = "BssOutput"
)

BssROIOutput <- setClass(
  "BssROIOutput",
  contains = "BssOutput"
)

#' @rdname save_out
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
           save_bss_rds(bss_model@pvalues, bss_model@main_effect, "pvalues", bss_data, bss_model, outdir) # Save pvalues as a rds file
         },
         bss_lm = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@main_effect, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@main_effect, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues_adjusted, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(bss_model@tvalues, bss_model@main_effect, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(bss_model@tvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           save_bss_rds(bss_model@pvalues, bss_model@main_effect, "pvalues", bss_data, bss_model, outdir) # Save pvalues as a rds file
           },
         bss_corr = {
           bss_model@corr_values[abs(log_pvalues) <= -1*log10(0.05)] <- 0
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@corr_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@corr_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues_adjusted, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(bss_model@corr_values, bss_model@corr_var, "corr_values", bss_data, bss_model, outdir)
           save_bss_out_surface(bss_model@corr_values, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           save_bss_rds(bss_model@pvalues, bss_model@corr_var, "pvalues", bss_data, bss_model, outdir) # Save pvalues as a rds file
           },
         pairedttest = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@group_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@group_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues_adjusted, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(bss_model@tvalues, bss_model@group_var, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(bss_model@tvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           save_bss_rds(bss_model@pvalues, bss_model@group_var, "pvalues", bss_data, bss_model, outdir) # Save pvalues as a rds file
         },
         unpairedttest = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@group_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@group_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_surface(log_pvalues_adjusted, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(bss_model@tvalues, bss_model@group_var, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_surface(bss_model@tvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           save_bss_rds(bss_model@pvalues, bss_model@group_var, "pvalues", bss_data, bss_model, outdir) # Save pvalues as a rds file
         }
  )

  # Copy modelspec file to the output directory
  file.copy(bss_model@mspec_file, bss_out@outdir)
  return(bss_out)
  }
)

#' @rdname save_out
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
           save_bss_rds(bss_model@pvalues, bss_model@main_effect, "pvalues", bss_data, bss_model, outdir) # Save pvalues as a rds file
         },

         bss_lm = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@main_effect, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@main_effect, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues_adjusted, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(tvalues, bss_model@main_effect, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(tvalues, bss_model@main_effect, bss_cmap, bss_data, bss_model, outdir)
           save_bss_rds(bss_model@pvalues, bss_model@main_effect, "pvalues", bss_data, bss_model, outdir) # Save pvalues as a rds file
         },
         bss_corr = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@corr_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@corr_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues_adjusted, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(corr_values, bss_model@corr_var, "corr_values", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(corr_values, bss_model@corr_var, bss_cmap, bss_data, bss_model, outdir)
           save_bss_rds(bss_model@pvalues, bss_model@corr_var, "pvalues", bss_data, bss_model, outdir) # Save pvalues as a rds file
         },
         pairedttest = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@group_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@group_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues_adjusted, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(tvalues, bss_model@group_var, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(tvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           save_bss_rds(bss_model@pvalues, bss_model@group_var, "pvalues", bss_data, bss_model, outdir) # Save pvalues as a rds file
         },
         unpairedttest = {
           bss_cmap <- save_bss_color_files(log_pvalues, bss_model@group_var, "log_pvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(log_pvalues_adjusted, bss_model@group_var, "log_pvalues_adjusted", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(log_pvalues_adjusted, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           bss_cmap <- save_bss_color_files(tvalues, bss_model@group_var, "tvalues", bss_data, bss_model, outdir)
           save_bss_out_nifti_image(tvalues, bss_model@group_var, bss_cmap, bss_data, bss_model, outdir)
           save_bss_rds(bss_model@pvalues, bss_model@group_var, "pvalues", bss_data, bss_model, outdir) # Save pvalues as a rds file
         }
  )

    # Copy modelspec file to the output directory
  file.copy(bss_model@mspec_file, bss_out@outdir)
  invisible(bss_out)
  }
)


#' @rdname save_out
setMethod("save_out", valueClass = "BssDBMOutput", signature = "BssDBMOutput", function(bss_out, bss_data, bss_model) {

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



#' @rdname save_out
setMethod("save_out", valueClass = "BssROIOutput", signature = "BssROIOutput", function(bss_out, bss_data, bss_model) {

  # Get the absolute path of outdir
  outdir <- tools::file_path_as_absolute(bss_out@outdir)

  # Copy demographics csv file to output directory
  csvfilename <- file.path(bss_out@outdir, basename(bss_data@csv))
  write.csv(bss_data@demographics, csvfilename)

  nb_header <- "### BrainSuite ROI statistical analysis report"

  nb_libraries <-"```{r librar_cmds, echo=FALSE}\n"
  nb_libraries <- paste(nb_libraries, "\nlibrary('bssr')\nlibrary('ggplot2')\n", sep="")
  nb_libraries <- paste(nb_libraries, "\n```\n", sep="")

  nb_data_header_one <- "The following command loads the data."

  # Add commands to load the data
  nb_data_command_one <- sprintf("\n```{r message=FALSE, warning=FALSE, results='hide', data_command_one}\n")
  nb_data_command_one <- paste(nb_data_command_one,"\n",
                            bss_data@load_data_command,"\n",sep = "")
  nb_data_command_one <- paste(nb_data_command_one, "\n```\n", sep = "")


  nb_load_data <- sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, load_data}\n")
  nb_load_data <- paste(nb_load_data, "\nDT::datatable(bss_data@demographics[,-(ncol(bss_data@demographics))], rownames = FALSE)\n", sep = "")
  nb_load_data <- paste(nb_load_data, "\n```\n", sep = "")


  nb_data_header_two <- "The following command creates the model used to analyze the data."

  nb_data_command_two <- sprintf("\n```{r message=FALSE, warning=FALSE, results='hide', data_command_two}\n")
  nb_data_command_two <- paste(nb_data_command_two, "\n",bss_model@load_data_command,"\n")
  nb_data_command_two <- paste(nb_data_command_two, "\n```\n", sep = "")

  nb_data_header_three <- "The final command (below) was used to render this document."

  command_save_bss_out <- sprintf("save_bss_out(bss_data, bss_model, outdir = '%s')", bss_out@outdir)

  nb_data_command_three <-  sprintf("\n```{r eval=FALSE, message=FALSE, warning=FALSE, data_command_three}\n")
  nb_data_command_three <- paste(nb_data_command_three, "\n", command_save_bss_out,"\n")
  nb_data_command_three <- paste(nb_data_command_three, "\n```\n", sep = "")


  nb_commands <- lapply(as.character(bss_data@roiids),function(x){sprintf("```{r warning=FALSE, run_command_%s}\n",x)})
  nb_plots <- lapply(as.character(bss_data@roiids),function(x){sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, plot_%s}\n",x)})
  nb_calculations <- lapply(as.character(bss_data@roiids),function(x){sprintf("```{r echo=FALSE, message=FALSE, warning=FALSE, pval_%s}\n",x)})


  selected_col <- rep(NA, length(bss_data@roiids))
  for (i in 1:length(bss_data@roiids)){
    selected_col[i] <- paste0(as.character(get_roi_tag(read_label_desc(),bss_data@roiids[i])), "(",bss_data@roiids[i],")")
    bss_data@demographics[,selected_col[i]]
  }

  for (m in 1:length(bss_data@roiids)){
    for (i in 1:length(bss_model@stats_commands[[m]])) {
      nb_commands[[m]] <- paste0(nb_commands[[m]], bss_model@stats_commands[[m]][i], "\n")
    }
    nb_commands[[m]] <- paste0(nb_commands[[m]], "```\n\n")
    nb_calculations[[m]] <- paste0(nb_calculations[[m]],"anova_table <- anova(lm_full_",
                                   bss_data@roiids[m],", lm_null_",bss_data@roiids[m],")\np_val_",
                                   bss_data@roiids[m]," <- round(anova_table$`Pr(>F)`[2],digits=4)\n",
                                   "pval_string <- paste('pvalue:', as.character(p_val_",bss_data@roiids[m],"))\n```\n\n")
    if (class(bss_data@demographics[,gsub("([A-Za-z]+).*", "\\1", bss_model@fullmodel)])=="integer"|class(bss_data@demographics[,gsub("([A-Za-z]+).*", "\\1", bss_model@fullmodel)])=="double"){
      nb_plots[[m]]<-paste0(nb_plots[[m]],"ggplot2::ggplot(data=bss_data@demographics, ggplot2::aes(x=",
                            gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),
                            ", y = `",selected_col[m],"`)) + ggplot2::geom_point() + ggplot2::ggtitle('",
                            as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bss_data@roiids[m])[[1]]),
                            " ", bss_data@roimeas,
                            " vs ", bss_model@main_effect, "') + ggplot2::labs(y='",
                            as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bss_data@roiids[m])[[1]]),
                            "') + ggplot2::geom_smooth(method=lm, se=TRUE) +
ggplot2::xlim(",min(bss_data@demographics[,gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel)]),", ",max(bss_data@demographics[,gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel)])+6,") +
ggplot2::theme(axis.title=ggplot2::element_text(size=16,face='bold')) +
ggplot2::theme(plot.title=ggplot2::element_text(size=18,face='bold')) +
ggplot2::annotate('label',x=",max(bss_data@demographics[,gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel)])+2,",y= max(bss_data@demographics$`",selected_col[m],"`)",
                            ",label= paste('pvalue:', as.character(p_val_",bss_data@roiids[m],")))\n
ggplot2::ggsave(filename='",
                            paste0(as.character(get_roi_tag(read_label_desc(),bss_data@roiids[m]))), "_roi",bss_data@roiids[m],
                            "_", bss_data@roimeas,
                            "_vs_", bss_model@main_effect, ".pdf',device='pdf')\n```\n")
    } else if (class(bss_data@demographics[,gsub("([A-Za-z]+).*", "\\1", bss_model@fullmodel)])=="factor"){
           nb_plots[[m]]<-paste0(nb_plots[[m]],
"mean_lengths <- rep(NA, length(levels(bss_data@demographics$",gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),")))
std_devs <- rep(NA, length(levels(bss_data@demographics$",gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),")))
for (i in 1:length(levels(bss_data@demographics$",gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),"))){
  mean_lengths[i] <-  mean(bss_data@demographics[bss_data@demographics$",gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),
" == levels(bss_data@demographics$",gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),")[i],]$`",
selected_col[m],"`)
  std_devs[i] <- sd (bss_data@demographics[bss_data@demographics$",gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),
" == levels(bss_data@demographics$",gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),")[i],]$`",
selected_col[m],"`)
} \n",
"modified_df <- data.frame(len = mean_lengths, sd = std_devs,",gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel)," = levels(bss_data@demographics$",
                gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),"))\n",
"ggplot2::ggplot(data=modified_df, ggplot2::aes(x = ",gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),
                 ", y = len, color = ",gsub('([A-Za-z]+).*', '\\1', bss_model@fullmodel),
                 ")) +
ggplot2::geom_bar(stat = 'identity') + ggplot2::ggtitle('",
                            as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bss_data@roiids[m])[[1]]),
                            " ", bss_data@roimeas,
                            " vs ", bss_model@main_effect,"') +
ggplot2::labs(y='",as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bss_data@roiids[m])[[1]]),
                            "') +
ggplot2::geom_errorbar(mapping=aes(ymin=len-sd, ymax=len+sd)) +
ggplot2::theme(axis.title=ggplot2::element_text(size=16,face='bold')) +
ggplot2::theme(plot.title=ggplot2::element_text(size=18,face='bold')) +
ggplot2::scale_color_discrete(name= pval_string)\n
ggplot2::ggsave(filename='",
                            paste0(as.character(get_roi_tag(read_label_desc(),bss_data@roiids[m]))), "_roi",bss_data@roiids[m],
                            "_", bss_data@roimeas,
                            "_vs_", bss_model@main_effect, ".pdf',device='pdf')\n```\n")
    }
    nb_commands[[m]] <- paste0(sprintf("\n#### Main effect of %s (%d) %s on %s controlling for %s \n",
                                      as.character(get_roi_name(label_desc_df = read_label_desc(),roiid=bss_data@roiids[m])[[1]]),
                                      bss_data@roiids[m],
                                      bss_data@roimeas, bss_model@main_effect,bss_model@covariates ),
                              nb_commands[[m]])
  }

  rmdfileconn<-file(file.path(outdir, "report.Rmd"))
  writeLines(c(nb_header, nb_libraries, nb_data_header_one, nb_data_command_one,
               nb_load_data, nb_data_header_two, nb_data_command_two,
               nb_data_header_three, nb_data_command_three,
               unlist(mapply(paste0, nb_commands, nb_calculations, nb_plots))), rmdfileconn)
  close(rmdfileconn)

  # Render the markdown
  rmarkdown::render(file.path(outdir, "report.Rmd"), output_file=file.path(outdir, "report.html"), quiet = TRUE)

  # Copy modelspec file to the output directory
  file.copy(bss_model@mspec_file, bss_out@outdir)
  return(bss_out)
  }
)

#' Save the statistical analysis output
#' @param bss_data object of type \code{BssData}
#' @param bss_model object of type \code{BssModel}
#' @param outdir output directory to save the results
#' @export
save_bss_out <- function(bss_data, bss_model, outdir="") {

  valid_types <- c("cbm", "tbm", "roi", "dbm", "nca")
  if (! bss_data@analysis_type %in% valid_types)
    stop(sprintf("Valid data types are %s.", paste(valid_types, collapse = ', ')), call. = FALSE)

  switch(bss_data@analysis_type,
         cbm = { bss_out <- new("BssCBMOutput", outdir)},
         tbm = { bss_out <- new("BssTBMOutput", outdir) },
         dbm = { bss_out <- new("BssDBMOutput", outdir) },
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

  # save ini colormap with ranges
  ini_fileprefix <- paste(paste(bss_model@model_type, var_name, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), bss_cmap@cmap_type, sep = '_'), '.ini', sep = '')
  save_colormap_to_ini(file.path(outdir, ini_fileprefix), bss_cmap)
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

save_bss_rds <- function(measure, var_name, label, bss_data, bss_model, outdir) {

  outprefix <- paste(paste(bss_model@model_type, var_name, tools::file_path_sans_ext(
    basename(bss_data@atlas_filename)), label, sep = '_'), ".rds", sep = '')
  saveRDS(measure, file=file.path(outdir, outprefix))
}
