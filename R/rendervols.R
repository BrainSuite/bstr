# BrainSuite Statistics Toolbox in R (bstr)
# Copyright (C) 2025 The Regents of the University of California
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


#' Render statistical map overlays on BrainSuite atlas
#' @param outdir path for the output directory to save the results
#' @param bstr_data object of type `BstrData`
#' @param bstr_model object of type `BstrModel`
#' @param var_name string specifying the statistical variable
#' @param stats_string list of strings specifying the statistical measure
#' @param voxelcoord list of voxel coordinates (peak cluster value)
#' @param alpha numeric value for render transparency. Default set to 120
#' @export
render_statmap_on_atlas <- function(outdir, bstr_data, bstr_model, var_name, stats_string, voxelcoord, alpha=120) {

  #create a folder to store png images in the output directory
  dir.create(paste0(outdir,"/png_images"))
  dir.create(paste0(outdir,"/png_images_crosshairs"))

  views <- c("sag", "cor", "ax")
  border = 384/32

  # For each cluster/voxelcoord, call statmap
  for (vxl in seq_along(voxelcoord)) {

    # Build a list of statmap commands to be executed on [no. of views (3) x stats_measures (4) ]
    render_cmd <- paste0(
      file.path(get_brainsuite_install_path(),bs_binary_files$statmap),
      " --atlas ",
      bstr_data@atlas_filename,
      " -o ",
      outdir, .Platform$file.sep, "png_images", .Platform$file.sep,
      apply(expand.grid(paste0(views, voxelcoord[[vxl]]), "_", stats_string, sprintf("_cluster%s.png", vxl)), 1, paste0, collapse = ""),
      " --stat ",
      outdir, .Platform$file.sep,
      rep(paste0(bstr_model@model_type, "_", var_name, "_", tools::file_path_sans_ext(
        basename(bstr_data@atlas_filename)), "_", stats_string,  bstr_data@data_type), times=3),
      " --xhaircolor 20 20 20 --xhair ",
      outdir, .Platform$file.sep, "png_images_crosshairs", .Platform$file.sep,
      apply(expand.grid(paste0(views, voxelcoord[[vxl]]), "_", stats_string, sprintf("_cluster%s.png", vxl)), 1, paste0, collapse = ""),
      " -p ",
      paste(as.character(voxelcoord[[vxl]]), collapse=" "),
      paste0(" --", rep(views, times=4)),
      " --isotropic ", " -a ", alpha,
      " --cbar ",
      outdir, .Platform$file.sep,
      rep(paste0(bstr_model@model_type, "_", var_name, "_", tools::file_path_sans_ext(
        basename(bstr_data@atlas_filename)), "_", stats_string,  ".cbar"), times=3)
    )
    for (cmd in render_cmd) {
      system_call_output <- system(cmd,intern=TRUE, ignore.stdout=FALSE, ignore.stderr=TRUE, wait=TRUE, input=NULL)
      if (!is.null(attributes(system_call_output))){
        warning(paste0("Statmap error. Status ",attributes(system_call_output)$status, " returned."),call. = FALSE)
      }
    }

    # Create image moscaics for each cluster
    for (jj in seq_along(stats_string)) {
      view_images <- apply(expand.grid(paste0(outdir, .Platform$file.sep, "png_images_crosshairs", .Platform$file.sep, views, voxelcoord[[vxl]]), "_", stats_string[jj], sprintf("_cluster%s.png", vxl)), 1, paste0, collapse = "")
      images <- magick::image_read(view_images)
      montage <- magick::image_transparent(magick::image_montage(images, bg = "transparent", tile="3x1", geometry = paste0("+", border, "+", border)), color = "black")

      cbar_filename <- paste(outdir, .Platform$file.sep, paste(bstr_model@model_type, var_name, tools::file_path_sans_ext(
        basename(bstr_data@atlas_filename)), stats_string[jj], sep = '_'), '_cbar.png', sep = '')
      cbar <- label_color_bar(cbar_filename, var_name)
      montage <- magick::image_append(c(montage, magick::image_resize(cbar,geometry=paste0("x",magick::image_info(montage)["height"]))))
      magick::image_write(montage, paste0(file.path(outdir,sprintf("cluster%s_", vxl)), stats_string[jj], "_figure.png"),format="png")
    }
  }
}

