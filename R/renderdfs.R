# BrainSuite Statistics Toolbox in R (bstr)
# Copyright (C) 2024 The Regents of the University of California
# renderdfs created by David W. Shattuck, Department of Neurology, Ahmanson Lovelace Brain Mapping Center, UCLA
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


#' Add text label and border to existing colorbar image
#'
#'
#' @param colorbar_file Name of file.
#' @param label Text label that goes under the colorbar (assumed to be two lines for scaling)
#'
#' @export
label_color_bar <- function(colorbar_file, label) {
  cbar <- image_read(colorbar_file)
  h <- image_info(cbar)["height"]
  cbar <- cbar |>
          image_transparent(color="white") |>
          image_trim() |>
          image_repage() |>
          image_border(color = "transparent", geometry = paste0("0x",floor(0.125*h)))
  h2 <- image_info(cbar)["height"]
  w <- image_info(cbar)["width"]
  cbar <- cbar |> image_extent(geometry = paste0(w, "x", h2 + floor(0.125*h)), gravity="North") |>
          image_border(color = "transparent", geometry = paste0(floor(w/5), "x0")) |>
          image_annotate(label, font="helvetica", gravity="South", location="+0+50", size=180)
  return (cbar)
}

#' Create a figure from six views of a pair of hemisphere surfaces (version for bstr).
#'
#'
#' @param surface_output_base path and prefix for output surface pngs
#' @param left_hemi_file filename for right hemisphere
#' @param right_hemi_file filename for right hemisphere
#' @param atlas_image a volumetric image from the atlas, e.g., the bfc file.
#' @param colorbar_file the colorbar file (png) output by the surface output program (optional)
#' @param colorbar_label text that will be displayed beneath the colorbar (optional)
#' @param image_pixels image height for initial renderings of dfs files, which will be cropped after rendering.
#'
#' @export
renderSurfaceFigure <- function(surface_output_base,
                                 left_hemi_file, right_hemi_file, atlas_image,
                                 colorbar_file = "", colorbar_label = "",
                                 image_pixels = 512) {

  renderdfs <- paste0("\"",file.path(get_brainsuite_install_path(),bs_binary_files$renderdfs),"\" ")
  # exit_code = suppressWarnings(system(renderdfs))
  # if (exit_code != 0 || exit_code != NULL)
  if (!file.exists(renderdfs))
  {
    warning("renderdfs is not part of your BrainSuite installation -- please visit https://brainsuite.org/bstr for information on how to obtain this program.")
    return(NULL);
  }
  zoom <- 0.45
  border <- floor(image_pixels / 32)
  callbase <- paste0(renderdfs," --vol ", atlas_image, " --zoom ", zoom, " -x ", image_pixels, " -y ", image_pixels)
  system_call_output<-system(paste0(callbase, " --left -s ", left_hemi_file, " -o ", surface_output_base, "_left.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  system_call_output<-system(paste0(callbase, " --right -s ", right_hemi_file, " -o ", surface_output_base, "_right.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  system_call_output<-system(paste0(callbase, " --right -s ", left_hemi_file, " -o ", surface_output_base, "_left_medial.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  system_call_output<-system(paste0(callbase, " --left -s ", right_hemi_file, " -o ", surface_output_base, "_right_medial.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  system_call_output<-system(paste0(callbase, " --ant -s ", left_hemi_file, " ", right_hemi_file, " -o ", surface_output_base, "_anterior.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  system_call_output<-system(paste0(callbase, " --pos -s ", left_hemi_file, " ", right_hemi_file, " -o ", surface_output_base, "_posterior.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)

  views <- c("_left.png",        "_anterior.png",  "_right.png",
             "_left_medial.png", "_posterior.png", "_right_medial.png")
  for (i in 1:6) {
    views[i] <- paste0(surface_output_base, views[i])
  }
  images <- magick::image_read(views)
  images <- magick::image_repage(magick::image_trim(images))
  for (i in 1:6) {
    magick::image_write(images[i], path = views[i])
  }
  montage <- magick::image_transparent(magick::image_montage(images, bg = "transparent", tile="3x2", geometry = paste0("+", border, "+", border)), color = "black")
  if (colorbar_file != "") {
    cbar <- label_color_bar(colorbar_file, colorbar_label)
    montage <- magick::image_append(c(montage, magick::image_resize(cbar,geometry=paste0("x",magick::image_info(montage)["height"]))))
  }
  return(montage)
}
