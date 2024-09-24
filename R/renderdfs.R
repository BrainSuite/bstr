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

renderSurfaceMontage <- function(studybase, studysuffix, surface_output_base, bfc_file,
                                 colorbar_file = "", colorbar_label = "",
                                 image_pixels = 1024) {
  zoom <- 0.45
  border <- floor(image_pixels / 32)
  left <- paste0(studybase, ".left.mid.cortex_", studysuffix, ".dfs")
  right <- paste0(studybase, ".right.mid.cortex_", studysuffix, ".dfs")

  dfs <- paste0("\"",file.path(get_brainsuite_install_path(),"bin/renderdfs"),"\" ")
  callbase <- paste0(dfs," --vol ", bfc_file, " --zoom ", zoom, " -x ", image_pixels, " -y ", image_pixels)

  system_call_output<-system(paste0(callbase, " --left -s ", left, " -o ", surface_output_base, ".left.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  system_call_output<-system(paste0(callbase, " --right -s ", right, " -o ", surface_output_base, ".right.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  system_call_output<-system(paste0(callbase, " --right -s ", left, " -o ", surface_output_base, ".left_mesial.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  system_call_output<-system(paste0(callbase, " --left -s ", right, " -o ", surface_output_base, ".right_mesial.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  system_call_output<-system(paste0(callbase, " --ant -s ", left, " ", right, " -o ", surface_output_base, ".anterior.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)
  system_call_output<-system(paste0(callbase, " --pos -s ", left, " ", right, " -o ", surface_output_base, ".posterior.png"),
    intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE, input = NULL)

  views <- c(".left.png",      ".anterior.png",  ".right.png",
             ".left_mesial.png", ".posterior.png", ".right_mesial.png")
  for (i in 1:6) {
    views[i] <- paste0(surface_output_base, views[i])
  }
  images <- magick::image_read(views)
  images <- image_repage(image_trim(images))
  for (i in 1:6) {
    magick::image_write(images[i], path = views[i])
  }
  montage <- image_transparent(image_montage(images, bg = "transparent", tile="3x2", geometry = paste0("+", border, "+", border)), color = "black")
  if (colorbar_file != "") {
    cbar <- label_color_bar(colorbar_file, colorbar_label)
    montage <- image_append(c(montage, image_resize(cbar,geometry=paste0("x",image_info(montage)["height"]))))
  }
  return(montage)
}
