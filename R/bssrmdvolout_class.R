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

#' R6 derived class for Rmd volume output functionality
#' @export
BssRmdVolumeOutput <-
  R6::R6Class("BssRmdVolumeOutput",
              inherit = BssRmdOutput,
              public = list(
                initialize = function(outdir = "./") {
                  initialize(outdir)
                },
                save_out = function(bss_data, bss_model, outdir, voxelcoord) {
                  get_custom_overlays = function(outdir) {
                    if (bss_model@model_type=="unpairedttest" | bss_model@model_type=="pairedttest"){
                      indep_var <- bss_model@group_var
                    } else {
                      indep_var <- bss_model@main_effect
                    }

                    adjp_overlay <- paste0(outdir, "/", bss_model@model_type, "_", indep_var,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$log_pvalues_adjusted, bss_data@data_type)
                    adjt_overlay <- paste0(outdir, "/", bss_model@model_type, "_", indep_var,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$tvalues_adjusted, bss_data@data_type)
                    p_overlay <- paste0(outdir, "/", bss_model@model_type, "_", indep_var,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$log_pvalues, bss_data@data_type)
                    t_overlay <- paste0(outdir, "/", bss_model@model_type, "_", indep_var,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$tvalues, bss_data@data_type)
                    corr_overlay <- paste0(outdir, "/", bss_model@model_type, "_", bss_model@corr_var,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$corr_values, bss_data@data_type)
                    adj_corr_overlay <- paste0(outdir, "/", bss_model@model_type, "_", bss_model@corr_var,"_",tools::file_path_sans_ext(
                       basename(bss_data@atlas_filename)),"_", bs_stat_overlays$corr_values_masked_adjusted, bss_data@data_type)

                    return(list("adjp_overlay" = adjp_overlay, "adjt_overlay" = adjt_overlay, "p_overlay" = p_overlay,
                                "t_overlay" = t_overlay, "corr_overlay" = corr_overlay, "adj_corr_overlay" = adj_corr_overlay))
                  }

                  get_custom_lut = function(outdir) {
                    switch(bss_model@model_type,
                           bss_anova = {var_name = bss_model@main_effect},
                           bss_lm = {var_name = bss_model@main_effect},
                           bss_lmer = {var_name = bss_model@main_effect},
                           bss_corr = {var_name = bss_model@corr_var},
                           pairedttest = {var_name = bss_model@group_var},
                           unpairedttest = {var_name = bss_model@group_var}
                    )
                    adjp_lut <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "log_pvalues_adjusted.lut")
                    adjt_lut <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "tvalues_adjusted.lut")
                    p_lut <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "log_pvalues.lut")
                    t_lut <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "tvalues.lut")
                    if (bss_model@model_type == "bss_corr"){
                      corr_lut <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                        basename(bss_data@atlas_filename)),"_", "corr_values.lut")
                      adjcorr_lut <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                        basename(bss_data@atlas_filename)),"_", "corr_values_masked_adjusted.lut")
                      return(list("adjp_lut" = adjp_lut, "adjt_lut" = adjt_lut, "p_lut" = p_lut, "t_lut" = t_lut,"corr_lut" = corr_lut, "adjcorr_lut" = adjcorr_lut))
                    }

                    return(list("adjp_lut" = adjp_lut, "adjt_lut" = adjt_lut, "p_lut" = p_lut, "t_lut" = t_lut))

                  }
                  get_min_vals = function(outdir) {
                    switch(bss_model@model_type,
                           bss_anova = {var_name = bss_model@main_effect},
                           bss_lm = {var_name = bss_model@main_effect},
                           bss_lmer = {var_name = bss_model@main_effect},
                           bss_corr = {var_name = bss_model@corr_var},
                           pairedttest = {var_name = bss_model@group_var},
                           unpairedttest = {var_name = bss_model@group_var}
                    )
                    adjp_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "log_pvalues_adjusted.ini")
                    adjp_min <- as.numeric(ini::read.ini(adjp_ini)$colormap$cnegmax)
                    adjt_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "tvalues_adjusted.ini")
                    adjt_min <- as.numeric(ini::read.ini(adjt_ini)$colormap$cnegmax)
                    p_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "log_pvalues.ini")
                    p_min <- as.numeric(ini::read.ini(p_ini)$colormap$cnegmax)
                    t_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "tvalues.ini")
                    t_min <- as.numeric(ini::read.ini(t_ini)$colormap$cnegmax)
                    if (bss_model@model_type == "bss_corr"){
                      corr_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                        basename(bss_data@atlas_filename)),"_", "corr_values.ini")
                      corr_min <- as.numeric(ini::read.ini(corr_ini)$colormap$cnegmax)
                      adjcorr_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                        basename(bss_data@atlas_filename)),"_", "corr_values_masked_adjusted.ini")
                      adjcorr_min <- as.numeric(ini::read.ini(adjcorr_ini)$colormap$cnegmax)
                      return(c(adjp_min,adjt_min,p_min,t_min,corr_min,adjcorr_min))
                    }
                    return(c(adjp_min,adjt_min,p_min,t_min))
                  }
                  get_max_vals = function(outdir) {
                    switch(bss_model@model_type,
                           bss_anova = {var_name = bss_model@main_effect},
                           bss_lm = {var_name = bss_model@main_effect},
                           bss_lmer = {var_name = bss_model@main_effect},
                           bss_corr = {var_name = bss_model@corr_var},
                           pairedttest = {var_name = bss_model@group_var},
                           unpairedttest = {var_name = bss_model@group_var}
                    )
                    adjp_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "log_pvalues_adjusted.ini")
                    adjp_max <- as.numeric(ini::read.ini(adjp_ini)$colormap$cposmax)
                    adjt_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "tvalues_adjusted.ini")
                    adjt_max <- as.numeric(ini::read.ini(adjt_ini)$colormap$cposmax)
                    p_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "log_pvalues.ini")
                    p_max <- as.numeric(ini::read.ini(p_ini)$colormap$cposmax)
                    t_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", "tvalues.ini")
                    t_max <- as.numeric(ini::read.ini(t_ini)$colormap$cposmax)
                    if (bss_model@model_type == "bss_corr"){
                      corr_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                        basename(bss_data@atlas_filename)),"_", "corr_values.ini")
                      corr_max <- as.numeric(ini::read.ini(corr_ini)$colormap$cposmax)
                      adjcorr_ini <- paste0(outdir, "/", bss_model@model_type, "_", var_name,"_",tools::file_path_sans_ext(
                        basename(bss_data@atlas_filename)),"_", "corr_values_masked_adjusted.ini")
                      adjcorr_max <- as.numeric(ini::read.ini(adjcorr_ini)$colormap$cposmax)
                      return(c(adjp_max,adjt_max,p_max,t_max,corr_max,adjcorr_max))
                    }
                    return(c(adjp_max,adjt_max,p_max,t_max))
                  }


                  #create a folder to store png images in
                  dir.create(paste0(outdir,"/png_images"))
                  dir.create(paste0(outdir,"/png_images_crosshairs"))
                  for(cluster_iter in 1:length(voxelcoord)) {
                    if (bss_model@model_type=="bss_corr"){
                      private$render_overlay(
                        bss_data,
                        bss_model,
                        cluster_iter,
                        voxelcoord,
                        atlaspath = bss_data@atlas_filename,
                        overlaypath = c(get_custom_overlays(outdir)[[5]],get_custom_overlays(outdir)[[6]]),
                        #stat_overlay
                        lutpaths = c(get_custom_lut(outdir)[[1]],get_custom_lut(outdir)[[2]],get_custom_lut(outdir)[[3]],get_custom_lut(outdir)[[4]],get_custom_lut(outdir)[[5]],get_custom_lut(outdir)[[6]]),
                        max_vals = c(get_max_vals(outdir)[1],get_max_vals(outdir)[2],get_max_vals(outdir)[3],get_max_vals(outdir)[4],get_max_vals(outdir)[5],get_max_vals(outdir)[6]),
                        min_vals = c(get_min_vals(outdir)[1],get_min_vals(outdir)[2],get_min_vals(outdir)[3],get_min_vals(outdir)[4],get_min_vals(outdir)[5],get_min_vals(outdir)[6]),
                        outdir,
                        name = c(bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues_adjusted,bs_stat_overlays$log_pvalues,bs_stat_overlays$tvalues,bs_stat_overlays$corr_values,bs_stat_overlays$corr_values_masked_adjusted), alpha = 120)
                    } else {
                    private$render_overlay(
                      bss_data,
                      bss_model,
                      cluster_iter,
                      voxelcoord,
                      atlaspath = bss_data@atlas_filename,
                      overlaypath = c(get_custom_overlays(outdir)[[1]],get_custom_overlays(outdir)[[2]],get_custom_overlays(outdir)[[3]],get_custom_overlays(outdir)[[4]]),
                      #stat_overlay
                      lutpaths = c(get_custom_lut(outdir)[[1]],get_custom_lut(outdir)[[2]],get_custom_lut(outdir)[[3]],get_custom_lut(outdir)[[4]]),
                      max_vals = c(get_max_vals(outdir)[1],get_max_vals(outdir)[2],get_max_vals(outdir)[3],get_max_vals(outdir)[4]),
                      min_vals = c(get_min_vals(outdir)[1],get_min_vals(outdir)[2],get_min_vals(outdir)[3],get_min_vals(outdir)[4]),
                      outdir,
                      name = c(bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues_adjusted,bs_stat_overlays$log_pvalues,bs_stat_overlays$tvalues), alpha = 120)
                    }
                    private$render_atlas(cluster_iter, voxelcoord,
                                         atlaspath = bss_data@atlas_filename,
                                         outdir)
                    # }

                    private$render_html(bss_data,
                                        bss_model,
                                        outdir,
                                        voxelcoord,
                                        overlay_name = c(bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues_adjusted,bs_stat_overlays$log_pvalues,bs_stat_overlays$tvalues,bs_stat_overlays$corr_values,bs_stat_overlays$corr_values_masked_adjusted))


                  }

                  private$save_rmd_preamble(bss_data, bss_model, file.path(outdir,  sprintf("report_%s_%s.Rmd", bss_model@model_type, switch(bss_model@model_type,
                                                                                                                        bss_corr = {bss_model@corr_var},
                                                                                                                        unpairedttest = {bss_model@group_var},
                                                                                                                        pairedttest = {bss_model@group_var},
                                                                                                                        bss_anova = {bss_model@main_effect},
                                                                                                                        bss_lm = {bss_model@main_effect},
                                                                                                                        bss_lmer = {bss_model@main_effect}))),
                                            outdir,
                                            voxelcoord,
                                            overlay_name = "c(bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues_adjusted,bs_stat_overlays$log_pvalues,bs_stat_overlays$tvalues,bs_stat_overlays$corr_values)")

                  rmarkdown::render(file.path(outdir,  sprintf("report_%s_%s.Rmd", bss_model@model_type, switch(bss_model@model_type,
                                                                                                                bss_corr = {bss_model@corr_var},
                                                                                                                unpairedttest = {bss_model@group_var},
                                                                                                                pairedttest = {bss_model@group_var},
                                                                                                                bss_anova = {bss_model@main_effect},
                                                                                                                bss_lm = {bss_model@main_effect},
                                                                                                                bss_lmer = {bss_model@main_effect}))))
                }
              ),

              private = list(

                render_overlay = function(bss_data, bss_model, voxelcoord_index,voxelcoord,atlaspath,overlaypath,outdir,name,alpha,lutpaths,min_vals,max_vals) {
                  view_order <- c("sag","cor","ax")
                  switch(bss_model@model_type,
                         bss_anova = {var_name = bss_model@main_effect},
                         bss_lm = {var_name = bss_model@main_effect},
                         bss_lmer = {var_name = bss_model@main_effect},
                         bss_corr = {var_name = bss_model@corr_var},
                         pairedttest = {var_name = bss_model@group_var},
                         unpairedttest = {var_name = bss_model@group_var}
                  )
                  # for (lut_index in 1:length(lutpaths)){
                  #   json_file <- paste0("lut2cbar -i ",lutpaths[lut_index]," -o ",outdir,"/",bss_model@model_type,"_",var_name,"_mri.bfc.nii_",name[lut_index],".cbar",
                  #                       " --min ",min_vals[lut_index], " --max ",max_vals[lut_index])
                  #   system(json_file,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                  # }
                  for (stats_measure_index in 1:length(overlaypath)) {
                    #if (min_vals[stats_measure_index]==0 & max_vals[stats_measure_index]==0){next}
                    #crosshair_length_multiplier <- 15
                    for (view in 1:3){
                       name_index <- ifelse(bss_model@model_type == "bss_corr",stats_measure_index+4, stats_measure_index)
                       # current_view_with_crosshairs <- paste0("statmap -i ",overlaypath[stats_measure_index], " -o ",outdir,"/png_images_crosshairs/",view_order[view],voxelcoord[[voxelcoord_index]][view], "_",name[name_index], "_cluster",voxelcoord_index,
                       #                                        ".png --atlas ", atlaspath," --slice ", voxelcoord[[voxelcoord_index]][view], " --", view_order[view],
                       #                                        " --max ", max_vals[name_index], " --min ", min_vals[name_index],
                       #                                        " -a ", alpha, " --cbar ", outdir,"/",bss_model@model_type,"_",var_name,"_", tools::file_path_sans_ext(basename(bss_data@atlas_filename)),"_",name[name_index],".cbar --isotropic")
                       # system(current_view_with_crosshairs,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                       current_view <- paste0("\"",file.path(get_brainsuite_install_path(),bs_binary_files$statmap),"\" --atlas \"",atlaspath, "\" -o \"",outdir,"/png_images/",view_order[view],voxelcoord[[voxelcoord_index]][view], "_",name[name_index], "_cluster",voxelcoord_index,
                                              ".png\" --stat \"", overlaypath[stats_measure_index], "\" --xhair \"",outdir,"/png_images_crosshairs/",view_order[view],voxelcoord[[voxelcoord_index]][view], "_",name[name_index], "_cluster",voxelcoord_index,
                                              ".png\" -p ", voxelcoord[[voxelcoord_index]][1]," ",voxelcoord[[voxelcoord_index]][2]," ",voxelcoord[[voxelcoord_index]][3]," --", view_order[view], " --isotropic",
                                              " -a ", alpha, " --cbar \"", outdir,"/",bss_model@model_type,"_",var_name,"_", tools::file_path_sans_ext(basename(bss_data@atlas_filename)),"_",name[name_index],".cbar\"")
                       system_call_output <- system(current_view,intern=TRUE, ignore.stdout=FALSE, ignore.stderr=TRUE, wait=TRUE, input=NULL)
                       if (!is.null(attributes(system_call_output))){
                         warning(paste0("Statmap error. Status ",attributes(system_call_output)$status, " returned."),call. = FALSE)
                       }
                       # coord_range <- c(dim(bss_data@atlas_image)[1],dim(bss_data@atlas_image)[2],dim(bss_data@atlas_image)[3])
                       # if (view == 1){
                       #   x_end = coord_range[2]/coord_range[1]
                       #   y_end = coord_range[3]/coord_range[1]
                       #   x0_center = x_end*(1-(voxelcoord[[voxelcoord_index]][2]*(1/coord_range[2])))
                       #   x_additive = x_end*(crosshair_length_multiplier/coord_range[2])
                       #   y0_center = y_end*(voxelcoord[[voxelcoord_index]][3]*(1/coord_range[3]))
                       #   y_additive = y_end*(crosshair_length_multiplier/coord_range[3])
                       # } else if (view == 2) {
                       #   x_end = coord_range[1]/coord_range[1]
                       #   y_end = coord_range[3]/coord_range[1]
                       #   x0_center = x_end*(1-(voxelcoord[[voxelcoord_index]][1]*(1/coord_range[1])))
                       #   x_additive = x_end*(crosshair_length_multiplier/coord_range[1])
                       #   y0_center = y_end*(voxelcoord[[voxelcoord_index]][3]*(1/coord_range[3]))
                       #   y_additive = y_end*(crosshair_length_multiplier/coord_range[3])
                       # } else {
                       #   x_end = coord_range[1]/coord_range[1]
                       #   y_end = coord_range[2]/coord_range[1]
                       #   x0_center = x_end*(1-(voxelcoord[[voxelcoord_index]][1]*(1/coord_range[1])))
                       #   x_additive = x_end*(crosshair_length_multiplier/coord_range[1])
                       #   y0_center = y_end*(voxelcoord[[voxelcoord_index]][2]*(1/coord_range[2]))
                       #   y_additive = y_end*(crosshair_length_multiplier/coord_range[2])
                       # }
                       # temp_image <- png::readPNG(get_render_image_filename(outdir,voxelcoord,name[name_index], view, voxelcoord_index))
                       # png(get_render_image_filename(outdir,voxelcoord,name[name_index], view, voxelcoord_index))
                       # plot(0:1, 0:1, type='n', axes = F, ann = F)
                       # par(mar = c(0,0,0,0))
                       # rasterImage(temp_image, 0, 0, x_end, y_end)
                       # # Add horizontal crosshair
                       # segments(x0 = x0_center - x_additive, y0 = y0_center, x1 = x0_center + x_additive , y1 = y0_center, col="black", lwd=3)
                       # # Add vertical crosshair
                       # segments(x0 = x0_center, y0 = y0_center - y_additive, x1 = x0_center, y1 = y0_center + y_additive, col="black", lwd=3)
                       # dev.off()
                    }
                  }
                  return(0)
                },
                render_atlas = function(voxelcoord_index,voxelcoord,atlaspath,outdir) {
                  view_name <- c("ax","cor","sag")
                  view_order <- c(3,2,1)
                  for (view_iter in 1:3){
                    current_view <- paste0("\"",file.path(get_brainsuite_install_path(),bs_binary_files$statmap),"\" --atlas \"",atlaspath,"\" -o \"", outdir,"/png_images/",view_name[view_iter],voxelcoord[[voxelcoord_index]][view_order[view_iter]],
                                           "_atlas.png\""," --xhair \"",outdir,"/png_images_crosshairs/",view_name[view_iter],voxelcoord[[voxelcoord_index]][view_order[view_iter]],
                                           ".png\" -p ", voxelcoord[[voxelcoord_index]][1]," ",voxelcoord[[voxelcoord_index]][2]," ",voxelcoord[[voxelcoord_index]][3]," --", view_name[view_iter], " --isotropic")
                    system(current_view,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                  }
                  return(0)
                },

                save_rmd_preamble = function(bss_data, bss_model, rmdfile, outdir, voxelcoord, overlay_name) {

                  load_library <- "library(bssr)"
                  data_command_1 <- paste0("bss_data <- load_bss_data(type = '",bss_data@analysis_type,"', subjdir = '",bss_data@subjdir,"', csv = '", bss_data@csv,"', measure = '", bss_data@measure,"', smooth = ",bss_data@smooth,")")
                  data_command_2 <- ifelse(bss_model@model_type=="bss_corr",
                                           paste0("bss_model <- ",bss_model@model_type,"(corr_var = '",bss_model@corr_var,"', bss_data = bss_data, mult_comp = '",bss_model@mult_comp,"')"),
                                           paste0("bss_model <- ",bss_model@model_type,"(main_effect = '",bss_model@main_effect,"', covariates = '", bss_model@covariates,"', bss_data = bss_data)"))
                  data_command_3 <- "bs_stat_overlays = list(log_pvalues_adjusted = 'log_pvalues_adjusted', tvalues_adjusted = 'tvalues_adjusted', log_pvalues = 'log_pvalues', tvalues = 'tvalues', pvalues = 'pvalues', corr_values = 'corr_values')"
                  user_input <- paste0(load_library,"\n",data_command_1,"\n",data_command_2,"\n",data_command_3,"\n")
                  templines <- deparse(private$render_html)
                  templines[1] <- "render_html = function(bss_data, bss_model, outdir, voxelcoord, overlay_name)"
                  sink(rmdfile, append=TRUE, type = "output")
                  cat("---\n")
                  cat("title: bssr report\n")
                  cat("output: html_document\n")
                  cat("runtime: shiny\n")
                  cat("---\n")
                  cat("<style>\n\n")
                  cat("table, td, th {\n")
                  cat("min-width: 50%;\n")
                  cat("margin-left: auto;\n")
                  cat("margin-right: auto;\n")
                  cat("margin-top: 1em;\n")
                  cat("margin-bottom: 1em;\n")
                  cat("}\n\n")
                  cat("</style>\n")
                  cat("```{r eval=FALSE, echo=FALSE, message=FALSE, results='hide', user_input_commands}\n")
                  cat("#Change the above argument eval to TRUE to knit the file.\n")
                  writeLines(user_input)
                  cat("```\n")
                  cat("```{r echo=FALSE, warning=FALSE}\n")
                  cat("vox_table <- read.table('", outdir,"/cluster.tsv', header = F, sep = '\t')\n", sep = "")
                  cat("colnames(vox_table) <- c('Cluster Number', 'Number of Voxels', ",ifelse(bss_model@model_type=='bss_corr',"'Corr Value'","'T-Value'"),", 'X Coord', 'Y Coord', 'Z Coord')\n")
                  cat("DT::datatable(vox_table, rownames = FALSE)\n")
                  writeLines(templines)
                  voxelcoord_char <- "list("
                  for (individual_voxelcoord in 1:length(voxelcoord)){
                    voxelcoord_char <- paste0(voxelcoord_char, "c(")
                    for (view in 1:length(voxelcoord[[individual_voxelcoord]])){
                      voxelcoord_char <- paste0(voxelcoord_char, voxelcoord[[individual_voxelcoord]][view],",")
                    }
                    voxelcoord_char <- paste0(substr(voxelcoord_char,1,nchar(voxelcoord_char)-1),"),")
                  }
                  voxelcoord_char <- paste0(substr(voxelcoord_char,1,nchar(voxelcoord_char)-1),")")
                  cat("render_html( bss_data, bss_model,  '",outdir,"', ", voxelcoord_char, ", ",overlay_name,")\n")
                  cat("```\n")
                  sink()

                },

                ## another function will generate the names for the pngs

                render_html = function(bss_data, bss_model, outdir, voxelcoord, overlay_name) {
                  overlay <- c(bs_stat_overlays$log_pvalues_adjusted, bs_stat_overlays$tvalues_adjusted, bs_stat_overlays$log_pvalues, bs_stat_overlays$tvalues, bs_stat_overlays$corr_values_masked_adjusted, bs_stat_overlays$corr_values)
                  if (bss_model@model_type=="bss_corr"){
                    cbar <- vector("list",2)
                    for (cbar_index in 1:length(cbar)){
                      cbar[[cbar_index]] <- paste0(paste(bss_model@model_type, bss_model@corr_var, tools::file_path_sans_ext(basename(bss_data@atlas_filename)), overlay[4+cbar_index], sep = '_'), '_cbar.png')
                    }
                  } else if (bss_model@model_type == "pairedttest" || bss_model@model_type == "unpairedttest"){
                    cbar <- vector("list",4)
                    for (cbar_index in 1:length(cbar)){
                      cbar[[cbar_index]] <- paste0(paste(bss_model@model_type, bss_model@group_var, tools::file_path_sans_ext(basename(bss_data@atlas_filename)), overlay[cbar_index], sep = '_'), '_cbar.png')
                    }
                  } else {
                    cbar <- vector("list",4)
                    for (cbar_index in 1:length(cbar)){
                      cbar[[cbar_index]] <- paste0(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(basename(bss_data@atlas_filename)), overlay[cbar_index], sep = '_'), '_cbar.png')
                    }
                  }
                  cluster_filepath <- paste0(outdir,"/cluster.tsv")
                  cluster_filepath <- gsub(" ", "", cluster_filepath, fixed = TRUE)
                  vox_table <- read.table(cluster_filepath, header = F, sep = '\t')
                  colnames(vox_table) <- c('Cluster Number', 'Number of Voxels', 'Maximum Value', 'X Coord', 'Y Coord', 'Z Coord')


                  #function to make rmd work
                  get_render_image_filename <- function(outdir, voxelcoord, overlay_name, brain_sector_index, voxelcoord_index) {
                    view_order <- c("sag","cor","ax")
                    return(paste0("./png_images_crosshairs/", view_order[brain_sector_index], voxelcoord[[voxelcoord_index]][brain_sector_index],"_",overlay_name,"_cluster",voxelcoord_index,".png"))
                  }

                  # Function to return a shiny image object
                  shiny_image = function(outdir, voxelcoord, index_brain_sector, overlay_type, width, voxelcoord_index){
                    align <- "left"
                    return(paste0("shiny::img(src='",get_render_image_filename(outdir, voxelcoord,overlay_type, index_brain_sector, voxelcoord_index),"', align = '", align, "', width = '", width,"')"))
                  }

                  # Function that creates a tabPanel
                  tab_panel = function(panel_type, voxelcoord_index){
                    dims_mm <- dim(atlas_image)[1:3]*RNifti::pixdim(atlas_image)[1:3]
                    # output PNG dimensions for each axis -- this is the dimension's length in mm divided by the smallest voxel edge
                    # this scales the image so that the width/height of one pixel is equal to the smallest dimension of a voxel
                    # so that it doesn't shrink to < 1 pixel, i.e., anisotropic voxels get stretched across two or more pixels rather
                    # than compressed
                    scaledPNGdimensions<-round(dims_mm/min(RNifti::pixdim(atlas_image)[1:3]))
                    # the rounding when the png changes the aspect ratio slightly, so this uses the new aspect ratio
                    # to compute the relative widths based on the height of the images so that images with the same height
                    # are displayed that way. scaling based on the width alone would not maintain that.
                    # so we take the physical height of the image and multiply it by the aspect ratio of the png
                    # to compute the physical width of the resampled image
                    sagWidth <- dims_mm[3] * scaledPNGdimensions[2]/scaledPNGdimensions[3];
                    corWidth <- dims_mm[3] * scaledPNGdimensions[1]/scaledPNGdimensions[3];
                    axWidth  <- dims_mm[2] * scaledPNGdimensions[1]/scaledPNGdimensions[2];
                    # we then rescale this physical width relative to the total width and multiply by 93 to leave 7% space for the colorbar
                    totalWidth <- sagWidth + corWidth + axWidth
                    sagWidth <- 93 * sagWidth / totalWidth
                    corWidth <- 93 * corWidth / totalWidth
                    axWidth <- 93 * axWidth / totalWidth
                    # This leaves 7% for the color bar.
                    width <- c(paste0(sagWidth,"%"),paste0(corWidth,"%"),paste0(axWidth,"%"),"7%")
                    
                    if (bss_model@model_type=="bss_corr"){
                      panel_names<- c("Adjusted Correlation Values","Correlation Values")
                      overlay <- c(bs_stat_overlays$corr_values_masked_adjusted,bs_stat_overlays$corr_values)
                    } else {
                      panel_names <- c("Adjusted P-Values","Adjusted T-Values","P-Values","T-Values")
                      overlay <- c(bs_stat_overlays$log_pvalues_adjusted, bs_stat_overlays$tvalues_adjusted,
                                   bs_stat_overlays$log_pvalues, bs_stat_overlays$tvalues)
                    }
                    images <- ""
                    for (inner_coord_index in 1:3) {
                      images <- paste0(images, shiny_image(outdir, voxelcoord, inner_coord_index,overlay[panel_type], width[inner_coord_index], voxelcoord_index), ",")
                    }
                    images <- paste0(images, "shiny::img(src=paste0('", cbar[panel_type], "'), align='left', width = '", width[4], "'),")
                    images <- substr(images,1,nchar(images)-1)
                    return(paste0("shiny::tabPanel(title_0 = '", panel_names[panel_type], "', value = c('", panel_names[panel_type], "'), shiny::p('", panel_names[panel_type], "'), ",images,")"))
                  }

                  # Function that creates clusters for each panel
                  cluster_panels = function(voxelcoord){
                    if (bss_model@model_type=="bss_corr"){
                      panel_names<- c("Adjusted Correlation Values","Correlation Values")
                    } else {
                      panel_names <- c("Adjusted P-Values","Adjusted T-Values","P-Values","T-Values")
                    }
                    panels <- paste0("shiny::tabPanel(title = shiny::h4(paste('Cluster 1: Voxel Coordinate (",voxelcoord[[1]][1],",",voxelcoord[[1]][2],",",voxelcoord[[1]][3],")'),shiny::tableOutput('data'),shiny::tabsetPanel(id = 'navbar',type = 'pills',")
                    for (tab_panel_type in 1:length(panel_names)){
                      panels <- paste0(panels,tab_panel(tab_panel_type, 1),",")
                    }
                    panels <- paste0(substr(panels,1,nchar(panels)-1),")))")
                    if (length(voxelcoord)>1){
                    for (voxelcoord_index in 2:length(voxelcoord)){
                      panels <- paste0(panels, ", shiny::tabPanel(title = shiny::h4(paste('Cluster ", voxelcoord_index, ": Voxel Coordinate (",voxelcoord[[voxelcoord_index]][1],",",voxelcoord[[voxelcoord_index]][2],",",voxelcoord[[voxelcoord_index]][3],")'),shiny::tableOutput('data'),shiny::tabsetPanel(id = 'navbar',type = 'pills',")
                      for (tab_panel_type in 1:length(panel_names)){
                        panels <- paste0(panels,tab_panel(tab_panel_type, voxelcoord_index),",")
                      }
                      panels <- paste0(substr(panels,1,nchar(panels)-1),")))")
                    }
                    }
                    return(substr(panels,1,nchar(panels)))
                  }


                  eval(parse(text= paste0("shiny::shinyUI(shiny::fluidPage(shinyjs::useShinyjs(),shiny::h3('Choose a cluster and an overlay below'),shiny::tabsetPanel(id = 'navbar',type = 'tabs',",
                        cluster_panels(voxelcoord),")))")))

                }

              )
  )

