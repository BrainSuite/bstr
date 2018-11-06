# BrainSuite Statistics Toolbox in R (bssr)
# Copyright (C) 2017 The Regents of the University of California
# Creator: Shantanu H. Joshi,e Department of Neurology, Ahmanson Lovelace Brain Mapping Center, UCLA
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
              #inherit = BssRmdOutput,
              public = list(
                initialize = function(outdir = "./") {
                  initialize(outdir)
                },
                save_out = function(bss_data, bss_model, outdir, voxelcoord) {
                  get_custom_tbm_overlays = function(outdir) {

                    p_overlay <- paste(outdir, bss_model@model_type, "_", bss_model@main_effect,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$log_pvalues, bss_data@data_type, sep="")
                    adjp_overlay <- paste(outdir, bss_model@model_type, "_", bss_model@main_effect,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$log_pvalues_adjusted, bss_data@data_type, sep="")
                    t_overlay <- paste(outdir, bss_model@model_type, "_", bss_model@main_effect,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$tvalues, bss_data@data_type, sep="")

                    return(list("p_overlay" = p_overlay, "adjp_overlay" = adjp_overlay, "t_overlay" = t_overlay))
                  }


                  #create a folder to store png images in
                  dir.create(paste0(outdir,"PNG_images"))
                  for(cluster_iter in 1:length(voxelcoord)) {
                    private$render_overlay(
                      cluster_iter,
                      voxelcoord,
                      atlaspath = bss_data@atlas_filename,
                      overlaypath = c(get_custom_tbm_overlays(outdir)[[1]],get_custom_tbm_overlays(outdir)[[2]],get_custom_tbm_overlays(outdir)[[3]]),
                      #stat_overlay
                      outdir,
                      name = c(bs_stat_overlays$log_pvalues,bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues), alpha = 120)

                    private$render_atlas(cluster_iter, voxelcoord,
                                         atlaspath = bss_data@atlas_filename,
                                         outdir)
                    # }

                    private$render_table()
                    private$render_html(outdir,
                                        voxelcoord,
                                        overlay_name = c(bs_stat_overlays$log_pvalues,bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues))


                  }

                  private$save_rmd_preamble(file.path(outdir, "save_rmd.Rmd"),
                                            outdir,
                                            voxelcoord,
                                            overlay_name = "c(bs_stat_overlays$log_pvalues,bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues)")

                  # save_rmd<-file(file.path(outdir, "save_rmd.Rmd"))
                  # writeLines(private$render_table(), save_rmd)
                  # close(save_rmd)
                  #close(save_rmd)
                  # file.append("/Users/sjoshi/Desktop/tbm_anova/save_rmd.Rmd", "/Users/sjoshi/Desktop/tbm_anova/justhtml.Rmd")
                  # file.append(file.path(outdir, "save_rmd.Rmd"), file.path(outdir, "justhtml.Rmd"))
                  rmarkdown::render(paste0(outdir, "save_rmd.Rmd"))
                }
              ),

              private = list(

                render_overlay = function(voxelcoord_index,voxelcoord,atlaspath,overlaypath,outdir,name,alpha) {
                  for (stats_measure_index in 1:3) {
                    view_ax <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[stats_measure_index], " -o ",outdir,"/PNG_images/","ax",voxelcoord[[voxelcoord_index]][1], "_",name[stats_measure_index], ".png --slice ", voxelcoord[[voxelcoord_index]][1], " --", "ax", " -a ", alpha)
                    view_cor <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[stats_measure_index], " -o ",outdir,"/PNG_images/","cor",voxelcoord[[voxelcoord_index]][2], "_",name[stats_measure_index], ".png --slice ", voxelcoord[[voxelcoord_index]][2], " --", "cor", " -a ", alpha)
                    view_sag <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[stats_measure_index], " -o ",outdir,"/PNG_images/","sag",voxelcoord[[voxelcoord_index]][3], "_",name[stats_measure_index], ".png --slice ", voxelcoord[[voxelcoord_index]][3], " --", "sag", " -a ", alpha)

                    system(view_cor,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                    system(view_sag,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                    system(view_ax,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                  }
                  return(0)
                },
                render_atlas = function(voxelcoord_index,voxelcoord,atlaspath,outdir) {

                  view_at_ax <- paste0("/usr/local/bin/volblend -i ",atlaspath," --view 1 --slice ",voxelcoord[[voxelcoord_index]][1]," --flop -o ", outdir,"/PNG_images/","ax",voxelcoord[[voxelcoord_index]][1],"_atlas.png")

                  view_at_cor <- paste0("/usr/local/bin/volblend -i ",atlaspath," --view 2 --slice ",voxelcoord[[voxelcoord_index]][1]," --flop -o ", outdir,"/PNG_images/","cor",voxelcoord[[voxelcoord_index]][1],"_atlas.png")
                  view_at_sag <- paste0("/usr/local/bin/volblend -i ",atlaspath," --view 3 --slice ",voxelcoord[[voxelcoord_index]][1]," --flop -o ", outdir,"/PNG_images/","sag",voxelcoord[[voxelcoord_index]][1],"_atlas.png")

                  system(view_at_cor,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                  system(view_at_sag,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                  system(view_at_ax,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)

                  return(0)
                },
                render_table = function() {

                  t <- c("#test1","#test2","#test3","#test4","#test5")
                  Cluster <- 1:5
                  table <- data.frame(Cluster = 1:5, Vol_Size = c(33, 22.3, 21, 25, 30), voxelcoord = c(3,3,4,5,3), T_val =c(8,7.2,6,9.8,7.8))
                  table$Cluster <- paste0("[", table$Cluster, "](", t, ")")
                  knitr::kable(table[1:4], align=c(rep('l', 4)))
                },

                save_rmd_preamble = function(rmdfile, outdir, voxelcoord, overlay_name) {


                  # file_rmd_preamble <-file(rmdfile)
                  load_library <- "library(bssr)"
                  data_command_1 <- paste0("bss_data <- load_bss_data(type = '",bss_data@analysis_type,"', subjdir = '",bss_data@subjdir,"', csv = '", bss_data@csv,"', smooth = ",bss_data@smooth,")")
                  data_command_2 <- paste0("bss_model <- bss_anova(main_effect = '",bss_model@main_effect,"', covariates = '", bss_model@covariates,"', bss_data = bss_data)")
                  data_command_3 <- "bs_stat_overlays = list(log_pvalues = 'log_pvalues', log_pvalues_adjusted = 'log_pvalues_adjusted', tvalues = 'tvalues', pvalues = 'pvalues')"
                  user_input <- paste0(load_library,"\n",data_command_1,"\n",data_command_2,"\n",data_command_3,"\n")
                  templines <- deparse(private$render_html)
                  templines[1] <- "render_html = function(outdir, voxelcoord, overlay_name)"
                  sink(rmdfile, append=TRUE, type = "output")
                  cat("---\n")
                  cat("title: BSSR Report\n")
                  cat("output: html_document\n")
                  cat("---\n")
                  cat("```{r eval=TRUE, echo=FALSE, message=FALSE, results='hide', user_input_commands}\n")
                  writeLines(user_input)
                  cat("```\n")
                  cat("```{r echo=FALSE, warning=FALSE}\n")
                  writeLines(templines)
                  voxelcoord_char <- "list("
                  for (i in 1:length(voxelcoord)){
                    voxelcoord_char <- paste0(voxelcoord_char, "c(")
                    for (m in 1:length(voxelcoord[[i]])){
                      voxelcoord_char <- paste0(voxelcoord_char, voxelcoord[[i]][m],",")
                    }
                    voxelcoord_char <- paste0(substr(voxelcoord_char,1,nchar(voxelcoord_char)-1),"),")
                  }
                  voxelcoord_char <- paste0(substr(voxelcoord_char,1,nchar(voxelcoord_char)-1),")")
                  cat("render_html('",outdir,"', ", voxelcoord_char, ", ",overlay_name,")\n")
                  cat("```\n")
                  sink()

                  # writeLines(templines, file_rmd_preamble)

                  # cmd_text <- sprintf("render_html(outdir=\"%s/\", voxelcoord = list(c(90,90,90),c(107,107,107),c(120,120,120)),overlay_name = c(bs_stat_overlays$log_pvalues,bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues))",
                  # dirname(rmdfile))
                  # write(file_rmd_preamble, cmd_text)
                },

                ## another function will generate the names for the pngs

                render_html = function(outdir, voxelcoord, overlay_name) {
                  cbar_filename <- paste0(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(basename(bss_data@atlas_filename)), bs_stat_overlays$log_pvalues, sep = '_'), '_cbar.png')
                  cbar_filename1 <- paste0(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(basename(bss_data@atlas_filename)), bs_stat_overlays$log_pvalues_adjusted, sep = '_'), '_cbar.png')
                  cbar_filename2 <- paste0(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(basename(bss_data@atlas_filename)), bs_stat_overlays$tvalues, sep = '_'), '_cbar.png')

                  #function to make rmd work
                  get_render_image_filename <- function(individual_voxelcoord, overlay_name, brain_sector_index) {
                    if (brain_sector_index == 1) {
                      return(paste0("./PNG_images/cor", individual_voxelcoord,"_",overlay_name,".png"))
                    } else if (brain_sector_index == 2){
                      return(paste0("./PNG_images/ax", individual_voxelcoord,"_",overlay_name,".png"))
                    } else {
                      return(paste0("./PNG_images/sag", individual_voxelcoord,"_",overlay_name,".png"))
                    }

                  }

                  # Function to return a shiny image object
                  shiny_image = function(specific_voxelcoord, index_brain_sector, overlay_type, width){
                    align <- "left"
                    return(paste0("shiny::img(src='",get_render_image_filename(specific_voxelcoord,overlay_type, index_brain_sector),"', align = '", align, "', width = '", width,"')"))
                  }

                  # Function that creates a tabPanel
                  tab_panel = function(panel_type, voxelcoord_index){
                    width <- c("28.8%","24%","34.5%","11%")
                    overlay <- c(bs_stat_overlays$log_pvalues, bs_stat_overlays$log_pvalues_adjusted, bs_stat_overlays$tvalues)
                    cbar <- c(cbar_filename, cbar_filename1,cbar_filename2)
                    if (panel_type == "P-Values"){
                      overlay = overlay[1]
                      cbar = cbar[1]
                    } else if (panel_type == "Adjusted P-Values"){
                      overlay = overlay[2]
                      cbar = cbar[2]
                    } else if (panel_type == "T-Values"){
                      overlay = overlay[3]
                      cbar = cbar[3]
                    }
                    images <- ""
                    if (panel_type == "All") {
                      for (overlay_index in 1:length(overlay)) {
                        for (inner_coord_index in 1:3) {
                          images <- paste0(images, shiny_image(voxelcoord[[voxelcoord_index]][inner_coord_index], inner_coord_index, overlay[overlay_index], width[inner_coord_index]), ",")
                        }
                        images <- paste0(images, paste0("shiny::img(src=paste0('", cbar[overlay_index], "'), align='left', width = '", width[4],"'),"))
                      }
                    } else {
                      for (inner_coord_index in 1:3) {
                        images <- paste0(images, shiny_image(voxelcoord[[voxelcoord_index]][inner_coord_index], inner_coord_index,overlay, width[inner_coord_index]), ",")
                      }
                      images <- paste0(images, "shiny::img(src=paste0('", cbar, "'), align='left', width = '", width[4], "'),")
                    }
                    images <- substr(images,1,nchar(images)-1)
                    return(paste0("shiny::tabPanel(title_0 = '", panel_type, "', value = c('", panel_type, "'), shiny::p('", panel_type, "'), ",images,")"))
                  }

                  # Function that creates clusters for each panel
                  cluster_panels = function(voxelcoord){
                    panels <- paste0("shiny::tabPanel(title = shiny::h4(paste('Cluster 1'),shiny::tableOutput('data'),shiny::tabsetPanel(id = 'navbar',type = 'pills',",
                                                                                   tab_panel('All', 1),", ",
                                                                                   tab_panel('P-Values', 1), ", ",
                                                                                   tab_panel('Adjusted P-Values', 1), ", ",
                                                                                   tab_panel('T-Values', 1),")))")
                    for (voxelcoord_index in 2:length(voxelcoord)){
                      panels <- paste0(panels, ", shiny::tabPanel(title = shiny::h4(paste('Cluster ", voxelcoord_index, "'),shiny::tableOutput('data'),shiny::tabsetPanel(id = 'navbar',type = 'pills',",
                                                                                               tab_panel('All', voxelcoord_index),",",
                                                                                               tab_panel('P-Values', voxelcoord_index),",",
                                                                                               tab_panel('Adjusted P-Values', voxelcoord_index),",",
                                                                                               tab_panel('T-Values', voxelcoord_index),")))")

                    }
                    return(panels)
                  }

                  eval(parse(text= paste0("shiny::shinyUI(shiny::fluidPage(shinyjs::useShinyjs(),shiny::h3('Choose a cluster and an overlay below'),shiny::tabsetPanel(id = 'navbar',type = 'tabs',",
                        cluster_panels(voxelcoord),")))")))

                }

              )
  )

                  # shiny::shinyUI(
                  #   shiny::fluidPage(
                      #shinyjs::useShinyjs(),
                      #shiny::h3("Choose a cluster and an overlay below"),
                      # shiny::tabsetPanel(
                      #   id = "navbar",
                      #   type = "tabs",
                        # shiny::tabPanel(title = shiny::h4("Cluster 1"),
                        #                 shiny::tableOutput("data"),
                        #
                        #                 shiny::tabsetPanel(
                        #                   id = "navbar",
                        #
                        #                   type = "pills",
                        #                   shiny::tabPanel(title_0 = "All",
                        #
                        #                                   value = c("All"),
                        #                                   shiny::p("All"),
                        #                                   ##CLUSTER 1, cluster_iter = 1, inner_iter = 1:3
                        #                                   # P-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename), align="left", width = "11%"),
                        #
                        #                                   # Adj P-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename1), align="left", width = "11.2%"),
                        #                                   # T-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename2), align="left", width = "11.2%"),
                        #                                   # Atlas
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")),
                        #
                        #                   shiny::tabPanel(title_0 = "P-Values",
                        #                                   value = c("P-Values"),
                        #                                   shiny::p("P-Values"),
                        #                                   # P-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename), align="left", width = "11%")),
                        #                   shiny::tabPanel(title_0 = "Adjusted P-Values",
                        #                                   value = c("Adjusted P-Values"),
                        #                                   shiny::p("Adjusted P-Values"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename1), align="left", width = "11.5%")),
                        #                   shiny::tabPanel(title_0 = "T-Values",
                        #                                   value = c("T-Values"),
                        #                                   shiny::p("T-Values"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename2), align="left", width = "11%")),
                        #                   shiny::tabPanel(title_0 = "Atlas",
                        #                                   value = c("Atlas"),
                        #                                   shiny::p("Atlas"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(1,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")))),
                        #


                        # shiny::tabPanel(title = shiny::h4("Cluster 2"),
                        #                 shiny::tableOutput("data"),
                        #                 shiny::tabsetPanel(
                        #                   id = "navbar",
                        #                   type = "pills",
                        #                   shiny::tabPanel(title_0 = "All",
                        #                                   value = c("All"),
                        #                                   shiny::p("All"),
                        #                                   # P-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename), align="left", width = "11%"),
                        #                                   # Adj P-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename1), align="left", width = "11.2%"),
                        #                                   # T-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename2), align="left", width = "11%"),
                        #                                   # Atlas
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")),
                        #
                        #                   shiny::tabPanel(title_0 = "P-Values",
                        #                                   value = c("P-Values"),
                        #                                   shiny::p("P-Values"),
                        #                                   # P-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename), align="left", width = "11%")),
                        #                   shiny::tabPanel(title_0 = "Adjusted P-Values",
                        #                                   value = c("Adjusted P-Values"),
                        #                                   shiny::p("Adjusted P-Values"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename1), align="left", width = "11.5%")),
                        #                   shiny::tabPanel(title_0 = "T-Values",
                        #                                   value = c("T-Values"),
                        #                                   shiny::p("T-Values"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename2), align="left", width = "11%")),
                        #                   shiny::tabPanel(title_0 = "Atlas",
                        #                                   value = c("Atlas"),
                        #                                   shiny::p("Atlas"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(2,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")))),


                        # shiny::tabPanel(title = shiny::h4("Cluster 3"),
                        #                 shiny::tableOutput("data"),
                        #                 shiny::tabsetPanel(
                        #                   id = "navbar",
                        #                   type = "pills",
                        #                   shiny::tabPanel(title_0 = "All",
                        #
                        #                                   value = c("All"),
                        #                                   shiny::p("All"),
                        #                                   # P-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename), align="left", width = "11%"),
                        #                                   # Adj P-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename1), align="left", width = "11.2%"),
                        #                                   # T-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename2), align="left", width = "11%"),
                        #                                   # Atlas
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")),
                        #
                        #                   shiny::tabPanel(title_0 = "P-Values",
                        #                                   value = c("P-Values"),
                        #                                   shiny::p("P-Values"),
                        #                                   # P-values
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename), align="left", width = "11%")),
                        #                   shiny::tabPanel(title_0 = "Adjusted P-Values",
                        #                                   value = c("Adjusted P-Values"),
                        #                                   shiny::p("Adjusted P-Values"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename1), align="left", width = "11.5%")),
                        #                   shiny::tabPanel(title_0 = "T-Values",
                        #                                   value = c("T-Values"),
                        #                                   shiny::p("T-Values"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                        #                                   shiny::img(src=paste0(outdir, cbar_filename2), align="left", width = "11%")),
                        #                   shiny::tabPanel(title_0 = "Atlas",
                        #                                   value = c("Atlas"),
                        #                                   shiny::p("Atlas"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                        #                                   shiny::img(src=paste0(get_render_image_filename(3,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")))
                        #
                        # )
                      # )
                  #   )
                  # )

  #               }
  #             )
  # )


