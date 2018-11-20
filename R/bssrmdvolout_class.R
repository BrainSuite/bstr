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

                    p_overlay <- paste0(outdir, bss_model@model_type, "_", bss_model@main_effect,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$log_pvalues, bss_data@data_type)
                    adjp_overlay <- paste0(outdir, bss_model@model_type, "_", bss_model@main_effect,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$log_pvalues_adjusted, bss_data@data_type)
                    t_overlay <- paste0(outdir, bss_model@model_type, "_", bss_model@main_effect,"_",tools::file_path_sans_ext(
                      basename(bss_data@atlas_filename)),"_", bs_stat_overlays$tvalues, bss_data@data_type)

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
                  view_order <- c("sag","cor","ax")
                  for (stats_measure_index in 1:3) {
                    for (view in 1:3){
                       current_view <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[stats_measure_index], " -o ",outdir,"/PNG_images/",view_order[view],voxelcoord[[voxelcoord_index]][view], "_",name[stats_measure_index], ".png --slice ", voxelcoord[[voxelcoord_index]][view], " --", view_order[view], " -a ", alpha)
                       system(current_view,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                    }
                  }
                  return(0)
                },
                render_atlas = function(voxelcoord_index,voxelcoord,atlaspath,outdir) {
                  view_name <- c("ax","cor","sag")
                  view_order <- c(3,2,1)
                  for (view_iter in 1:3){
                    current_view <- paste0("/usr/local/bin/volblend -i ",atlaspath," --view ", view_iter," --slice ",voxelcoord[[voxelcoord_index]][view_order[view_iter]]," --flop -o ", outdir,"/PNG_images/",view_name[view_iter], voxelcoord[[voxelcoord_index]][view_order[view_iter]],"_atlas.png")
                    system(current_view,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                  }
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
                  data_command_1 <- paste0("bss_data <- load_bss_data(type = '",bss_data@analysis_type,"', subjdir = '",bss_data@subjdir,"', csv = '", bss_data@csv,"', measure = '", bss_data@measure,"', smooth = ",bss_data@smooth,")")
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
                  cat("<style>\n\n")
                  cat("table, td, th {")
                  cat("border: none;")
                  cat("padding-left: 1em;")
                  cat("padding-right: 1em;")
                  cat("min-width: 50%;")
                  cat("margin-left: auto;")
                  cat("margin-right: auto;")
                  cat("margin-top: 1em;")
                  cat("margin-bottom: 1em;")
                  cat("}\n\n")
                  cat("</style>\n")
                  cat("```{r eval=TRUE, echo=FALSE, message=FALSE, results='hide', user_input_commands}\n")
                  writeLines(user_input)
                  cat("```\n")
                  cat("```{r echo=FALSE, warning=FALSE}\n")
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
                  cat("render_html('",outdir,"', ", voxelcoord_char, ", ",overlay_name,")\n")
                  cat("```\n")
                  sink()

                },

                ## another function will generate the names for the pngs

                render_html = function(outdir, voxelcoord, overlay_name) {
                  cbar <- vector("list",3)
                  overlay <- c(bs_stat_overlays$log_pvalues, bs_stat_overlays$log_pvalues_adjusted, bs_stat_overlays$tvalues)
                  for (cbar_index in 1:3){
                    cbar[[cbar_index]] <- paste0(paste(bss_model@model_type, bss_model@main_effect, tools::file_path_sans_ext(basename(bss_data@atlas_filename)), overlay[cbar_index], sep = '_'), '_cbar.png')

                  }
                  vox_table <- read.table("~/cluster.tsv",header=F,sep="\t")
                  vox_table[,1] <- 1:nrow(vox_table)
                  colnames(vox_table) <- c("Cluster Number","Number of Voxels","Maximum Value","X Coord","Y Coord","Z Coord")

                  #function to make rmd work
                  get_render_image_filename <- function(individual_voxelcoord, overlay_name, brain_sector_index) {
                    view_order <- c("sag","cor","ax")
                    return(paste0("./PNG_images/", view_order[brain_sector_index], individual_voxelcoord,"_",overlay_name,".png"))
                  }

                  # Function to return a shiny image object
                  shiny_image = function(specific_voxelcoord, index_brain_sector, overlay_type, width){
                    align <- "left"
                    return(paste0("shiny::img(src='",get_render_image_filename(specific_voxelcoord,overlay_type, index_brain_sector),"', align = '", align, "', width = '", width,"')"))
                  }

                  # Function that creates a tabPanel
                  tab_panel = function(panel_type, voxelcoord_index){
                    width <- c("34.5%","28.8%","24%","11%")
                    overlay <- c(bs_stat_overlays$log_pvalues, bs_stat_overlays$log_pvalues_adjusted, bs_stat_overlays$tvalues)
                    panel_names <- c("P-Values","Adjusted P-Values","T-Values")
                    images <- ""
                    for (inner_coord_index in 1:3) {
                      images <- paste0(images, shiny_image(voxelcoord[[voxelcoord_index]][inner_coord_index], inner_coord_index,overlay[panel_type], width[inner_coord_index]), ",")
                    }
                    images <- paste0(images, "shiny::img(src=paste0('", cbar[panel_type], "'), align='left', width = '", width[4], "'),")
                    images <- substr(images,1,nchar(images)-1)
                    return(paste0("shiny::tabPanel(title_0 = '", panel_names[panel_type], "', value = c('", panel_names[panel_type], "'), shiny::p('", panel_names[panel_type], "'), ",images,")"))
                  }

                  # Function that creates clusters for each panel
                  cluster_panels = function(voxelcoord){
                    panels <- ""
                    for (voxelcoord_index in 1:length(voxelcoord)){
                      panels <- paste0(panels, ", shiny::tabPanel(title = shiny::h4(paste('Cluster ", voxelcoord_index, ": Voxel Coordinate (",voxelcoord[[voxelcoord_index]][1],",",voxelcoord[[voxelcoord_index]][2],",",voxelcoord[[voxelcoord_index]][3],")'),shiny::tableOutput('data'),shiny::tabsetPanel(id = 'navbar',type = 'pills',")
                      for (tab_panel_type in 1:3){
                        panels <- paste0(panels,tab_panel(tab_panel_type, voxelcoord_index),",")
                      }
                      panels <- paste0(substr(panels,1,nchar(panels)-1),")))")
                    }
                    return(substr(panels,1,nchar(panels)))
                  }

                  # Function that creates the table tab
                  table_tab <- function(){
                    table <- "shiny::tags$tr(shiny::tags$th('"
                    for (table_header in 1:length(colnames(vox_table))){
                      table <- paste0(table,colnames(vox_table)[table_header],"    '),shiny::tags$th('")
                    }
                    table <- paste0(substr(table,1,nchar(table)-17),"),")
                    for (colwise_data in 1:nrow(vox_table)){
                      table <- paste0(table, "shiny::tags$tr(")
                      for (rowwise_data in 1:length(colnames(vox_table))){
                        table <- paste0(table,"shiny::tags$td(vox_table[",colwise_data,",",rowwise_data,"]),")
                      }
                      table <- paste0(substr(table,1,nchar(table)-1),"),")
                    }
                    table <- paste0(substr(table,1,nchar(table)-1),")")
                    panel <- paste0("shiny::tabPanel(title = shiny::h4(paste(''),shiny::tableOutput('data'),shiny::tabsetPanel(id = 'navbar',type = 'pills',")
                    panel <- paste0(panel,"shiny::tabPanel(title_0 = 'Voxel Coordinate Table' , value = c('Voxel Coordinate Table'), shiny::p('Voxel Coordinate Table'), shiny::tags$table(",table,"))))")
                    return(panel)
                  }

                  eval(parse(text= paste0("shiny::shinyUI(shiny::fluidPage(shinyjs::useShinyjs(),shiny::h3('Choose a cluster and an overlay below'),shiny::tabsetPanel(id = 'navbar',type = 'tabs',",
                        table_tab(), cluster_panels(voxelcoord),")))")))

                }

              )
  )

