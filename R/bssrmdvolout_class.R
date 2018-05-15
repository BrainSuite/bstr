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
                save_out = function(bss_data, bss_model, voxelcoord, outdir) {


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
                      view = c("ax", "cor", "sag"), name = c(bs_stat_overlays$log_pvalues,bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues), alpha = 120)

                    private$render_atlas(cluster_iter, voxelcoord,
                                         atlaspath = bss_data@atlas_filename,
                                         outdir,
                                         view = c("ax", "cor", "sag"))
                    # }

                    private$render_table()    #these need to know previous names
                    rmdout1 <- private$render_table()
                    private$render_html(outdir,
                                        voxelcoord,
                                        view = c("ax", "cor", "sag"),
                                        overlay_name = c(bs_stat_overlays$log_pvalues,bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues))

                  rmdout2 <-private$render_html(outdir,
                                                  voxelcoord,
                                                  view = c("ax", "cor", "sag"),
                                                  overlay_name = c(bs_stat_overlays$log_pvalues,bs_stat_overlays$log_pvalues_adjusted,bs_stat_overlays$tvalues))

                  }

                  save_rmd<-file(file.path(outdir, "save_rmd.Rmd"))
                  writeLines(text = "this text will show up", save_rmd)
                  close(save_rmd)
                  rmarkdown::render("/Users/sarapesavento/Desktop/tbm_anova/save_rmd.Rmd")

                }


              ),

              private = list(

                render_overlay = function(cluster_iter,voxelcoord,atlaspath,overlaypath,outdir,view,name,alpha) {
                  for (inner_iter in 1:3) {
                    view_ax <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[inner_iter], " -o ",outdir,"/PNG_images/",view[1],voxelcoord[[cluster_iter]][1], "_",name[inner_iter], ".png --slice ", voxelcoord[[cluster_iter]][1], " --", view[1], " -a ", alpha)
                    view_cor <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[inner_iter], " -o ",outdir,"/PNG_images/",view[2],voxelcoord[[cluster_iter]][2], "_",name[inner_iter], ".png --slice ", voxelcoord[[cluster_iter]][2], " --", view[2], " -a ", alpha)
                    view_sag <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[inner_iter], " -o ",outdir,"/PNG_images/",view[3],voxelcoord[[cluster_iter]][3], "_",name[inner_iter], ".png --slice ", voxelcoord[[cluster_iter]][3], " --", view[3], " -a ", alpha)

                    system(view_cor,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                    system(view_sag,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                    system(view_ax,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                  }
                  return(0)
                  #return message for error. Test if returns code for an error (i.e. pstatmap0)
                },
                render_atlas = function(cluster_iter,voxelcoord,atlaspath,outdir,view) {

                  view_at_ax <- paste0("/usr/local/bin/volblend -i ",atlaspath," --view 1 --slice ",voxelcoord[[cluster_iter]][1]," --flop -o ", outdir,"/PNG_images/",view[1],voxelcoord[[cluster_iter]][1],"_atlas.png")

                  view_at_cor <- paste0("/usr/local/bin/volblend -i ",atlaspath," --view 2 --slice ",voxelcoord[[cluster_iter]][1]," --flop -o ", outdir,"/PNG_images/",view[2],voxelcoord[[cluster_iter]][1],"_atlas.png")
                  view_at_sag <- paste0("/usr/local/bin/volblend -i ",atlaspath," --view 3 --slice ",voxelcoord[[cluster_iter]][1]," --flop -o ", outdir,"/PNG_images/",view[3],voxelcoord[[cluster_iter]][1],"_atlas.png")

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
                ## another function will generate the names for the pngs

                render_html = function(outdir, voxelcoord, view,overlay_name) {

                  image <- image_read("/Users/sarapesavento/Desktop/tbm_anova/bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf")
                  ppng <- image_convert(image, format = "PNG", type = NULL, colorspace = NULL,
                                         depth = NULL, antialias = NULL)
                  image_write(ppng, path ="/Users/sarapesavento/Desktop/tbm_anova/pvalues.png", format = "png")

                  image <- image_read("/Users/sarapesavento/Desktop/tbm_anova/bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf")
                  adppng <- image_convert(image, format = "PNG", type = NULL, colorspace = NULL,
                                         depth = NULL, antialias = NULL)
                  image_write(adppng, path ="/Users/sarapesavento/Desktop/tbm_anova/adjpvalues.png", format = "png")

                  image <- image_read("/Users/sarapesavento/Desktop/tbm_anova/bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf")
                  tpng <- image_convert(image, format = "PNG", type = NULL, colorspace = NULL,
                                         depth = NULL, antialias = NULL)
                  image_write(tpng, path ="/Users/sarapesavento/Desktop/tbm_anova/adjpvalues.png", format = "png")


                  #function to make rmd work
                    get_render_image_filename <- function(view, cluster_iter, inner, voxelcoord, overlay_name) {
                      return(paste0("./PNG_images/",view, voxelcoord[[cluster_iter]][inner],"_",overlay_name,".png"))
                    }


                  shiny::shinyUI(
                    shiny::fluidPage(
                      shinyjs::useShinyjs(),
                      shiny::h3("Choose a cluster and an overlay below"),
                      shiny::tabsetPanel(
                        id = "navbar",
                        type = "tabs",
                        shiny::tabPanel(title = shiny::h4("Cluster 1"),
                                        shiny::tableOutput("data"),

                                        shiny::tabsetPanel(
                                          id = "navbar",

                                          type = "pills",
                                          shiny::tabPanel(title_0 = "All",

                                                          value = c("All"),
                                                          shiny::p("All"),
                                                          bs_stat_overlays <- list(
                                                            log_pvalues = "log_pvalues",
                                                            log_pvalues_adjusted = "log_pvalues_adjusted",
                                                            tvalues = "tvalues",
                                                            pvalues = "pvalues"
                                                          ),

                                                          ##CLUSTER 1, cluster_iter = 1, inner_iter = 1:3
                                                          view = c("ax", "cor", "sag"),
                                                          # P-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],1,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],1,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],1,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%"),

                                                          # Adj P-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],1,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],1,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],1,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.2%"),
                                                          # T-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],1,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],1,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],1,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Atlas
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],1,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],1,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],1,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")),

                                          shiny::tabPanel(title_0 = "P-Values",
                                                          value = c("P-Values"),
                                                          shiny::p("P-Values"),
                                                          # P-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],1,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],1,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],1,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Adjusted P-Values",
                                                          value = c("Adjusted P-Values"),
                                                          shiny::p("Adjusted P-Values"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],1,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],1,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],1,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%")),
                                          shiny::tabPanel(title_0 = "T-Values",
                                                          value = c("T-Values"),
                                                          shiny::p("T-Values"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],1,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],1,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],1,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Atlas",
                                                          value = c("Atlas"),
                                                          shiny::p("Atlas"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],1,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],1,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],1,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")))),



                        shiny::tabPanel(title = shiny::h4("Cluster 2"),
                                        shiny::tableOutput("data"),
                                        shiny::tabsetPanel(
                                          id = "navbar",
                                          type = "pills",
                                          shiny::tabPanel(title_0 = "All",
                                                          value = c("All"),
                                                          shiny::p("All"),
                                                          # P-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],2,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],2,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],2,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Adj P-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],2,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],2,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],2,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.2%"),
                                                          # T-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],2,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],2,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],2,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Atlas
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],2,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],2,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],2,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")),

                                          shiny::tabPanel(title_0 = "P-Values",
                                                          value = c("P-Values"),
                                                          shiny::p("P-Values"),
                                                          # P-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],2,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],2,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],2,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Adjusted P-Values",
                                                          value = c("Adjusted P-Values"),
                                                          shiny::p("Adjusted P-Values"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],2,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],2,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],2,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%")),
                                          shiny::tabPanel(title_0 = "T-Values",
                                                          value = c("T-Values"),
                                                          shiny::p("T-Values"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],2,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],2,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],2,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Atlas",
                                                          value = c("Atlas"),
                                                          shiny::p("Atlas"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],2,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],2,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],2,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")))),


                        shiny::tabPanel(title = shiny::h4("Cluster 3"),
                                        shiny::tableOutput("data"),
                                        shiny::tabsetPanel(
                                          id = "navbar",
                                          type = "pills",
                                          shiny::tabPanel(title_0 = "All",

                                                          value = c("All"),
                                                          shiny::p("All"),
                                                          # P-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],3,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],3,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],3,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Adj P-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],3,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],3,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],3,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.2%"),
                                                          # T-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],3,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],3,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],3,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Atlas
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],3,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],3,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],3,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")),

                                          shiny::tabPanel(title_0 = "P-Values",
                                                          value = c("P-Values"),
                                                          shiny::p("P-Values"),
                                                          # P-values
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],3,inner=1,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],3,inner=2,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],3,inner=3,voxelcoord,bs_stat_overlays$log_pvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Adjusted P-Values",
                                                          value = c("Adjusted P-Values"),
                                                          shiny::p("Adjusted P-Values"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],3,inner=1,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],3,inner=2,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],3,inner=3,voxelcoord,bs_stat_overlays$log_pvalues_adjusted)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%")),
                                          shiny::tabPanel(title_0 = "T-Values",
                                                          value = c("T-Values"),
                                                          shiny::p("T-Values"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],3,inner=1,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],3,inner=2,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],3,inner=3,voxelcoord,bs_stat_overlays$tvalues)), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Atlas",
                                                          value = c("Atlas"),
                                                          shiny::p("Atlas"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[2],3,inner=1,voxelcoord,"atlas")), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[1],3,inner=2,voxelcoord,"atlas")), align="left", width = "24%"),
                                                          shiny::img(src=paste0(get_render_image_filename(view[3],3,inner=3,voxelcoord,"atlas")), align="left", width = "34.5%")))

                        )
                      )
                    )
                  )

                }
              )
  )

