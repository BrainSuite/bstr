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
              #inherit = BssRmdOutput,
              public = list(
                initialize = function(outdir ="/Users/sarapesavento/Desktop/tbm_anova") {
                  #super$
                  initialize(outdir)
                },
                save_out = function(bss_data, bss_model, voxelcoord, outdir) {

                  #create a folder to store png images in
                  dir.create("PNG_images")

                  for(i in 1:length(voxelcoord)) {
                    private$render_overlay(
                      ii=i,
                      voxelcoord,
                      atlaspath = bss_data@atlas_filename,
                      overlaypath = c("/Users/sarapesavento/Desktop/tbm_anova/bss_anova_age_mri.bfc.nii_log_pvalues.nii.gz","/Users/sarapesavento/Desktop/tbm_anova/bss_anova_age_mri.bfc.nii_log_pvalues_adjusted.nii.gz","/Users/sarapesavento/Desktop/tbm_anova/bss_anova_age_mri.bfc.nii_tvalues.nii.gz"),
                      outdir,
                      view = c("cor","sag","ax"), name = c("p","adjp","t"), alpha = 120)

                    private$render_atlas(ii = i, voxelcoord,
                                         filePath = bss_data@atlas_filename,
                                         outdir,
                                         view = c("cor","sag","ax"))
                    # }

                    private$render_table()    #these need to know previous names
                    private$render_html(outdir)
                  }
                }
              ),
              private = list(
                render_overlay = function(ii,voxelcoord,atlaspath,overlaypath,outdir,view,name,alpha) {

                  for (j in 1:3) {
                    #if (check error) { message, break}
                    #P VALUE
                    view_ax <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[j], " -o ",outdir,"/PNG_images/",view[1],voxelcoord[[ii]][1], "_",name[j], ".png --slice ", voxelcoord[[ii]][1], " --", view[1], " -a ", alpha)
                    view_cor <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[j], " -o ",outdir,"/PNG_images/",view[2],voxelcoord[[ii]][2], "_",name[j], ".png --slice ", voxelcoord[[ii]][2], " --", view[2], " -a ", alpha)
                    view_sag <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[j], " -o ",outdir,"/PNG_images/",view[3],voxelcoord[[ii]][3], "_",name[j], ".png --slice ", voxelcoord[[ii]][3], " --", view[3], " -a ", alpha)

                    system(view_cor,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                    system(view_sag,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                    system(view_ax,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                  }
                  return(0)
                  #return message for error. Test if returns code for an error (i.e. pstatmap0)
                },
                render_atlas = function(ii,voxelcoord,filePath,outdir,view) {

                  view_at_ax <- paste0("volblend -i ",filePath," --view 1 --slice ",voxelcoord[[ii]][1]," --flop -o ", outdir,"/PNG_images/",view[3],voxelcoord[[ii]][1],"_atlas.png")
                  view_at_cor <- paste0("volblend -i ",filePath," --view 2 --slice ",voxelcoord[[ii]][1]," --flop -o ", outdir,"/PNG_images/",view[1],voxelcoord[[ii]][1],"_atlas.png")
                  view_at_sag <- paste0("volblend -i ",filePath," --view 3 --slice ",voxelcoord[[ii]][1]," --flop -o ", outdir,"/PNG_images/",view[2],voxelcoord[[ii]][1],"_atlas.png")

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
                render_html = function(outdir) {

                  store_png_names <- list.files(path = paste0(outdir,"/PNG_images"))


                  library(png)
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

                                                          #automatically generate path like before, same function as before, can make seperate function

                                                          # P-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[19]), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[7]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[31]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Adj P-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[17]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[5]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[29]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.2%"),
                                                          # T-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[20]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[8]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[32]), align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Atlas
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[18]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[6]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[30]), align="left",width = "34.5%")),
                                          shiny::tabPanel(title_0 = "P-Values",
                                                          value = c("P-Values"),
                                                          shiny::p("P-Values"),
                                                          # P-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[19]), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[7]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[31]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Adjusted P-Values",
                                                          value = c("Adjusted P-Values"),
                                                          shiny::p("Adjusted P-Values"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[17]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[5]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[29]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%")),
                                          shiny::tabPanel(title_0 = "T-Values",
                                                          value = c("T-Values"),
                                                          shiny::p("T-Values"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[20]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[8]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[32]), align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Atlas",
                                                          value = c("Atlas"),
                                                          shiny::p("Atlas"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[18]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[6]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[26]), align="left",width = "34.5%")))),
                        shiny::tabPanel(title = shiny::h4("Cluster 2"),
                                        shiny::tableOutput("data"),

                                        shiny::tabsetPanel(
                                          id = "navbar",

                                          type = "pills",
                                          shiny::tabPanel(title_0 = "All",

                                                          value = c("All"),
                                                          shiny::p("All"),
                                                          # P-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[15]), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[3]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[27]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Adj P-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[13]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[1]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[25]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.2%"),
                                                          # T-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[16]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[4]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[28]), align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Atlas
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[14]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[2]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[30]), align="left",width = "34.5%")),
                                          shiny::tabPanel(title_0 = "P-Values",
                                                          value = c("P-Values"),
                                                          shiny::p("P-Values"),
                                                          # P-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[15]), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[3]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[27]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Adjusted P-Values",
                                                          value = c("Adjusted P-Values"),
                                                          shiny::p("Adjusted P-Values"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[13]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[1]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[25]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%")),
                                          shiny::tabPanel(title_0 = "T-Values",
                                                          value = c("T-Values"),
                                                          shiny::p("T-Values"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[16]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[4]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[28]), align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Atlas",
                                                          value = c("Atlas"),
                                                          shiny::p("Atlas"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[14]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[2]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[30]), align="left",width = "34.5%")))),


                        shiny::tabPanel(title = shiny::h4("Cluster 3"),
                                        shiny::tableOutput("data"),

                                        shiny::tabsetPanel(
                                          id = "navbar",

                                          type = "pills",
                                          shiny::tabPanel(title_0 = "All",

                                                          value = c("All"),
                                                          shiny::p("All"),
                                                          # P-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[23]), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[11]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[35]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Adj P-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[21]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[9]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[33]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.2%"),
                                                          # T-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[24]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[12]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[36]), align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Atlas
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[22]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[10]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[34]), align="left",width = "34.5%")),
                                          shiny::tabPanel(title_0 = "P-Values",
                                                          value = c("P-Values"),
                                                          shiny::p("P-Values"),
                                                          # P-values
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[23]), align="left", width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[11]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[35]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Adjusted P-Values",
                                                          value = c("Adjusted P-Values"),
                                                          shiny::p("Adjusted P-Values"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[21]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[9]), align="left", width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[33]), align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%")),
                                          shiny::tabPanel(title_0 = "T-Values",
                                                          value = c("T-Values"),
                                                          shiny::p("T-Values"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[24]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[12]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[36]), align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Atlas",
                                                          value = c("Atlas"),
                                                          shiny::p("Atlas"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[22]), align="left",width = "28.8%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[10]), align="left",width = "24%"),
                                                          shiny::img(src=paste0('./PNG_images/',store_png_names[34]), align="left",width = "34.5%")))

                        )
                      )
                    )
                  )

                }

              )
  )

