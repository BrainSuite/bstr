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

                  get_custom_tbm_overlays <- function(outdir) {
                    p_overlay <- paste(outdir, "/bss_anova_age_mri.bfc.nii_",bs_stat_overlays$log_pvalues,bs_data_types$nifti_image, sep="") #log_pvalues.nii.gz
                    adjp_overlay <- paste(outdir, "/bss_anova_age_mri.bfc.nii_",bs_stat_overlays$log_pvalues_adjusted,bs_data_types$nifti_image, sep="") #log_pvalues.nii.gz
                    t_overlay <- paste(outdir, "/bss_anova_age_mri.bfc.nii_",bs_stat_overlays$tvalues,bs_data_types$nifti_image, sep="") #log_pvalues.nii.gz

                    return(list("p_overlay" = p_overlay, "adjp_overlay" = adjp_overlay, "t_overlay" = t_overlay))
                  }

                  #create a folder to store png images in
                  dir.create(paste0(outdir,"/PNG_images"))

                  for(out_iter in 1:length(voxelcoord)) {
                    private$render_png_names(outdir,
                                             view = c("cor","sag","ax"),
                                             voxelcoord,
                                             out_iter,
                                             name = c("p","adjp","t"))
                    private$render_overlay(
                      out_iter,
                      voxelcoord,
                      atlaspath = bss_data@atlas_filename,
                      overlaypath = c(get_custom_tbm_overlays(outdir)[[1]],get_custom_tbm_overlays(outdir)[[2]],get_custom_tbm_overlays(outdir)[[3]]),
                      #stat_overlay
                      outdir,
                      view = c("cor","sag","ax"), name = c("p","adjp","t"), alpha = 120)

                    private$render_atlas(out_iter, voxelcoord,
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

                  render_png_names = function(outdir, view, voxelcoord, out_iter, name) {
                    r = 1:9
                    png_names <- c()
                    for (i in 1:3) {
                      png_names[r[1]]<- paste0(outdir,"/PNG_images/",view[1],voxelcoord[[out_iter]][1], "_",name[i], ".png")
                      png_names[r[2]] <- paste0(outdir,"/PNG_images/",view[2],voxelcoord[[out_iter]][2], "_",name[i], ".png")
                      png_names[r[3]] <- paste0(outdir,"/PNG_images/",view[3],voxelcoord[[out_iter]][3], "_",name[i], ".png")
                      r + 3
                    }
                    png_names[10] <- paste0(outdir,"/PNG_images/",view[3],voxelcoord[[out_iter]][1],"_atlas.png")
                    png_names[11] <- paste0(outdir,"/PNG_images/",view[1],voxelcoord[[out_iter]][1],"_atlas.png")
                    png_names[12] <- paste0(outdir,"/PNG_images/",view[2],voxelcoord[[out_iter]][1],"_atlas.png")

                    paste(0)
                  },
                render_overlay = function(out_iter,voxelcoord,atlaspath,overlaypath,outdir,view,name,alpha) {
                  for (inner_iter in 1:3) {
                    #if (check error) { message, break}
                    #P VALUE
                    view_ax <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[inner_iter], " -o ",outdir,"/PNG_images/",view[2],voxelcoord[[out_iter]][1], "_",name[inner_iter], ".png --slice ", voxelcoord[[out_iter]][1], " --", view[1], " -a ", alpha)
                    view_cor <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[inner_iter], " -o ",outdir,"/PNG_images/",view[2],voxelcoord[[out_iter]][2], "_",name[inner_iter], ".png --slice ", voxelcoord[[out_iter]][2], " --", view[2], " -a ", alpha)
                    view_sag <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[inner_iter], " -o ",outdir,"/PNG_images/",view[3],voxelcoord[[out_iter]][3], "_",name[inner_iter], ".png --slice ", voxelcoord[[out_iter]][3], " --", view[3], " -a ", alpha)


                    system(view_cor,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                    system(view_sag,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                    system(view_ax,intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE, input=NULL)
                  }
                  return(0)
                  #return message for error. Test if returns code for an error (i.e. pstatmap0)
                },
                render_atlas = function(out_iter,voxelcoord,filePath,outdir,view) {

                  view_at_ax <- paste0("/usr/local/bin/volblend -i ",filePath," --view 1 --slice ",voxelcoord[[out_iter]][1]," --flop -o ", outdir,"/PNG_images/",view[3],voxelcoord[[out_iter]][1],"_atlas.png")
                  view_at_cor <- paste0("/usr/local/bin/volblend -i ",filePath," --view 2 --slice ",voxelcoord[[out_iter]][1]," --flop -o ", outdir,"/PNG_images/",view[1],voxelcoord[[out_iter]][1],"_atlas.png")
                  view_at_sag <- paste0("/usr/local/bin/volblend -i ",filePath," --view 3 --slice ",voxelcoord[[out_iter]][1]," --flop -o ", outdir,"/PNG_images/",view[2],voxelcoord[[out_iter]][1],"_atlas.png")

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

                 render_html = function(outdir, voxelcoord) {
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
                                                          #shiny::img(src=paste0('./PNG_images/cor',voxelcoord[[1]][1],'p.png'), align="left", width = "28.8%"),
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

bssrmd_volout <- BssRmdVolumeOutput$new(outdir="/Users/sarapesavento/Desktop/tbm_anova")
bssrmd_volout

library(bssr)
bss_data <- load_bss_data(type="tbm",
                          subjdir = "/Users/sarapesavento/Desktop/RocklandSample25OHBM",
                          csv = "/Users/sarapesavento/Desktop/RocklandSample25OHBM/OHBM_workshop_demographics25.tsv",smooth=2)
bss_model <- bss_anova(main_effect = "age",
                       covariates = "sex",
                       bss_data = bss_data)
bssrmd_volout$save_out(bss_data, bss_model, voxelcoord = list(c(90,90,90),c(107,107,107),c(120,120,120)), outdir="/Users/sarapesavento/Desktop/tbm_anova")

