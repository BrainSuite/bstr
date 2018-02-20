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
                initialize = function(outdir = NULL) {
                  #super$
                  initialize(outdir)
                  },
                save_out = function(bss_data, bss_model) {
                  out_list = list()
                  voxel_list <- list(c(90,90,90), c(107,107,107), c(120,120,120))
                  #ex list_voxel_centers < get_clusters_from_stats(bss_data, bss_model)
                  for(i in 1:length(voxel_list)) {
                       out_list[[i]] <- render_overlay(voxelcoord = voxel_list[[i]],
                                                       firstpath = "/Users/sarapesavento/Desktop/tbm_anova/mri.bfc.nii.gz",
                                                       overlaypath = c("/Users/sarapesavento/Desktop/tbm_anova/bss_anova_age_mri.bfc.nii_log_pvalues.nii.gz","/Users/sarapesavento/Desktop/tbm_anova/bss_anova_age_mri.bfc.nii_log_pvalues_adjusted.nii.gz","/Users/sarapesavento/Desktop/tbm_anova/bss_anova_age_mri.bfc.nii_tvalues.nii.gz"),
                                                       view = c("cor","sag","ax"),
                                                       name = c("p","adjp","t"),
                                                       alpha = 120)
                       render_atlas(voxelcoord = voxel_list[[i]],
                                    filePath = "/Applications/BrainSuite17a/svreg/BCI-DNI_brain_atlas/BCI-DNI_brain.bfc.nii.gz",
                                    folderPath = "/Users/sarapesavento/Desktop/tbm_anova/",
                                    atlasName="atlas",
                                    view = c("cor","sag","ax"))
                    }

                  render_table()    #these need to know previous names
                  render_html()

                    render_atlas(voxelcoord = c(107,107,107), filePath = "/Applications/BrainSuite17a/svreg/BrainSuiteAtlas1/mri.bfc.nii.gz",
                                 folderPath = "/Users/sarapesavento/Desktop/tbm_anova/", atlasName="atlas")
                }
              ),
              private = list(
                render_overlay = function(voxelcoord, atlaspath,overlaypath,view,name,alpha) {
                  intern = FALSE
                  ignore.stdout = FALSE
                  ignore.stderr = FALSE
                  wait = TRUE
                  input = NULL

                  #P VALUE
                  view1 <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[1], " -o ",view[1],voxelcoord[1], "_",name[1], ".png --slice ", voxelcoord[1], " --", view[1], " -a ", alpha)
                  view2 <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[1], " -o ",view[2],voxelcoord[1], "_",name[1], ".png --slice ", voxelcoord[1], " --", view[2], " -a ", alpha)
                  view3 <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[1], " -o ",view[3],voxelcoord[1], "_",name[1], ".png --slice ", voxelcoord[1], " --", view[3], " -a ", alpha)

                  #ADJUSTED P
                  view4 <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[2], " -o ",view[1],voxelcoord[1], "_",name[2], ".png --slice ", voxelcoord[1], " --", view[1], " -a ", alpha)
                  view5 <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[2], " -o ",view[2],voxelcoord[1], "_",name[2], ".png --slice ", voxelcoord[1], " --", view[2], " -a ", alpha)
                  view6 <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[2], " -o ",view[3],voxelcoord[1], "_",name[2], ".png --slice ", voxelcoord[1], " --", view[3], " -a ", alpha)

                  # T VALUE
                  view7 <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[3], " -o ",view[1],voxelcoord[1], "_",name[3], ".png --slice ", voxelcoord[1], " --", view[1], " -a ", alpha)
                  view8 <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[3], " -o ",view[2],voxelcoord[1], "_",name[3], ".png --slice ", voxelcoord[1], " --", view[2], " -a ", alpha)
                  view9 <- paste0("pstatmap --atlas ", atlaspath, " -i ",overlaypath[3], " -o ",view[3],voxelcoord[1], "_",name[3], ".png --slice ", voxelcoord[1], " --", view[3], " -a ", alpha)
#store all png's into their own seperate folder
#take system calls out of return
#return()
                  return(c(system(view1, intern, ignore.stdout, ignore.stderr, wait, input),
                           system(view2, intern, ignore.stdout, ignore.stderr, wait, input),
                           system(view3, intern, ignore.stdout, ignore.stderr, wait, input)))
#return message for error. Test if returns code for an error (i.e. pstatmap0)
                },
                render_atlas = function(voxelcoord,filePath,folderPath,atlasName,view) {
                  intern = FALSE
                  ignore.stdout = FALSE
                  ignore.stderr = FALSE
                  wait = TRUE
                  input = NULL

                  view1 <- sprintf("volblend -i ",filePath," --view 1 --slice ",voxelcoord[3]," --flop -o ", folderPath, atlasName, "atlas",view[1],voxelcoord[1],".png")
                  view2 <- sprintf("volblend -i ",filePath," --view 2 --slice ",voxelcoord[2]," --flop -o ", folderPath, atlasName, "atlas",view[2],voxelcoord[2],".png")
                  view3 <- sprintf("volblend -i ",filePath," --view 3 --slice ",voxelcoord[1]," --flop -o ", folderPath, atlasName, "atlas",view[3],voxelcoord[3],".png")
                  return(c(system(view1, intern, ignore.stdout, ignore.stderr, wait, input), system(view2, intern, ignore.stdout, ignore.stderr, wait, input), system(view3, intern, ignore.stdout, ignore.stderr, wait, input)))
                },
                render_table = function() {

                  t <- c("#test1","#test2","#test3","#test4","#test5")
                  Cluster <- 1:5
                  table <- data.frame(Cluster = 1:5, Vol_Size = c(33, 22.3, 21, 25, 30), voxelcoord = c(3,3,4,5,3), T_val =c(8,7.2,6,9.8,7.8))
                  table$Cluster <- paste0("[", table$Cluster, "](", t, ")")
                  knitr::kable(table[1:4], align=c(rep('l', 4)))
                },
                render_html = function() {
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

                                                          # P-values
                                                          shiny::img(src='./cor120_p.png', align="left", width = "28.8%"),
                                                          shiny::img(src='./ax120_p.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag120_p.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Adj P-values
                                                          shiny::img(src='./cor120_adjp.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax120_adjp.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag120_adjp.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.2%"),
                                                          # T-values
                                                          shiny::img(src='./cor120_t.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax120_t.png', align="left",width = "24%"),
                                                          shiny::img(src='./sag120_t.png', align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Atlas
                                                          shiny::img(src='./atlas.107.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./atlas.106.png', align="left",width = "24%"),
                                                          shiny::img(src='./atlas.108.png', align="left",width = "34.5%")),
                                          shiny::tabPanel(title_0 = "P-Values",
                                                          value = c("P-Values"),
                                                          shiny::p("P-Values"),
                                                          # P-values
                                                          shiny::img(src='./cor120_p.png', align="left", width = "28.8%"),
                                                          shiny::img(src='./ax120_p.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag120_p.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Adjusted P-Values",
                                                          value = c("Adjusted P-Values"),
                                                          shiny::p("Adjusted P-Values"),
                                                          shiny::img(src='./cor120_adjp.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax120_adjp.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag120_adjp.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%")),
                                          shiny::tabPanel(title_0 = "T-Values",
                                                          value = c("T-Values"),
                                                          shiny::p("T-Values"),
                                                          shiny::img(src='./cor120_t.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax120_t.png', align="left",width = "24%"),
                                                          shiny::img(src='./sag120_t.png', align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Atlas",
                                                          value = c("Atlas"),
                                                          shiny::p("Atlas"),
                                                          shiny::img(src='./atlas.107.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./atlas.106.png', align="left",width = "24%"),
                                                          shiny::img(src='./atlas.108.png', align="left",width = "34.5%")))),
                        shiny::tabPanel(title = shiny::h4("Cluster 2"),
                                        shiny::tableOutput("data"),

                                        shiny::tabsetPanel(
                                          id = "navbar",

                                          type = "pills",
                                          shiny::tabPanel(title_0 = "All",

                                                          value = c("All"),
                                                          shiny::p("All"),
                                                          # P-values
                                                          shiny::img(src='./cor100_p.png', align="left", width = "28.8%"),
                                                          shiny::img(src='./ax100_p.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag100_p.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Adj P-values
                                                          shiny::img(src='./cor100_adjp.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax100_adjp.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag100_adjp.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.2%"),
                                                          # T-values
                                                          shiny::img(src='./cor100_t.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax100_t.png', align="left",width = "24%"),
                                                          shiny::img(src='./sag100_t.png', align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Atlas
                                                          shiny::img(src='./atlas.107.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./atlas.106.png', align="left",width = "24%"),
                                                          shiny::img(src='./atlas.108.png', align="left",width = "34.5%")),
                                          shiny::tabPanel(title_0 = "P-Values",
                                                          value = c("P-Values"),
                                                          shiny::p("P-Values"),
                                                          # P-values
                                                          shiny::img(src='./cor100_p.png', align="left", width = "28.8%"),
                                                          shiny::img(src='./ax100_p.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag100_p.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Adjusted P-Values",
                                                          value = c("Adjusted P-Values"),
                                                          shiny::p("Adjusted P-Values"),
                                                          shiny::img(src='./cor100_adjp.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax100_adjp.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag100_adjp.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%")),
                                          shiny::tabPanel(title_0 = "T-Values",
                                                          value = c("T-Values"),
                                                          shiny::p("T-Values"),
                                                          shiny::img(src='./cor100_t.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./sag100_t.png', align="left",width = "28.8%"),
                                                          hiny::img(src='./ax100_t.png', align="left",width = "24%"),
                                                          shiny::img(src='./sag100_t.png', align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Atlas",
                                                          value = c("Atlas"),
                                                          shiny::p("Atlas"),
                                                          shiny::img(src='./atlas.107.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./atlas.106.png', align="left",width = "24%"),
                                                          shiny::img(src='./atlas.108.png', align="left",width = "34.5%")))),


                        shiny::tabPanel(title = shiny::h4("Cluster 3"),
                                        shiny::tableOutput("data"),

                                        shiny::tabsetPanel(
                                          id = "navbar",

                                          type = "pills",
                                          shiny::tabPanel(title_0 = "All",

                                                          value = c("All"),
                                                          shiny::p("All"),
                                                          # P-values
                                                          shiny::img(src='./cor90_p.png', align="left", width = "28.8%"),
                                                          shiny::img(src='./ax90_p.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag90_p.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Adj P-values
                                                          shiny::img(src='./cor90_adjp.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax90_adjp.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag90_adjp.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.2%"),
                                                          # T-values
                                                          shiny::img(src='./cor90_t.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax90_t.png', align="left",width = "24%"),
                                                          shiny::img(src='./sag90_t.png', align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%"),
                                                          # Atlas
                                                          shiny::img(src='./atlas.107.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./atlas.106.png', align="left",width = "24%"),
                                                          shiny::img(src='./atlas.108.png', align="left",width = "34.5%")),
                                          shiny::tabPanel(title_0 = "P-Values",
                                                          value = c("P-Values"),
                                                          shiny::p("P-Values"),
                                                          # P-values
                                                          shiny::img(src='./cor90_p.png', align="left", width = "28.8%"),
                                                          shiny::img(src='./ax90_p.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag90_p.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Adjusted P-Values",
                                                          value = c("Adjusted P-Values"),
                                                          shiny::p("Adjusted P-Values"),
                                                          shiny::img(src='./cor90_adjp.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax90_adjp.png', align="left", width = "24%"),
                                                          shiny::img(src='./sag90_adjp.png', align="left", width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%")),
                                          shiny::tabPanel(title_0 = "T-Values",
                                                          value = c("T-Values"),
                                                          shiny::p("T-Values"),
                                                          shiny::img(src='./cor90_t.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./ax90_t.png', align="left",width = "24%"),
                                                          shiny::img(src='./sag90_t.png', align="left",width = "34.5%"),
                                                          shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%")),
                                          shiny::tabPanel(title_0 = "Atlas",
                                                          value = c("Atlas"),
                                                          shiny::p("Atlas"),
                                                          shiny::img(src='./atlas.107.png', align="left",width = "28.8%"),
                                                          shiny::img(src='./atlas.106.png', align="left",width = "24%"),
                                                         shiny::img(src='./atlas.108.png', align="left",width = "34.5%")))

                        )
                      )
                    )
                  )

                }
              )
  )

testR6 <- BssRmdVolumeOutput$new(outdir = "/Users/sarapesavento/Desktop/tbm_anova/")
testR6$save_out(bss_data, bss_model)
