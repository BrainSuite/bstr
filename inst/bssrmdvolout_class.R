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
                                     initialize = function(outdir = NULL) {
                                         super$initialize(outdir)
                                       },
                                     save_out = function(bss_data, bss_model) {
                       
                                        }
                                   ),
                                 private = list(
                                       render_overlay = function(voxelcoord,filePath,overlay,folderPath,fileName) {
                                         intern = FALSE
                                         ignore.stdout = FALSE
                                         ignore.stderr = FALSE
                                         wait = TRUE
                                         input = NULL
                                         
                                         view1 <- sprintf("volblend -r ",filePath[3]," --view 1 --slice ",voxelcoord[3]," --flop -o atlas.",voxelcoord,".p.png -c jet -i ",overlay," -o ",folderPath, fileName,"cor_1.png")
                                         view2 <- sprintf("volblend -r ",filePath[2]," --view 2 --slice ",voxelcoord[2]," --flop -o atlas.",voxelcoord,".p.png -c jet -i ",overlay," -o ",folderPath, fileName,"sag_2.png")
                                         view3 <- sprintf("volblend -r ",filePath[1]," --view 3 --slice ",voxelcoord[1]," --flop -o atlas.",voxelcoord,".p.png -c jet -i ",overlay," -o ",folderPath, fileName,"ax_3.png")
                                         
                                         return(c(system(view1, intern, ignore.stdout, ignore.stderr, wait, input), system(view2, intern, ignore.stdout, ignore.stderr, wait, input), system(view3, intern, ignore.stdout, ignore.stderr, wait, input)))
                                       },
                                       render_atlas = function(voxelcoord,filePath,folderPath,atlasName) {
                                         intern = FALSE
                                         ignore.stdout = FALSE
                                         ignore.stderr = FALSE
                                         wait = TRUE
                                         input = NULL
                                         
                                         view1 <- sprintf("volblend -i ",filePath," --view 1 --slice ",voxelcoord[3]," --flop -o ", folderPath, atlasName,"cor_1.png")
                                         view2 <- sprintf("volblend -i ",filePath," --view 2 --slice ",voxelcoord[2]," --flop -o ", folderPath, atlasName,"sag_2.png")
                                         view3 <- sprintf("volblend -i ",filePath," --view 3 --slice ",voxelcoord[1]," --flop -o ", folderPath, atlasName,"ax_3.png")
                                         return(c(system(view1, intern, ignore.stdout, ignore.stderr, wait, input), system(view2, intern, ignore.stdout, ignore.stderr, wait, input), system(view3, intern, ignore.stdout, ignore.stderr, wait, input)))
                                         },
                                       render_table = function() {
                                         t <- c("#test3","#test4","#test5","#test6","#test7")
                                         Cluster <- 1:5
                                         table <- data.frame(Cluster = 1:5, Vol_Size = c(33,22.3, 21, 25, 30), MNI_coord = c(3,3,4,5,3), T_val =c(8,7.2,6,9.8,7.8))
                                         table$Cluster <- sprintf("[", table$Cluster, "](", t, ")")
                                         knitr::kable(table[1:4], align=c(rep('l', 4)))
                                       },
                                       render_vol = function() {
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
                                                                                 shiny::img(src='./Pview1.png', align="left", width = "28.8%"),
                                                                                 shiny::img(src='./Pview3.png', align="left", width = "24%"),
                                                                                 shiny::img(src='./Pview2.png', align="left", width = "34.5%"),
                                                                                 shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%"),
                                                                                 # Adj P-values 
                                                                                 shiny::img(src='./AdjPvalue1.png', align="left",width = "28.8%"),
                                                                                 shiny::img(src='./AdjPvalue3.png', align="left", width = "24%"),
                                                                                 shiny::img(src='./AdjPvalue2.png', align="left", width = "34.5%"),
                                                                                 shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%"),
                                                                                 # T-values
                                                                                 shiny::img(src='./Tvalue1.png', align="left",width = "28.8%"),
                                                                                 shiny::img(src='./Tvalue3.png', align="left",width = "24%"),
                                                                                 shiny::img(src='./Tvalue2.png', align="left",width = "34.5%"),
                                                                                 shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%"),
                                                                                 # Atlas
                                                                                 shiny::img(src='./atlas.107.png', align="left",width = "28.8%"),
                                                                                 shiny::img(src='./atlas.106.png', align="left",width = "24%"),
                                                                                 shiny::img(src='./atlas.108.png', align="left",width = "34.5%")),
                                                                 shiny::tabPanel(title_0 = "P-Values",
                                                                                 value = c("P-Values"),
                                                                                 shiny::p("P-Values"),
                                                                                 # P-values
                                                                                 shiny::img(src='./Pview1.png', align="left", width = "28.8%"),
                                                                                 shiny::img(src='./Pview3.png', align="left", width = "24%"),
                                                                                 shiny::img(src='./Pview2.png', align="left", width = "34.5%"),
                                                                                 shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_cbar.pdf', align="left", width = "11%")),
                                                                 shiny::tabPanel(title_0 = "Adjusted P-Values",
                                                                                 value = c("Adjusted P-Values"),
                                                                                 shiny::p("Adjusted P-Values"),
                                                                                 shiny::img(src='./AdjPvalue1.png', align="left",width = "28.8%"),
                                                                                 shiny::img(src='./AdjPvalue3.png', align="left", width = "24%"),
                                                                                 shiny::img(src='./AdjPvalue2.png', align="left", width = "34.5%"),
                                                                                 shiny::img(src='./bss_anova_age_mri.bfc.nii_log_pvalues_adjusted_cbar.pdf', align="left", width = "11.5%")),
                                                                 shiny::tabPanel(title_0 = "T-Values",
                                                                                 value = c("T-Values"),
                                                                                 shiny::p("T-Values"),
                                                                                 shiny::img(src='./Tvalue1.png', align="left",width = "28.8%"),
                                                                                 shiny::img(src='./Tvalue3.png', align="left",width = "24%"),
                                                                                 shiny::img(src='./Tvalue2.png', align="left",width = "34.5%"),
                                                                                 shiny::img(src='./bss_anova_age_mri.bfc.nii_tvalues_cbar.pdf', align="left", width = "11%")),
                                                                 shiny::tabPanel(title_0 = "Atlas",
                                                                                 value = c("Atlas"),
                                                                                 shiny::p("Atlas"),
                                                                                 shiny::img(src='./atlas.107.png', align="left",width = "28.8%"),
                                                                                 shiny::img(src='./atlas.106.png', align="left",width = "24%"),
                                                                                 shiny::img(src='./atlas.108.png', align="left",width = "34.5%")))),
                                               shiny::tabPanel(title = shiny::h4("Cluster 2"),
                                                               
                                                               value = "tab2",
                                                               shiny::h2("Cluster 2")
                                               ),
                                               shiny::tabPanel(title = shiny::h4("Cluster 3"),
                                                               value = "tab3",
                                                               shiny::h2("Cluster 3:"),
                                                               shiny::mainPanel("hi")
                                                               
                                               ),
                                             ),
                                           ),
                                         )}         
                                         
                                 )
      )
                  
                  