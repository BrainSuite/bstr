# BrainSuite Statistics Toolbox in R (bstr)
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

#' R6 super class for Rmd output functionality
#' @export
BssRmdOutput <-
  R6::R6Class("BssRmdOutput",
              public = list(
                #' @field outdir output directory
                outdir = NULL,
                #' @description initialize function
                #' @param outdir path to the output directory
                initialize = function(outdir = NULL){
                  self$outdir <- outdir
                },
                #' @description Save results as Rmd
                #' @param bss_data object of type \code{BssData}
                #' @param bss_model object of type \code{BssModel}
                save_out = function(bss_data, bss_model){

                },
                #' @description finalize function
                finalize = function(){

                }
              )
  )



