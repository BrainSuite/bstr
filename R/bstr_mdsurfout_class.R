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

#' R6 derived class for Rmd surface output functionality
#' @export
BstrRmdVolumeOutput <-
  R6::R6Class("BstrRmdSurfaceOutput",
              inherit = BstrRmdOutput,
              public = list(
                initialize = function(outdir = NULL) {
                  super$initialize(outdir)
                },
                save_out = function(bstr_data, bstr_model) {

                }
              ),
              private = list(
                render_overlay = function() {

                },
                render_atlas = function() {

                },
                render_html = function() {

                }
              )
  )
