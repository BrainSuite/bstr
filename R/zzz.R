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

#' @importFrom grDevices col2rgb colorRamp colorRampPalette dev.new
#' @importFrom graphics axis plot rect
#' @importFrom methods .valueClassTest new as
#' @importFrom stats formula model.matrix p.adjust pf pt p.adjust.methods
#' @importFrom utils read.csv read.table write.csv setTxtProgressBar txtProgressBar
#' @importFrom DT datatable
#' @importFrom pander pander
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom bit as.bit
#'
NULL

.onLoad <- function(libname, pkgname) {

  if(!is_brainsute_installed(quiet = FALSE, raise_error = FALSE)) {
    setup(quiet = FALSE, raise_error = FALSE)
    invisible()
  }
}
