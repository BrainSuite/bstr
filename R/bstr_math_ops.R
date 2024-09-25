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

#' log10 transform
#' @param values numeric vector
#' @return numeric vector containing the log10 transformed `values`
#' @export
log10_transform <- function(values) {
  eps <- .Machine$double.eps
  sgn <- sign(values)
  sgn[sgn == 0] <- 1
  values[abs(values) < eps] <- eps
  logvalues <- - 1*sgn*log10(abs(values))
  return(logvalues)
}

# TBD
image_to_shape <- function(imagefile, shapefile, outputshapefile, resample) {

}

#' This function calculates the sign of t values. Differently from the sign() function, sign_tvalues(0) = 1
#' @param tvalues numeric vector
sign_tvalues <- function(tvalues) {
  tvalues_sign <- sign(tvalues)
  tvalues_sign[tvalues_sign == 0] <- 1
  return(tvalues_sign)
}
