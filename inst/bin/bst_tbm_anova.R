#!/usr/bin/env Rscript
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
library(methods)
library(bstr)

"usage:
bss_tbm_anova.R --subjdir=<subjdir> --csv=<csv> --main_effect=<main_effect> --covariates=<covariates> --odir=<odir> [--smooth=<smooth>]
-h --help   show this help message and exit
options:
--subjdir subjdir subject directory
--csv csv demographics csv
--main_effect main_effect main effect (single variable from the demographcis csv)
--covariates covariates covariates (multiple variables separated by a + sign, no spaces)
--odir odir output directory
--smooth smooth smoothing level

" -> doc

opt <- docopt::docopt(doc)
bss_data <- load_bss_data(type="tbm", subjdir = opt$subjdir,
                          csv = opt$csv, smooth = as.numeric(opt$smooth))
bss_model <- bss_anova(main_effect = opt$main_effect, covariates = opt$covariates, bss_data = bss_data)
save_bss_out(bss_data, bss_model, outdir=opt$odir)
