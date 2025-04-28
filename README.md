# README
bstr: BrainSuite Statistics Toolbox in R

Copyright (C) 2025 The Regents of the University of California

Created by Shantanu H. Joshi, Yeun Kim, Kayla A. Schroeder, and David W. Shattuck

bstr is licensed under an GPLv2-only license (https://spdx.org/licenses/GPL-2.0-only.html).
Please see the enclosed LICENSE file for more details.

---

The BrainSuite Statistics toolbox in R (bstr) is a software package developed in R that performs statistical analysis of population-level neuroimaging data processed using BrainSuite [1]. Specifically, it provides statistical tools for conducting cortical thickness analysis, tensor based morphometry, and analysis of diffusion measures.

## Installation

### Prerequisites
* Install [R](https://cran.r-project.org), [RStudio](https://posit.co/products/open-source/rstudio/), and [Rools](https://cran.r-project.org/bin/windows/Rtools) (users on MS Windows only)
* Install [BrainSuite](https://brainsuite.org)

### Install from GitHub (recommended)
* Open RStudio and enter the following commands to install the latest version of [bstr](https://github.com/BrainSuite/bstr). 
```
install.packages('remotes')
remotes::install_github("BrainSuite/bstr/")
```

### Check your installation
Type
```
library(bstr)
get_brainsuite_install_path()
```
This should display the BrainSuite installation path.

### Statistical analysis using *bstr*
Please follow detailed instructions for data preparation, workflows, and usage examples at [brainsuite.org/bstr](http://brainsuite.org/bstr/).


## Types of Analyses 
Bstr performs statistical analysis on the outputs of the BrainSuite structural workflow, which performs cortical surface extraction [1], alignment to a reference atlas using surface-constrained volumetric registration (SVReg) [2], and, optionally, processing of diffusion MRI data using the BrainSuite diffusion pipeline (BDP) [3]. SVReg performs surface registration of triangular meshes based on curvature and volumetric registration based on image intensities. BDP performs distortion correction, alignment of diffusion MRI to T1-weighted MRI, and fitting of various diffusion models to the corrected diffusion data. Bstr is used to perform population-level statistical analysis of various neuroimaging measures produced by these components. Statistical analysis of voxel-wise and surface-based data is performed in the common coordinate space of the atlas by resampling the data from subject coordinates to a the atlas space using SVReg.

Bstr supports the following analysis methods:

* tensor based morphometry (TBM) analysis of voxel-wise magnitudes of the 3D deformation fields of MRI images registered to the atlas
* cortical surface analysis (SBA) of the vertex-wise thickness in the atlas space
* diffusion parameter maps analysis (DBA) of fractional anisotropy, mean diffusivity, radial diffusivity
* region of interest (ROI)-based analysis of average gray matter thickness, surface area, and gray matter volume within cortical ROIs
* correction for multiple comparisons using false discovery rate (FDR) or permutation testing methods

Bstr is cross-platform and is available on macOS, Windows,and Linux based systems (all platforms with R support). Bstr is distributed under an open source license (GPLv2-only). Bstr supports functionality for automated report generation to visualize statistical results using R-shiny and R markdown. The volumetric analysis report contains the cluster table, visualizations of clusters on image slices, and shows both the unadjusted and the adjusted versions of p-values and t statistics, respectively. The ROI analysis report shows the demographic spreadsheet, automatic bar plots for ANOVA and regressions, and scatter plot for correlation analyses. Bstr also exports an R markdown report that contains reproducible R commands in both the Rmd file and in the html document [4]. This enables complete reproducibility of statistical results and only requires packaging the R markdown file along with the data.

---

### References
1. Shattuck DW et al. (2002) BrainSuite: An Automated Cortical Surface Identication Tool Medical Image Analysis, 8(2):129-142.
2. Joshi AA et al. (2007) Surface-Constrained Volumetric Brain Registration Using Harmonic Mappings IEEE Trans. on Medical Imaging 26(12):1657-1669.
3. Bhushan C et al. (2015) Co-registration and distortion correction of diffusion and anatomical images based on inverse contrast normalization. Neuroimage (7):115:269-80.
4. Xie Y (2017) Dynamic Documents with R and knitr. Chapman and Hall/CRC.
