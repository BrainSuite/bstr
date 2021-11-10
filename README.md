# README
bssr: BrainSuite Statistics Toolbox in R

Copyright (C) 2021 The Regents of the University of California

Created by Shantanu H. Joshi, Yeun Kim, Kayla A. Schroeder, and David W. Shattuck

bssr is licensed under an GPLv2-only license (https://spdx.org/licenses/GPL-2.0-only.html).
Please see the enclosed LICENSE file for more details.

---

The BrainSuite Statistics toolbox in R (bssr) is a software package developed in R that performs statistical analysis of population-level neuroimaging data processed using BrainSuite [1]. Specifically, it provides statistical tools for conducting cortical thickness analysis, tensor based morphometry, and analysis of diffusion measures.

## bssr Installation
For more detailed installation instructions, usage examples, or to check for updated versions of bssr, please visit the bssr website: http://brainsuite.org/bssr.

### Prerequisites
* Ensure BrainSuite is installed on your computer: <http://brainsuite.org/download>
* Ensure R is installed on your computer: https://cran.r-project.org/
* Install RStudio: https://www.rstudio.com/products/rstudio/#Desktop
* Windows only: Install Rtools, available at https://cran.r-project.org/bin/windows/Rtools/

### Steps for installation
* Open RStudio and enter the following two commands to install **bssr version 0.3.2**.
* You can install from the bssr_0.3.2.tar.gz file in the BrainSuite21a folder, or you can pull directly from the BrainSuite server.
* You do not need to untar bssr_0.3.2.tar.gz -- the installer needs it in this format.
* The following examples assume you have installed BrainSuite21a in the default location for your operating system. If you have used a different path, you will need to specify it in the installation call.

### Mac OS - Install from folder
```
install.packages('devtools')
devtools::install_local('/Applications/BrainSuite21a/bssr/bssr_0.3.2.tar.gz')
```

### Windows - Install from folder
Note that on Windows, you will need to use double backslashes (\\) in the path because backslash is an escape character. You can also replace the backslashes with forward slashes.
```
install.packages('devtools')
devtools::install_local('C:\\Program Files\\BrainSuite21a\\bssr\\bssr_0.3.2.tar.gz')
```
or
```
install.packages('devtools')
devtools::install_local('C:/Program Files/BrainSuite21a/bssr/bssr_0.3.2.tar.gz')
```


### Linux - Install from folder
```
install.packages('devtools')
devtools::install_local('/opt/BrainSuite21a/bssr/bssr_0.3.2.tar.gz')
```


### Online Installation
```
install.packages('devtools')
devtools::install_url('http://brainsuite.org/wp-content/uploads/2021/05/bssr_0.3.2.tar.gz')
```

### Check your installation
Type
```
library(bssr)
```
then 
```
get_brainsuite_install_path()
```
This should display the BrainSuite installation path.


## Methods
Bssr performs statistical analysis on the outputs of the BrainSuite structural workflow, which performs cortical surface extraction [1], alignment to a reference atlas using surface-constrained volumetric registration (SVReg) [2], and, optionally, processing of diffusion MRI data using the BrainSuite diffusion pipeline (BDP) [3]. SVReg performs surface registration of triangular meshes based on curvature and volumetric registration based on image intensities. BDP performs distortion correction, alignment of diffusion MRI to T1-weighted MRI, and fitting of various diffusion models to the corrected diffusion data. Bssr is used to perform population-level statistical analysis of various neuroimaging measures produced by these components. Statistical analysis of voxel-wise and surface-based data is performed in the common coordinate space of the atlas by resampling the data from subject coordinates to a the atlas space using SVReg.

Bssr supports the following analysis methods:

* tensor based morphometry (TBM) analysis of voxel-wise magnitudes of the 3D deformation fields of MRI images registered to the atlas
* cortical surface analysis of the vertex-wise thickness in the atlas space
* diffusion parameter maps analysis (e.g., fractional anisotropy, mean diffusivity, radial diffusivity) 
* region of interest (ROI)-based analysis of average gray matter thickness, surface area, and gray matter volume within cortical ROIs
* correction for multiple comparisons using false discovery rate (FDR) or permutation testing methods

Bssr is cross-platform and is available on macOS, Windows,and Linux based systems (all platforms with R support). Bssr is distributed under an open source license (GPLv2-only). Bssr supports functionality for automated report generation to visualize statistical results using R-shiny and R markdown. The volumetric analysis report contains the cluster table, visualizations of clusters on image slices, and shows both the unadjusted and the adjusted versions of p-values and t statistics, respectively. The ROI analysis report shows the demographic spreadsheet, automatic bar plots for ANOVA and regressions, and scatter plot for correlation analyses. Bssr also exports an R markdown report that contains reproducible R commands in both the Rmd file and in the html document [4]. This enables complete reproducibility of statistical results and only requires packaging the R markdown file along with the data.

---

### References
1. Shattuck DW et al. (2002) BrainSuite: An Automated Cortical Surface Identication Tool Medical Image Analysis, 8(2):129-142.
2. Joshi AA et al. (2007) Surface-Constrained Volumetric Brain Registration Using Harmonic Mappings IEEE Trans. on Medical Imaging 26(12):1657-1669.
3. Bhushan C et al. (2015) Co-registration and distortion correction of diffusion and anatomical images based on inverse contrast normalization. Neuroimage (7):115:269-80.
4. Xie Y (2017) Dynamic Documents with R and knitr. Chapman and Hall/CRC.
