library(methods)
library(bssr)

"usage:
bss_cbm.R <mspec> <outdir>

options:
--mspec=<mspec> ini file for model specification
--outdir=<outdir> output directory

-h --help   show this help message and exit
" -> doc

opt <- docopt::docopt(doc)

start.time <- Sys.time()

bss_cbm <- function(subjdir, csv, lh_surf_atlas, rh_surf_atlas, sigma_smooth, main_effect, covariates, outdir) {

  # Left hemisphere
  bss_cbm_data <- new("BssCBMData", subjdir, csv)
  bss_cbm_data <- load_data(bss_cbm_data, atlas_filename = lh_surf_atlas,
                            hemi = 'left', smooth=sigma_smooth)
  bss_model <- new("BssModel", model_type="anova", main_effect = main_effect, covariates,
                   demographics = bss_cbm_data@demographics)
  bss_model <- run(bss_model, bss_cbm_data)
  bss_output <- new("BssCBMOutput", outdir)
  bss_output <- save_out(bss_output, bss_cbm_data, bss_model)

  # Right hemisphere
  bss_cbm_data <- new("BssCBMData", subjdir, csv)
  bss_cbm_data <- load_data(bss_cbm_data, atlas_filename = rh_surf_atlas,
                            hemi = 'right', smooth=sigma_smooth)
  bss_model <- new("BssModel", model_type="anova", main_effect = main_effect, covariates,
                   demographics = bss_cbm_data@demographics)
  bss_model <- run(bss_model, bss_cbm_data)
  bss_output <- new("BssCBMOutput", outdir)
  bss_output <- save_out(bss_output, bss_cbm_data, bss_model)


  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
}

if ( !is.null(opt$mspec) ) {
  # modelspec.ini is provided. Parse it and get the options
  mspec <- read_modelspec(opt$mspec)
  bss_cbm(subjdir = mspec$subjdir, csv = mspec$csv, lh_surf_atlas = mspec$lh_surf_atlas, rh_surf_atlas = mspec$rh_surf_atlas,
          sigma_smooth = mspec$smooth, main_effect = mspec$main_effect, covariates = mspec$covariates, outdir = opt$outdir)
}


