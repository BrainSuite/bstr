library(methods)
library(bssr)

"usage:
bss_tbm.R <mspec> <outdir>

options:
--mspec=<mspec> ini file for model specification
--outdir=<outdir> output directory

-h --help   show this help message and exit
" -> doc

opt <- docopt::docopt(doc)

start.time <- Sys.time()

bss_tbm <- function(subjdir, csv, atlas, maskfile, sigma_smooth, main_effect, covariates, outdir) {

  bss_tbm_data <- new("BssTBMData", subjdir, csv)
  bss_tbm_data <- load_data(bss_tbm_data, atlas_filename = atlas, maskfile = maskfile, smooth = sigma_smooth)
  bss_model <- new("BssModel", model_type="anova", main_effect = main_effect, covariates,
                   demographics = bss_tbm_data@demographics)
  bss_model <- run(bss_model, bss_tbm_data)
  bss_output <- new("BssTBMOutput", outdir)
  bss_output <- save_out(bss_output, bss_tbm_data, bss_model)

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
}

if ( !is.null(opt$mspec) ) {
  # modelspec.ini is provided. Parse it and get the options
  mspec <- read_modelspec(opt$mspec)
  bss_tbm(subjdir = mspec$subjdir, csv = mspec$csv, atlas = mspec$nii_atlas, maskfile = mspec$maskfile,
          sigma_smooth = mspec$smooth, main_effect = mspec$main_effect, covariates = mspec$covariates, outdir = opt$outdir)
}


