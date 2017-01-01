library(methods)

"usage:
bss_tbm.R <mspec> <outdir>

options:
--mspec=<mspec> ini file for model specification
--outdir=<outdir> output directory

-h --help   show this help message and exit
" -> doc

opt <- docopt::docopt(doc)

start.time <- Sys.time()

if ( !is.null(opt$mspec) ) {
  # modelspec.ini is provided. Parse it and get the options
  mspec <- bssr::read_modelspec(opt$mspec)
  bssr::bss_tbm(subjdir = mspec$subjdir, csv = mspec$csv, atlas = mspec$nii_atlas, maskfile = mspec$maskfile,
          sigma_smooth = mspec$smooth, main_effect = mspec$main_effect, covariates = mspec$covariates, outdir = opt$outdir)
}


