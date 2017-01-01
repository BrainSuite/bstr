library(methods)

"usage:
bss_roi.R <mspec> <outdir>

options:
--mspec=<mspec> ini file for model specification
--outdir=<outdir> output directory

-h --help   show this help message and exit
" -> doc

opt <- docopt::docopt(doc)

if ( !is.null(opt$mspec) ) {
  # modelspec.ini is provided. Parse it and get the options
  mspec <- bssr::read_modelspec(opt$mspec)
  bssr::bss_roi(subjdir = mspec$subjdir, csv = mspec$csv, mspec$roiid, mspec$roimeasure,
          main_effect = mspec$main_effect, covariates = mspec$covariates, outdir = opt$outdir)
}


