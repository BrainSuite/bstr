library(methods)

"usage:
bss_roi.R <mspec> <outdir>

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
  bssr::bss_roi(subjdir = mspec$subjdir, csv = mspec$csv, mspec$roiid, mspec$roimeasure,
          main_effect = mspec$main_effect, covariates = mspec$covariates, mspec_file = mspec$mspec_file,
          outdir = opt$outdir)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

