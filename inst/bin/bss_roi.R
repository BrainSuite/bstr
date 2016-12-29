library(methods)
library(bssr)

"usage:
bss_roi.R <mspec> <outdir>

options:
--mspec=<mspec> ini file for model specification
--outdir=<outdir> output directory

-h --help   show this help message and exit
" -> doc

opt <- docopt::docopt(doc)



bss_roi <- function(subjdir, csv, roiid, roimeas, main_effect, covariates, outdir) {

  start.time <- Sys.time()

  bss_roi_data <- new("BssROIData", subjdir, csv)
  bss_roi_data <- load_data(bss_roi_data, roiid = roiid, roimeas = roimeas, outdir)
  bss_model <- new("BssModel", model_type="anova", main_effect = main_effect, covariates,
                   demographics = bss_roi_data@demographics)
  bss_model <- run(bss_model, bss_roi_data)
  bss_output <- new("BssROIOutput", outdir)
  bss_output <- save_out(bss_output, bss_roi_data, bss_model)

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
}

if ( !is.null(opt$mspec) ) {
  # modelspec.ini is provided. Parse it and get the options
  mspec <- read_modelspec(opt$mspec)
  bss_roi(subjdir = mspec$subjdir, csv = mspec$csv, mspec$roiid, mspec$roimeasure,
          main_effect = mspec$main_effect, covariates = mspec$covariates, outdir = opt$outdir)
}


