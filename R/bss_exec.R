#' Bss executable functions for cbm, tbm, roi etc.

#' @export
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

#' @export
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

#' @export
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

