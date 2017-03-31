#' Bss statistical functions for regression (lm), correlation etc.

#' linear regression
#'
#' @export
bss_lm <- function(main_effect="", covariates="", bss_data) {

  if (class(bss_data) == "BssROIData") {
    return(bss_roi_lm(main_effect = main_effect, covariates = covariates, bss_data = bss_data))
  }

  # Check the model type and call the appropriate method
  bss_model <- new("BssModel", model_type="bss_lm", main_effect = main_effect, covariates = covariates,
                   demographics = bss_data@demographics, mspec_file="")
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  Xtemp <- solve(t(bss_model@X_design_full) %*% bss_model@X_design_full) %*% t(bss_model@X_design_full) # Pre Hat matrix
  beta_full <- Xtemp %*% bss_data@data_array  # beta coefficients
  y_full <- bss_model@X_design_full %*% beta_full  # Predicted response
  RSS_full <- colSums((bss_data@data_array - y_full)^2)

  Xtemp <- solve(t(bss_model@X_design_null) %*% bss_model@X_design_null) %*% t(bss_model@X_design_null) # Pre Hat matrix
  beta_null <- Xtemp %*% bss_data@data_array  # beta coefficients
  y_null <- bss_model@X_design_null %*% beta_null  # Predicted response
  RSS_null <- colSums((bss_data@data_array - y_null)^2)

  N <- nrow(bss_data@data_array)
  Fstat <- (RSS_null - RSS_full)/RSS_full * (N - bss_model@Npfull - 1)/(bss_model@Npfull - bss_model@Npnull)  # F statistic
  model_unique_idx <- which(bss_model@unique %in% bss_model@fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept

  se_full_unique <- sqrt(diag(solve(t(bss_model@X_design_full) %*% bss_model@X_design_full)))

  se_full_unique <- sqrt(diag(solve(t(bss_model@X_design_full) %*% bss_model@X_design_full)))[model_unique_idx] *
    sqrt(RSS_full / (N - bss_model@Npfull - 1))

  tvalue_sign <- (beta_full[model_unique_idx, ] + .Machine$double.eps)/(abs(beta_full[model_unique_idx, ]) + .Machine$double.eps)

  pvalues <- 1 - pf(Fstat, bss_model@Npfull - bss_model@Npnull, N - bss_model@Npfull - 1)

  pvalues[is.nan(pvalues)] <- 1

  pvalues <- pvalues*tvalue_sign
  tvalues <- beta_full[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  bss_model@pvalues <- pvalues
  bss_model@tvalues <- tvalues
  bss_model@tvalues[abs(pvalues) >= 0.05] <- 0
  bss_model@pvalues_adjusted <- p.adjust(bss_model@pvalues, 'BH')
  message('Done.')
  return(bss_model)
}

bss_roi_lm <- function(main_effect="", covariates="", bss_data=bss_data) {
  # Check the model type and call the appropriate method
  bss_model <- new("BssModel", model_type="bss_lm", main_effect = main_effect, covariates = covariates,
                   demographics = bss_data@demographics, mspec_file="")
  message('Running the statistical model. This may take a while...', appendLF = FALSE)

  bss_data@demographics[paste('ROI_', as.character(bss_data@roiid), sep = '' )]

  cmd1 <- sprintf("lm_full <- lm(%s, data = bss_data@demographics)",
                  paste('ROI_', as.character(bss_data@roiid), ' ~ ', bss_model@fullmodel, sep = ''))
  cmd2 <- sprintf("lm_null <- lm(%s, data = bss_data@demographics)",
                  paste('ROI_', as.character(bss_data@roiid), ' ~ ', bss_model@nullmodel, sep = ''))
  cmd3 <- "pander::pander(anova(lm_full, lm_null))"

  stats_commands <- c(cmd1, cmd2, cmd3)

  for (cmd in stats_commands) {
    eval(parse(text = cmd))
  }
  bss_model@stats_commands <- stats_commands

  return(bss_model)

}

#' Correlation
#' @export
bss_corr <- function(corr_var, bss_data) {

  message('Running correlations...', appendLF = FALSE)
  bss_model <- new("BssModel", model_type="bss_corr", corr_var = corr_var,
                   demographics = bss_data@demographics, mspec_file="")
  X <- sweep(bss_data@data_array, 2, colMeans(bss_data@data_array))
  Y <- bss_data@demographics[[corr_var]] - mean(bss_data@demographics[[corr_var]])
  corr_coeff  <- as.numeric((Y %*% X)/sqrt(colSums(X^2)*sum(Y^2)))
  num_subjects <- length(Y)
  tvalues <- corr_coeff * sqrt((num_subjects-2)/(1-corr_coeff^2 + .Machine$double.eps))
  bss_model@tvalues <- tvalues
  bss_model@pvalues <- 1 - pt(abs(tvalues), num_subjects-2)

  bss_model@pvalues <- sign(corr_coeff)*bss_model@pvalues
  bss_model@pvalues[is.na(bss_model@pvalues)] <- 1 # Set the p-values with the NA correlations to 1
  bss_model@corr_values <- corr_coeff
  bss_model@pvalues_adjusted <- p.adjust(bss_model@pvalues, 'BH')
  bss_model@corr_values[abs(bss_model@pvalues) >= 0.05] <- 0
  message('Done.')
  return(bss_model)
}



