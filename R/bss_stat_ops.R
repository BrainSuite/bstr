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

  corr_result <- corr_vec(bss_data@data_array, bss_data@demographics[[corr_var]])
  corr_coeff <- corr_result$corr_coeff
  bss_model@tvalues <- corr_result$tvalues
  bss_model@pvalues <- corr_result$pvalues

  bss_model@pvalues <- sign(corr_coeff)*bss_model@pvalues
  bss_model@pvalues[is.na(bss_model@pvalues)] <- 1 # Set the p-values with the NA correlations to 1
  bss_model@corr_values <- corr_coeff
  bss_model@pvalues_adjusted <- p.adjust(bss_model@pvalues, 'BH')
  bss_model@corr_values[abs(bss_model@pvalues) >= 0.05] <- 0
  message('Done.')
  return(bss_model)
}

#' Vectorized correlation
corr_vec <- function(X, Y) {

  N <- length(Y)
  X_dev <- sweep(X, 2, colMeans(X))
  Y_dev <- Y - mean(Y)
  corr_coeff <- as.numeric((Y_dev %*% X_dev)/sqrt(colSums(X_dev^2)*sum(Y_dev^2)))
  tvalues <- corr_coeff * sqrt((N-2)/(1-corr_coeff^2 + .Machine$double.eps))
  pvalues <- 2*pt(abs(tvalues), N-2, lower.tail = FALSE)
  return(list("tvalues"=tvalues, "pvalues"=pvalues, "corr_coeff"=corr_coeff))
}

#' Welch's t-test
#' @export
bss_ttest <- function(group_var, bss_data, paired = FALSE) {

  if (paired == FALSE)
    bss_model <- new("BssModel", model_type="unpairedttest", group_var = group_var,
                     demographics = bss_data@demographics, mspec_file="")
  else
    bss_model <- new("BssModel", model_type="pairedttest", group_var = group_var,
                     demographics = bss_data@demographics, mspec_file="")

  group1 <- levels(bss_data@demographics[[group_var]])[1]
  group2 <- levels(bss_data@demographics[[group_var]])[2]
  idx_group1 <- which(bss_data@demographics[[group_var]] == group1)
  idx_group2 <- which(bss_data@demographics[[group_var]] == group2)

  message(sprintf("The 2 groups are %s and %s", group1, group2), appendLF = TRUE)
  message('Running t-tests...', appendLF = FALSE)

  test_result <- ttest_vec(bss_data@data_array[idx_group1,], bss_data@data_array[idx_group2,], paired)
  pvalues <- test_result$pvalues
  tvalues <- test_result$tvalues

  pvalues[is.na(pvalues)] <- 1
  pvalues <- pvalues*sign(tvalues)
  bss_model@pvalues <- pvalues
  bss_model@tvalues <- tvalues
  bss_model@tvalues[abs(pvalues) >= 0.05] <- 0
  bss_model@pvalues_adjusted <- p.adjust(abs(bss_model@pvalues), 'BH')

  message('Done.')
  return(bss_model)
}

#' Vectorized Welch's t-test
ttest_vec <- function(X1, X2, paired=FALSE) {

  n1 <- dim(X1)[1]
  n2 <- dim(X2)[1]

  if (paired == TRUE) {
    D = X1 - X2
    D_mean <- colMeans(D)
    D_dev <- sweep(D, 2, D_mean)
    s1 <- sqrt(colSums(D_dev^2)/(n1-1))
    tvalues <- D_mean/(s1/sqrt(n1))
    pvalues <- 2*pt(abs(tvalues),n1-1, lower.tail = FALSE)
  }
  else {
    X1_mean <- colMeans(X1)
    X2_mean <- colMeans(X2)

    X1_dev <- sweep(X1, 2, X1_mean)
    X2_dev <- sweep(X2, 2, X2_mean)

    SS1 <- colSums(X1_dev^2) # sum of squares of group 1
    SS2 <- colSums(X2_dev^2) # sum of squares of group 2

    s1_sq <- SS1/(n1-1)
    s2_sq <- SS2/(n2-1)

    se_diff <- sqrt(s1_sq/n1 + s2_sq/n2)
    tvalues <- (X1_mean - X2_mean + .Machine$double.eps)/(se_diff + .Machine$double.eps)
    # Calculate the degrees of freedom using the Welch–Satterthwaite approximation
    deg <- (s1_sq/n1 + s2_sq/n2)^2/(s1_sq^2/(n1^2*(n1-1)) + s2_sq^2/(n2^2*(n2-1)))
    pvalues <- 2*pt(abs(tvalues), deg, lower.tail = FALSE)
  }
  return(list("tvalues"=tvalues, "pvalues"=pvalues))
}
