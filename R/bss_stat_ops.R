# BrainSuite Statistics Toolbox in R (bstr)
# Copyright (C) 2017 The Regents of the University of California
# Creator: Shantanu H. Joshi, Department of Neurology, Ahmanson Lovelace Brain Mapping Center, UCLA
#
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation; version 2.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License version 2 for more details.
#
# You should have received a copy of the GNU General Public License along with this program;
# if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

#' Perform analysis of variance (ANOVA) for brain imaging data.
#'
#' This function accepts a `main_effect` and a set of covariates (using the R formula notation) and uses an
#' F-test to compare the full model including the `main_effect + covariates` with the reduced (null) model
#' that only includes the `covariates`.
#'
#' Slightly different from the standard R anova function, `bstr_anova` currently does not directly accept the
#' results from `lm_vec`. This could be accomodated in the future versions.

#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a `+` operator similar to an R formula.
#' @param  bstr_data Object of type [BstrData()]
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See [bstr_p_adjust()] for valid values.
#' @seealso [lm_vec()] for linear regression, [bstr_ttest()] for independent sample and paired t-tests.
#'
#' @export

# param niter numeric variable for the number of iterations for permutations test. Will be ignored if mult_comp="fdr"
bstr_anova <- function(main_effect="", covariates="", bstr_data, mult_comp="fdr") { #, niter=5000

  if (inherits(bstr_data, "BstrROIData")) {
    return(bstr_roi_anova(main_effect = main_effect, covariates = covariates, bstr_data = bstr_data))
  }
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  bstr_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bstr_data = bstr_data)
  bstr_lm_null <- lm_vec(main_effect = "", covariates = covariates, bstr_data = bstr_data)
  bstr_model <- anova_vec(bstr_lm_full, bstr_lm_null, bstr_data)

  switch(mult_comp,
         # perm={
         #   pvalue_and_nulldist <- maxTperm(main_effect = main_effect, covariates = covariates, bstr_data = bstr_data, niter)
         #   bstr_model@pvalues <- pvalue_and_nulldist[[1]]
         #   bstr_model@pvalues[is.nan(bstr_model@pvalues)] <- 1
         #   bstr_model@pvalues <- bstr_model@pvalues*bstr_model@tvalues_sign
         #   bstr_model@tvalues[abs(bstr_model@pvalues) >= 0.05] <- 0
         #   nulldist <- pvalue_and_nulldist[[2]]
         #   bstr_model@pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, bstr_data, nulldist)
         #   bstr_model@tvalues_adjusted <- bstr_model@tvalues
         #   bstr_model@tvalues_adjusted[abs(bstr_model@pvalues_adjusted) >= 0.05] <- 0
         #   bstr_model@pvalues_adjusted <- bstr_model@pvalues_adjusted*bstr_model@tvalues_sign
         # },

         fdr={
           bstr_model@pvalues[is.nan(bstr_model@pvalues)] <- 1
           bstr_model@pvalues <- bstr_model@pvalues*bstr_model@tvalues_sign
           bstr_model@tvalues[abs(bstr_model@pvalues) >= 0.05] <- 0
           bstr_model@pvalues_adjusted <- bstr_p_adjust(bstr_model@pvalues, mult_comp)*bstr_model@tvalues_sign
           bstr_model@tvalues_adjusted <- bstr_model@tvalues
           bstr_model@tvalues_adjusted[abs(bstr_model@pvalues_adjusted) >= 0.05] <- 0
         }
 )

 message('Done.')
 return(bstr_model)


}

#' A vectorized version of analysis of variance (ANOVA).
#'
#' This function compares results of model fitting after [lm_vec()].
#' It accepts a full model and and a reduced model and compares them using an F-test.
#' For most scenarios, the user does not need to call this function directly. This function
#' will be called internally from [bstr_anova()]
#'
#' @param bstr_lm_full An object of type [BstrModel()] returned from [lm_vec()].
#' This is a full model including both the main effect and covariates.
#' @param bstr_lm_null An object of type [BstrModel()] returned from [lm_vec()].
#' This is a null model including only the covariates.
#' @param  bstr_data Object of type [BstrData()]
#'
#' @seealso [bstr_anova()] for most commonly used function for ANOVA, [lm_vec()] for vectorized linear regression, [ttest_vec()] for
#' vectorized independent sample and paired t-tests.
#'
#' @export
anova_vec <- function(bstr_lm_full, bstr_lm_null, bstr_data) {

  N <- nrow(bstr_data@data_array)
  Fstat <- (bstr_lm_null@rss - bstr_lm_full@rss)/bstr_lm_full@rss * (N - bstr_lm_full@Npfull - 1)/(bstr_lm_full@Npfull - bstr_lm_null@Npnull)  # F statistic

  model_unique_idx <- which(bstr_lm_null@unique %in% bstr_lm_null@fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept

  se_full_unique <- sqrt(diag(solve(t(bstr_lm_full@X_design_full) %*% bstr_lm_full@X_design_full)))[model_unique_idx] *
    sqrt(bstr_lm_full@rss / (N - bstr_lm_full@Npfull - 1))

  tvalues <- bstr_lm_full@beta_coeff[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  pvalues <- 1 - pf(Fstat, bstr_lm_full@Npfull - bstr_lm_null@Npnull, N - bstr_lm_full@Npfull - 1)
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  tvalues_sign <- sign_tvalues(tvalues)

  bstr_model <- new("BstrModel", model_type="bstr_anova", main_effect = bstr_lm_full@main_effect, covariates = bstr_lm_full@covariates,
                   demographics = bstr_data@demographics, mspec_file="")
  bstr_model@pvalues <- pvalues
  bstr_model@tvalues <- tvalues
  bstr_model@tvalues_sign <- tvalues_sign
  bstr_model@se <- se_full_unique
  bstr_model@Fstat <- Fstat
  return(bstr_model)
}

#' Linear regression for brain imaging data.
#'
#' This function accepts a `main_effect` and a set of covariates (using the R formula notation) and performs
#' a linear regression including `main_effect + covariates`.
#'
#' Slightly different from the standard R `lm` function, `bstr_lm` currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the `main_effect`
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a `+` operator similar to an R formula.
#' @param  bstr_data Object of type [BstrData()]
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See [bstr_p_adjust()] for valid values.
#'
#' @export

# param  niter numeric variable for the number of iterations for permutations test. Will be ignored if mult_comp="fdr"
bstr_lm <- function(main_effect="", covariates="", bstr_data, mult_comp = "fdr") { #, niter=5000

  if (inherits(bstr_data, "BstrROIData")) {
    return(bstr_roi_anova(main_effect = main_effect, covariates = covariates, bstr_data = bstr_data))
  }
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  bstr_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bstr_data = bstr_data)
  bstr_lm_null <- lm_vec(main_effect = "", covariates = covariates, bstr_data = bstr_data)
  bstr_model <- anova_vec(bstr_lm_full, bstr_lm_null, bstr_data)

  switch(mult_comp,
         # perm={
         #   pvalue_and_nulldist <- maxTperm(main_effect = main_effect, covariates = covariates, bstr_data = bstr_data, niter)
         #   bstr_model@pvalues <- pvalue_and_nulldist[[1]]
         #   bstr_model@pvalues[is.nan(bstr_model@pvalues)] <- 1
         #   bstr_model@pvalues <- bstr_model@pvalues*bstr_model@tvalues_sign
         #   #bstr_model@tvalues[abs(bstr_model@pvalues) >= 0.05] <- 0
         #   nulldist <- pvalue_and_nulldist[[2]]
         #   bstr_model@pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, bstr_data, nulldist)
         #   bstr_model@tvalues_adjusted <- bstr_model@tvalues
         #   bstr_model@tvalues_adjusted[abs(bstr_model@pvalues_adjusted) >= 0.05] <- 0
         #   bstr_model@pvalues_adjusted <- bstr_model@pvalues_adjusted*bstr_model@tvalues_sign
         # },

         fdr={
           bstr_model@pvalues[is.nan(bstr_model@pvalues)] <- 1
           bstr_model@pvalues <- bstr_model@pvalues*bstr_model@tvalues_sign
           bstr_model@tvalues[abs(bstr_model@pvalues) >= 0.05] <- 0
           bstr_model@pvalues_adjusted <- bstr_p_adjust(bstr_model@pvalues, mult_comp)*bstr_model@tvalues_sign
           bstr_model@tvalues_adjusted <- bstr_model@tvalues
           bstr_model@tvalues_adjusted[abs(bstr_model@pvalues_adjusted) >= 0.05] <- 0
         }
 )

 message('Done.')
 return(bstr_model)

  # # Check the model type and call the appropriate method
  # Xtemp <- solve(t(bstr_model@X_design_full) %*% bstr_model@X_design_full) %*% t(bstr_model@X_design_full) # Pre Hat matrix
  # beta_full <- Xtemp %*% bstr_data@data_array  # beta coefficients
  # y_full <- bstr_model@X_design_full %*% beta_full  # Predicted response
  # RSS_full <- colSums((bstr_data@data_array - y_full)^2)
  #
  # Xtemp <- solve(t(bstr_model@X_design_null) %*% bstr_model@X_design_null) %*% t(bstr_model@X_design_null) # Pre Hat matrix
  # beta_null <- Xtemp %*% bstr_data@data_array  # beta coefficients
  # y_null <- bstr_model@X_design_null %*% beta_null  # Predicted response
  # RSS_null <- colSums((bstr_data@data_array - y_null)^2)
  #
  # N <- nrow(bstr_data@data_array)
  # Fstat <- (RSS_null - RSS_full)/RSS_full * (N - bstr_model@Npfull - 1)/(bstr_model@Npfull - bstr_model@Npnull)  # F statistic
  # model_unique_idx <- which(bstr_model@unique %in% bstr_model@fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept
  #
  # se_full_unique <- sqrt(diag(solve(t(bstr_model@X_design_full) %*% bstr_model@X_design_full)))[model_unique_idx] *
  #   sqrt(RSS_full / (N - bstr_model@Npfull - 1))
  #
  # tvalue_sign <- (beta_full[model_unique_idx, ] + .Machine$double.eps)/(abs(beta_full[model_unique_idx, ]) + .Machine$double.eps)
  #
  # pvalues <- 1 - pf(Fstat, bstr_model@Npfull - bstr_model@Npnull, N - bstr_model@Npfull - 1)
  #
  # pvalues[is.nan(pvalues)] <- 1
  #
  # pvalues <- pvalues*tvalue_sign
  # tvalues <- beta_full[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  # bstr_model@pvalues <- pvalues
  # bstr_model@tvalues <- tvalues
  # bstr_model@tvalues[abs(pvalues) >= 0.05] <- 0
  # bstr_model@pvalues_adjusted <- p.adjust(bstr_model@pvalues, 'BH')

}

#' Vectorized linear regression for brain imaging phenotypes.
#'
#' This function accepts a `main_effect` and a set of covariates (using the R formula notation) and performs
#' a linear regression including `main_effect + covariates`.
#'
#' Slightly different from the standard R `lm` function, `lm_vec` currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the `main_effect`
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a `+` operator similar to an R formula.
#' @param  bstr_data Object of type [BstrData()]
#'
#' @export
lm_vec <- function(main_effect = "", covariates = "", bstr_data) {

  bstr_model <- new("BstrModel", model_type="bstr_lm", main_effect = main_effect, covariates = covariates,
                   demographics = bstr_data@demographics, mspec_file="")

  # lm_formula <- formula(sprintf('~ %s', paste(main_effect, '+', covariates)))
  lm_formula <- bstr_model@lm_formula

  # Fit model
  N <- dim(bstr_data@data_array)[1]
  Np <- length(unlist(strsplit(as.character(lm_formula)[2], '\\+')))

  X <- model.matrix(lm_formula, data = bstr_data@demographics)
  X_hat <- solve(t(X) %*% X) %*% t(X) # pre hat matrix
  beta_coeff <- X_hat %*% bstr_data@data_array  # beta coefficients
  Y <- X %*% beta_coeff  # predicted response
  rss <- colSums((bstr_data@data_array - Y)^2) # residual sum of squares

  if (main_effect == "")
    main_effect = "(Intercept)" # If main_effect is empty, return the parameters of the Intercept

  se <- sqrt(diag(solve(t(X) %*% X)))[[main_effect]] * sqrt(rss / (N-Np-1)) # standard error
  tvalues <- as.numeric(beta_coeff[main_effect, ]/(se + .Machine$double.eps)) # tvalues
  pvalues <- 2*pt(abs(tvalues), N-Np-1, lower.tail = FALSE) # pvalues
  residuals <- bstr_data@data_array - Y
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  bstr_model@pvalues <- pvalues
  bstr_model@tvalues <- tvalues
  bstr_model@beta_coeff <- beta_coeff
  bstr_model@rss <- rss
  bstr_model@residuals <- residuals
  return(bstr_model)
}


#' This function performs the ANOVA analysis for ROIs
#'
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a `+` operator similar to an R formula.
#' @param  bstr_data Object of type [BstrData()]
#'
#' @export

bstr_roi_anova <- function(main_effect="", covariates="", bstr_data=bstr_data) {
  # Check the model type and call the appropriate method
  bstr_model <- new("BstrModel", model_type="bstr_lm", main_effect = main_effect, covariates = covariates,
                   demographics = bstr_data@demographics, mspec_file="")
  message('Running the statistical model. This may take a while...', appendLF = FALSE)

  selected_col <- rep(NA, length(bstr_data@roiids))
  for (i in 1:length(bstr_data@roiids)){
    selected_col[i] <- paste0(as.character(get_roi_tag(read_label_desc(),bstr_data@roiids[i])), "(",bstr_data@roiids[i],")")
    bstr_data@demographics[,selected_col[i]]
  }

  cmd1 <- list()
  cmd2 <- list()
  cmd3 <- list()
  stats_commands <- list()

  for (i in 1:length(bstr_data@roiids)){
    cmd1[[i]] <- sprintf("lm_full_%s <- lm(%s, data = bstr_data@demographics)",as.character(bstr_data@roiids[i]),
                         paste('`',as.character(selected_col[i]),'`', ' ~ ', bstr_model@fullmodel, sep = ''))
    cmd2[[i]] <- sprintf("lm_null_%s <- lm(%s, data = bstr_data@demographics)",as.character(bstr_data@roiids[i]),
                         paste('`',as.character(selected_col[i]),'`', ' ~ ', bstr_model@nullmodel, sep = ''))
    cmd3[[i]] <- sprintf("pander::pander(anova(lm_full_%s, lm_null_%s))",as.character(bstr_data@roiids[i]),as.character(bstr_data@roiids[i]))
    stats_commands[[i]] <- c(cmd1[[i]], cmd2[[i]], cmd3[[i]])
  }

  for (i in 1:length(bstr_data@roiids)){
    for (cmd in stats_commands[[i]]) {
      eval(parse(text = cmd))
    }
  }
  bstr_model@stats_commands <- stats_commands

  bstr_model@load_data_command <- sprintf("bstr_model <- bstr_anova( main_effect = '%s', covariates = '%s', bstr_data = bstr_data) ",
                                         main_effect, covariates)

  return(bstr_model)

}


# bstr_roi_lm <- function(main_effect="", covariates="", bstr_data=bstr_data) {
#   # Check the model type and call the appropriate method
#   bstr_model <- new("BstrModel", model_type="bstr_lm", main_effect = main_effect, covariates = covariates,
#                    demographics = bstr_data@demographics, mspec_file="")
#   message('Running the statistical model. This may take a while...', appendLF = FALSE)
#
#   bstr_data@demographics[,as.character(bstr_data@roiids)]
#
#   cmd1 <- sprintf("lm_full_%s <- lm(%s, data = bstr_data@demographics)",as.character(bstr_data@roiids),
#                   paste('`',as.character(bstr_data@roiids),'`', ' ~ ', bstr_model@fullmodel, sep = ''))
#   cmd2 <- sprintf("lm_null_%s <- lm(%s, data = bstr_data@demographics)",as.character(bstr_data@roiids),
#                   paste('`',as.character(bstr_data@roiids),'`', ' ~ ', bstr_model@nullmodel, sep = ''))
#   cmd3 <- sprintf("pander::pander(anova(lm_full_%s, lm_null_%s))",as.character(bstr_data@roiids),as.character(bstr_data@roiids))
#
#   stats_commands <- c(cmd1, cmd2, cmd3)
#
#   for (cmd in stats_commands) {
#     eval(parse(text = cmd))
#   }
#   bstr_model@stats_commands <- stats_commands
#
#   return(bstr_model)
#
# }

#' Test for Correlation between a variable `corr_var` and a brain imaging phenotype.
#'
#' Test for correlation between a brain imaging phenotype (cortical thickness, determinant
#' of the jacobian matrix) and `corr_var` using the Pearson's product moment correlation
#' coefficient. The brain imaging phenotype is automatically selected from the type of `bstr_data`.
#' @param corr_var Character variable name. This should be present in the demographics csv file associated
#' with `bstr_data`.
#' @param  bstr_data Object of type [BstrData()].
#' @param  group_var Character variable for groups in the data. Currently this argument is only used in plotting for ROI analysis.
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See [bstr_p_adjust()] for valid values.
#' @details
#' `bstr_data` can be of the type "sba", "tbm", or "roi".
#'
#' @export
bstr_corr <- function(corr_var, bstr_data, group_var = "", mult_comp="fdr") {

  message('Running correlations...', appendLF = FALSE)
  bstr_model <- new("BstrModel", model_type="bstr_corr", corr_var = corr_var, group_var = group_var,
                   demographics = bstr_data@demographics, mspec_file="")

  corr_result <- corr_vec(bstr_data@data_array, bstr_data@demographics[[corr_var]])
  corr_coeff <- corr_result$corr_coeff
  bstr_model@tvalues <- corr_result$tvalues
  bstr_model@pvalues <- corr_result$pvalues
  tvalues_sign <- sign_tvalues(bstr_model@tvalues)
  bstr_model@tvalues_sign <- tvalues_sign

  bstr_model@pvalues <- sign(corr_coeff)*bstr_model@pvalues
  bstr_model@pvalues[abs(bstr_model@pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  bstr_model@pvalues[is.na(bstr_model@pvalues)] <- 1 # Set the p-values with the NA correlations to 1
  bstr_model@corr_values <- corr_coeff
  bstr_model@pvalues_adjusted <- bstr_p_adjust(bstr_model@pvalues, mult_comp)*bstr_model@tvalues_sign
  # bstr_model@pvalues_adjusted <- p.adjust(bstr_model@pvalues, 'BH')
  bstr_model@corr_values[abs(bstr_model@pvalues) >= 0.05] <- 0
  bstr_model@corr_values_masked_adjusted <- bstr_model@corr_values
  bstr_model@corr_values_masked_adjusted[abs(bstr_model@pvalues_adjusted) >= 0.05] <- 0
  bstr_model@tvalues_adjusted <- bstr_model@tvalues
  bstr_model@tvalues_adjusted[abs(bstr_model@pvalues_adjusted) >= 0.05] <- 0

  bstr_model@load_data_command <- sprintf("bstr_model <- bstr_corr(corr_var = '%s', bstr_data = bstr_data, mult_comp= '%s')",
                                         corr_var, mult_comp)

  if (bstr_data@analysis_type == "roi") {
    rois_string <- paste(as.character(bstr_data@roiids), collapse = ', ')
    corr_values_string <- sprintf("%.2f", bstr_model@corr_values)

    paste(as.character(bstr_model@corr_values), collapse = ', ')
    message(sprintf("Correlation(s) of %s with %s for ROI(s) %s are %s: ", bstr_data@roimeas, corr_var, rois_string, corr_values_string), appendLF = TRUE)
  }
  message('Done.')
  return(bstr_model)

}

#' Vectorized correlation between a variable and brain imaging data.
#'
#' For most scenarios, the user does not need to call this function directly.
#' Instead call [bstr_corr()] which calls this function internally.
#' @param X matrix of dimensions (\eqn{N x T}), where \eqn{N} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param Y vector of length \eqn{N}.
#' @export
corr_vec <- function(X, Y) {

  N <- length(Y)
  X_dev <- sweep(X, 2, colMeans(X))
  Y_dev <- Y - mean(Y)
  corr_coeff <- as.numeric((Y_dev %*% X_dev)/sqrt(colSums(X_dev^2)*sum(Y_dev^2)))
  tvalues <- corr_coeff * sqrt((N-2)/(1-corr_coeff^2 + .Machine$double.eps))
  pvalues <- 2*pt(abs(tvalues), N-2, lower.tail = FALSE)
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  return(list("tvalues"=tvalues, "pvalues"=pvalues, "corr_coeff"=corr_coeff))
}

#' T-test for for brain imaging phenotypes
#'
#' Perform independent sample and paired sample t-tests for differences between means of brain
#' imaging phenotypes for a categorical variable.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with `bstr_data`.
#' @param  bstr_data Object of type [BstrData()].
#' @param  paired logical; is TRUE if `group_var` contains matching (dependent) samples. The default value is `FALSE`.
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See [bstr_p_adjust()] for valid values.
#' @details
#' The degrees of freedom are calculated using the Welch–Satterthwaite approximation by default.
#' `bstr_data` can be of the type "sba", "tbm", or "roi".
#'
#' @export
bstr_ttest <- function(group_var, bstr_data, paired = FALSE, mult_comp="fdr") {

  if (paired == FALSE)
    bstr_model <- new("BstrModel", model_type="unpairedttest", group_var = group_var,
                     demographics = bstr_data@demographics, mspec_file="")
  else
    bstr_model <- new("BstrModel", model_type="pairedttest", group_var = group_var,
                     demographics = bstr_data@demographics, mspec_file="")

  group1 <- levels(as.factor(bstr_data@demographics[[group_var]]))[1]
  group2 <- levels(as.factor(bstr_data@demographics[[group_var]]))[2]
  idx_group1 <- which(bstr_data@demographics[[group_var]] == group1)
  idx_group2 <- which(bstr_data@demographics[[group_var]] == group2)

  message(sprintf("The group variable is %s with two levels (%s, %s).", group_var, group1, group2), appendLF = TRUE)
  message(sprintf("Testing for the difference between the %s groups (%s - %s).", group_var, group1, group2), appendLF = TRUE)
  message('Running t-tests...', appendLF = FALSE)

  test_result <- ttest_vec(bstr_data@data_array[idx_group1,], bstr_data@data_array[idx_group2,], group_var, paired)

  bstr_model@pvalues <- test_result$pvalues
  bstr_model@tvalues <- test_result$tvalues
  bstr_model@tvalues_sign <- sign_tvalues(bstr_model@tvalues)

  bstr_model@pvalues[is.na(bstr_model@pvalues)] <- 1
  bstr_model@pvalues <- bstr_model@pvalues*bstr_model@tvalues_sign
  bstr_model@tvalues[abs(bstr_model@pvalues) >= 0.05] <- 0
  bstr_model@pvalues_adjusted <- bstr_p_adjust(bstr_model@pvalues, mult_comp)*bstr_model@tvalues_sign
  # bstr_model@pvalues_adjusted <- p.adjust(abs(bstr_model@pvalues), 'BH')
  bstr_model@tvalues_adjusted <- bstr_model@tvalues
  bstr_model@tvalues_adjusted[abs(bstr_model@pvalues_adjusted) >= 0.05] <- 0

  bstr_model@load_data_command <- sprintf("bstr_model <- bstr_ttest(group_var = '%s', bstr_data = bstr_data, paired = %s, mult_comp= '%s')",
                                         group_var, paired,mult_comp)

  message('Done.')
  return(bstr_model)
}

#' Vectorized t-test for for brain imaging phenotypes
#'
#' Perform independent sample and paired sample t-tests between two numerfor differences between means of brain
#' imaging phenotypes for a categorical variable.
#' For most scenarios, the user does not need to call this function directly. This function
#' will be called internally from [bstr_ttest()]
#' @param X1 matrix of dimensions (\eqn{N1 x T}), where \eqn{N1} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param  X2 matrix of dimensions (\eqn{N2 x T}), where \eqn{N2} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with `bstr_data`.
#' @param  paired logical; is TRUE if `group_var` contains matching (dependent) samples. The default value is `FALSE`.
#' @details
#' For an independent samples t-test \eqn{N1} not equal to \eqn{N2}.
#' For a dependent (paired) samples t-test, \eqn{N1 = N2}.
#' The degrees of freedom are calculated using the Welch–Satterthwaite approximation by default.
#' `bstr_data` can be of the type "sba", "tbm", or "roi".
#'
#' @export
ttest_vec <- function(X1, X2, group_var, paired=FALSE) {

  n1 <- dim(X1)[1]
  n2 <- dim(X2)[1]

  if (paired == TRUE) {
    D = X1 - X2
    D_mean <- colMeans(D)
    D_dev <- sweep(D, 2, D_mean)
    s1 <- sqrt(colSums(D_dev^2)/(n1-1))
    tvalues <- D_mean/(s1/sqrt(n1))
    pvalues <- 2*pt(abs(tvalues),n1-1, lower.tail = FALSE)
    pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
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
    tvalues <- (X1_mean - X2_mean)/(se_diff + .Machine$double.eps)
    # Calculate the degrees of freedom using the Welch–Satterthwaite approximation
    deg <- (s1_sq/n1 + s2_sq/n2)^2/(s1_sq^2/(n1^2*(n1-1)) + s2_sq^2/(n2^2*(n2-1)))
    pvalues <- 2*pt(abs(tvalues), deg, lower.tail = FALSE)
    pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  }
  return(list("tvalues"=tvalues, "pvalues"=pvalues))
}

#' Adjust p-values for multiple comparisons testing
#'
#' Perform multiple comparisons correction for mass univariate tests.
#' @param pvalues numeric vector of p values.
#' @param method character string specifying the method for correction. The default method is 'fdr'. Supported methods are
#' all types given in [p.adjust.methods()]
#'
#' @export
bstr_p_adjust <- function(pvalues, method='fdr') {

  valid_methods <- c(p.adjust.methods) #, "perm"
  if ( !(method  %in% valid_methods) ) {
    warning(sprintf("%s is not a valid multiple comparisons method. Using no correction.", method), call. = FALSE)
    return(pvalues)
  }

  if (method %in% p.adjust.methods && method != "perm") {
    return(p.adjust(abs(pvalues), method))
  }

  # if (method == "perm") {
  #   tvalues_null <- bstr_perm()
  # }


}

# Calculate p-values for vanilla permutation test and determine null distribution for max-t permutation test.
#
# Perform permutation tests using Freedman-Lane method.
# @param main_effect Character string containing an independent variable whose effect you want to measure.
# It could be disease status, age, gender etc. This should strictly be a single variable. This can be
# either a categorical or a continuous variable.
# @param covariates Character string containing a set of other predictors (variables) in the model. If more than
# one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
# @param bstr_data Object of type \code{\link{BstrData}}
# @param num_of_perm Number of iterations/shuffles for permutation test.
# @details
# The permutation test handles the exchangeability assumption with Freedman-Lane Method. This function also
# utilizes multiprocessing to reduce computing time.
# \code{bstr_data} can be of the type "sba", "tbm", or "roi".
#
# @export
# maxTperm <- function(main_effect = "", covariates = "", bstr_data, num_of_perm){
#
#   N <- dim(bstr_data@data_array)[1]
#
#   bstr_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bstr_data = bstr_data)
#   bstr_lm_null <- lm_vec(main_effect = "", covariates = covariates, bstr_data = bstr_data)
#
#   # (1) compute full model t-statistic
#   T0 <- bstr_lm_full@tvalues
#   maxT0 <- T0[ which.max( abs(T0) ) ]
#
#   # (2) compute estimated gamma_hat and estimated residuals from reduced model
#   # gamma_hat <- bstr_lm_null@beta_coeff
#   # residuals_null <- bstr_lm_null@residuals
#
#   # (3) compute a set of permuted data Y
#   bstr_data_cp <- bstr_data
#
#   t_bin_int <- rep(0, dim(bstr_data@data_array)[2])
#   t_max_per_perm <- rep(0, num_of_perm-1)
#   # rewrote into for loop, single core
#   for (j in 1:(num_of_perm-1)){
#     set.seed(j)
#     pmatrix <- as(sample(N), "pMatrix")
#     Y_j <- (pmatrix %*% bstr_lm_null@residuals) + (bstr_lm_null@X_design_null %*% bstr_lm_null@beta_coeff)
#
#     # (4) regress permuted data Y_j against the full model
#     bstr_data_cp@data_array <- Y_j
#     bstr_lm_full_perm <- lm_vec(main_effect = main_effect, covariates = covariates, bstr_data = bstr_data_cp )
#
#     ## binarize vector after comparing permuted T and observed T, where T0 is the vector containing
#     # all t values from the observed model
#     idx <- which( abs(T0) <= abs(bstr_lm_full_perm@tvalues) )
#
#     t_bin_int[idx] <- t_bin_int[idx] + 1
#
#     # t_bin <- as.bit(rep(FALSE, dim(bstr_data_cp@data_array)[2]))
#     # t_bin[idx] <- TRUE
#
#     t_max_per_perm[j] <- bstr_lm_full_perm@tvalues[ which.max( abs(bstr_lm_full_perm@tvalues) ) ] # max tvalue from perm
#
#   }
#
#   pvalues <- rep(0, dim(bstr_data@data_array)[2])
#   for (n in 1:dim(bstr_data@data_array)[2]){
#     count <- t_bin_int[n]
#     pvalues[n] <- as.double((count+1)/num_of_perm)
#   }
#
#   return(list(pvalues, t_max_per_perm))
#
# }


# Adjust p-values for multiple comparisons testing
#
# Perform multiple comparisons correction for mass univariate tests using the max-t method.
# @param main_effect Character string containing an independent variable whose effect you want to measure.
# It could be disease status, age, gender etc. This should strictly be a single variable. This can be
# either a categorical or a continuous variable.
# @param covariates Character string containing a set of other predictors (variables) in the model. If more than
# one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
# @param bstr_data Object of type \code{\link{BstrData}}
# @param tvalues_null Null distribution output from \code{\link{maxTperm}}. Statistics drawn from each
# shuffle are the maximum t-statistic within ROI.
#
# @export
# perm_p_adjust <- function(main_effect = "", covariates = "", bstr_data, tvalues_null){
#   bstr_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bstr_data = bstr_data)
#
#   tvalues <- bstr_lm_full@tvalues
#   pvalues_adj <- rep(0, length(tvalues))
#
#   for (i in 1:length(tvalues)){
#     p <- (sum(abs(tvalues_null) >= abs(tvalues[i]))+1) / (length(tvalues_null)+1)
#     pvalues_adj[i] <- p
#   }
#   return(pvalues_adj)
# }

#' Linear mixed-effects model for brain imaging data.
#'
#' Linear regression for brain imaging data.
#'
#' This function accepts a `main_effect` and a set of covariates (using the R formula notation) and performs
#' a linear regression including `main_effect + covariates`.
#'
#' Slightly different from the standard R `lmer` function, `bstr_lmer` currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the `main_effect`
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with `bstr_data`.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a `+` operator similar to an R formula.
#' @param  bstr_data Object of type [BstrData()]
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See [bstr_p_adjust()] for valid values.
#'
#' @export
bstr_lmer <- function(group_var, main_effect="", covariates="", bstr_data, mult_comp = "fdr") {

  if (inherits(bstr_data, "BstrROIData")) {
    return(bstr_roi_lmer_anova(group_var, main_effect = main_effect, covariates = covariates, bstr_data = bstr_data))
  }
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  bstr_model <- lmer_vec(group_var, main_effect = main_effect, covariates = covariates, bstr_data = bstr_data)
  # switch(mult_comp,
  #        perm={
  #          cl <- parallel::makeCluster(parallel::detectCores())
  #          registerDoParallel(cl)
  #          options(warn=-1)
  #          pvalue_and_nulldist <- maxTperm(main_effect = main_effect, covariates = covariates, bstr_data = bstr_data, niter)
  #          bstr_model@pvalues <- pvalue_and_nulldist[[1]]
  #          bstr_model@pvalues[is.nan(bstr_model@pvalues)] <- 1
  #          bstr_model@pvalues <- bstr_model@pvalues*bstr_model@tvalues_sign
  #         #bstr_model@tvalues[abs(bstr_model@pvalues) >= 0.05] <- 0
  #          nulldist <- pvalue_and_nulldist[[2]]
  #          bstr_model@pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, bstr_data, nulldist)
  #          bstr_model@tvalues_adjusted <- bstr_model@tvalues
  #          bstr_model@tvalues_adjusted[abs(bstr_model@pvalues_adjusted) >= 0.05] <- 0
  #          stopCluster(cl)
  #          options(warn=0)
  #        },
  #        fdr={
  #          bstr_model@pvalues[is.nan(bstr_model@pvalues)] <- 1
  #          bstr_model@pvalues <- bstr_model@pvalues*bstr_model@tvalues_sign
  #          #bstr_model@tvalues[abs(bstr_model@pvalues) >= 0.05] <- 0
  #          bstr_model@pvalues_adjusted <- bstr_p_adjust(bstr_model@pvalues, mult_comp)
  #          bstr_model@tvalues_adjusted <- bstr_model@tvalues
  #          bstr_model@tvalues_adjusted[abs(bstr_model@pvalues_adjusted) >= 0.05] <- 0
  #        }
  # )

  message('Done.')
  return(bstr_model)
}

#' Vectorized linear mixed-effects regression for brain imaging phenotypes.
#'
#' This function accepts a `main_effect` and a set of covariates (using the R formula notation) and performs
#' a linear regression including `main_effect + covariates`.
#'
#' Slightly different from the standard R `lmer` function, `lmer_vec` currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the `main_effect`
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with `bstr_data`.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a `+` operator similar to an R formula.
#' @param  bstr_data Object of type [BstrData()]
#'
#' @export
lmer_vec <- function(group_var, main_effect = "", covariates = "", bstr_data) {

  bstr_model <- new("BstrModel", model_type="bstr_lmer", main_effect = main_effect, covariates = covariates,
                   group_var = group_var, demographics = bstr_data@demographics, mspec_file="")

  # Mass Univariate Linear Mixed Effects Analysis (still in progress)
  #lm_formula <- formula(sprintf('~ %s', paste(main_effect, '+', covariates)))
  #lmer_formula <- bstr_model@lm_formula

  # Fit model
  N <- dim(bstr_data@data_array)[1]
  # Np <- length(unlist(strsplit(as.character(covariates), '\\+'))) + 1
  #
  # X <- model.matrix(lm_formula, data = bstr_data@demographics)
  # X_hat <- solve(t(X) %*% X) %*% t(X) # pre hat matrix
  # beta_coeff <- X_hat %*% bstr_data@data_array  # beta coefficients
  # Y <- X %*% beta_coeff  # predicted response
  # rss <- colSums((bstr_data@data_array - Y)^2) # residual sum of squares
  #
  # if (main_effect == "")
  #   main_effect = "(Intercept)" # If main_effect is empty, return the parameters of the Intercept
  #
  # se <- sqrt(diag(solve(t(X) %*% X)))[[main_effect]] * sqrt(rss / (N-Np-1)) # standard error
  # tvalues <- as.numeric(beta_coeff[main_effect, ]/(se + .Machine$double.eps)) # tvalue
  # pvalues <- 2*pt(abs(tvalues), N-Np-1, lower.tail = FALSE) # pvalu
  # residuals <- bstr_data@data_array - Y
  # pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  # bstr_model@pvalues <- pvalues
  # bstr_model@tvalues <- tvalues
  # bstr_model@beta_coeff <- beta_coeff
  # bstr_model@rss <- rss
  # bstr_model@residuals <- residuals

  return(bstr_model)
}

#' This function performs the ANOVA Linear Mixed Effects analysis for ROIs
#'
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with `bstr_data`.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a `+` operator similar to an R formula.
#' @param  bstr_data Object of type [BstrData()]
#'
#' @export

bstr_roi_lmer_anova <- function(group_var, main_effect="", covariates="", bstr_data){

  bstr_model <- new("BstrModel", model_type="bstr_lmer", main_effect = main_effect, covariates = covariates,
                   group_var = group_var, demographics = bstr_data@demographics, mspec_file="")
  message('Running the statistical model. This may take a while...', appendLF = FALSE)

  # For univariate lmer
  selected_col <- rep(NA, length(bstr_data@roiids))
  for (i in 1:length(bstr_data@roiids)){
    selected_col[i] <- paste0(as.character(get_roi_tag(read_label_desc(),bstr_data@roiids[i])), "(",bstr_data@roiids[i],")")
  }

  cmd1 <- list()
  cmd2 <- list()
  cmd3 <- list()
  stats_commands <- list()

  for (i in 1:length(bstr_data@roiids)){
    bstr_lmer_full_formula <- paste0("`",as.character(get_roi_tag(read_label_desc(),bstr_data@roiids[i])), "(",bstr_data@roiids[i],")`"," ~ ",bstr_model@main_effect," + (",bstr_model@main_effect,"|",bstr_model@group_var,") + ",bstr_model@covariates)
    bstr_lmer_null_formula <- paste0("`",as.character(get_roi_tag(read_label_desc(),bstr_data@roiids[i])), "(",bstr_data@roiids[i],")`"," ~ (1|",bstr_model@group_var,") + ",bstr_model@covariates)

    cmd1[[i]] <- sprintf("lmer_full_%s <- lme4::lmer(%s, data = bstr_data@demographics,control=lme4::lmerControl(check.nobs.vs.nRE='ignore'))",
                         as.character(bstr_data@roiids[i]),bstr_lmer_full_formula)
    cmd2[[i]] <- sprintf("lmer_null_%s <- lme4::lmer(%s, data = bstr_data@demographics,control=lme4::lmerControl(check.nobs.vs.nRE='ignore'))",
                         as.character(bstr_data@roiids[i]),bstr_lmer_null_formula)
    cmd3[[i]] <- sprintf("pander::pander(anova(lmer_full_%s, lmer_null_%s))",as.character(bstr_data@roiids[i]),as.character(bstr_data@roiids[i]))
    stats_commands[[i]] <- c(cmd1[[i]], cmd2[[i]], cmd3[[i]])
  }

  eval(parse(text = stats_commands))
  bstr_model@stats_commands <- stats_commands
  bstr_model@load_data_command <- sprintf("bstr_model <- bstr_lmer(group_var = '%s', main_effect = '%s', covariates = '%s', bstr_data = bstr_data) ",
                                         group_var, main_effect, covariates)

  message('Done.')
  return (bstr_model)
}
