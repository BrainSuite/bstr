# BrainSuite Statistics Toolbox in R (bssr)
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
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and uses an
#' F-test to compare the full model including the \code{main_effect + covariates} with the reduced (null) model
#' that only includes the \code{covariates}.
#'
#' Slightly different from the standard R anova function, \code{bss_anova} currently does not directly accept the
#' results from \code{lm_vec}. This could be accomodated in the future versions.

#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  bss_data Object of type \code{\link{BssData}}
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{bss_p_adjust}} for valid values.
#' @param  niter numeric variable for the number of iterations for permutations test. Will be ignored if mult_comp="fdr"
#' @seealso \code{\link{lm_vec}} for linear regression, \code{\link{bss_ttest}} for independent sample and paired t-tests.
#'
#' @export
bss_anova <- function(main_effect="", covariates="", bss_data, mult_comp="fdr", niter=5000) {

  if (class(bss_data) == "BssROIData") {
    return(bss_roi_anova(main_effect = main_effect, covariates = covariates, bss_data = bss_data))
  }
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  bss_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data)
  bss_lm_null <- lm_vec(main_effect = "", covariates = covariates, bss_data = bss_data)
  bss_model <- anova_vec(bss_lm_full, bss_lm_null, bss_data)

  switch(mult_comp,
         perm={

           pvalue_and_nulldist <- maxTperm(main_effect = main_effect, covariates = covariates, bss_data = bss_data, niter)
           bss_model@pvalues <- pvalue_and_nulldist[[1]]
           bss_model@pvalues[is.nan(bss_model@pvalues)] <- 1
           bss_model@pvalues <- bss_model@pvalues*bss_model@tvalues_sign
           bss_model@tvalues[abs(bss_model@pvalues) >= 0.05] <- 0
           nulldist <- pvalue_and_nulldist[[2]]
           bss_model@pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, bss_data, nulldist)
           bss_model@tvalues_adjusted <- bss_model@tvalues
           bss_model@tvalues_adjusted[abs(bss_model@pvalues_adjusted) >= 0.05] <- 0
           bss_model@pvalues_adjusted <- bss_model@pvalues_adjusted*bss_model@tvalues_sign

         },
         fdr={
           bss_model@pvalues[is.nan(bss_model@pvalues)] <- 1
           bss_model@pvalues <- bss_model@pvalues*bss_model@tvalues_sign
           bss_model@tvalues[abs(bss_model@pvalues) >= 0.05] <- 0
           bss_model@pvalues_adjusted <- bss_p_adjust(bss_model@pvalues, mult_comp)*bss_model@tvalues_sign
           bss_model@tvalues_adjusted <- bss_model@tvalues
           bss_model@tvalues_adjusted[abs(bss_model@pvalues_adjusted) >= 0.05] <- 0
         }
 )

 message('Done.')
 return(bss_model)


}

#' A vectorized version of analysis of variance (ANOVA).
#'
#' This function compares results of model fitting after \code{\link{lm_vec}}.
#' It accepts a full model and and a reduced model and compares them using an F-test.
#' For most scenarios, the user does not need to call this function directly. This function
#' will be called internally from \code{\link{bss_anova}}
#'
#' @param bss_lm_full An object of type \code{\link{BssModel}} returned from \code{\link{lm_vec}}.
#' This is a full model including both the main effect and covariates.
#' @param bss_lm_null An object of type \code{\link{BssModel}} returned from \code{\link{lm_vec}}.
#' This is a null model including only the covariates.
#' @param  bss_data Object of type \code{\link{BssData}}
#'
#' @seealso \code{\link{bss_anova}} for most commonly used function for ANOVA, \code{\link{lm_vec}} for vectorized linear regression, \code{\link{ttest_vec}} for
#' vectorized independent sample and paired t-tests.
#'
#' @export
anova_vec <- function(bss_lm_full, bss_lm_null, bss_data) {

  N <- nrow(bss_data@data_array)
  Fstat <- (bss_lm_null@rss - bss_lm_full@rss)/bss_lm_full@rss * (N - bss_lm_full@Npfull - 1)/(bss_lm_full@Npfull - bss_lm_null@Npnull)  # F statistic

  model_unique_idx <- which(bss_lm_null@unique %in% bss_lm_null@fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept

  se_full_unique <- sqrt(diag(solve(t(bss_lm_full@X_design_full) %*% bss_lm_full@X_design_full)))[model_unique_idx] *
    sqrt(bss_lm_full@rss / (N - bss_lm_full@Npfull - 1))

  tvalues <- bss_lm_full@beta_coeff[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  pvalues <- 1 - pf(Fstat, bss_lm_full@Npfull - bss_lm_null@Npnull, N - bss_lm_full@Npfull - 1)
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  tvalues_sign <- sign_tvalues(tvalues)

  bss_model <- new("BssModel", model_type="bss_anova", main_effect = bss_lm_full@main_effect, covariates = bss_lm_full@covariates,
                   demographics = bss_data@demographics, mspec_file="")
  bss_model@pvalues <- pvalues
  bss_model@tvalues <- tvalues
  bss_model@tvalues_sign <- tvalues_sign
  bss_model@se <- se_full_unique
  bss_model@Fstat <- Fstat
  return(bss_model)
}

#' Linear regression for brain imaging data.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and performs
#' a linear regression including \code{main_effect + covariates}.
#'
#' Slightly different from the standard R \code{lm} function, \code{bss_lm} currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the \code{main_effect}
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  bss_data Object of type \code{\link{BssData}}
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{bss_p_adjust}} for valid values.
#' @param  niter numeric variable for the number of iterations for permutations test. Will be ignored if mult_comp="fdr"
#'
#' @export
bss_lm <- function(main_effect="", covariates="", bss_data, mult_comp = "fdr", niter=5000) {

  if (class(bss_data) == "BssROIData") {
    return(bss_roi_anova(main_effect = main_effect, covariates = covariates, bss_data = bss_data))
  }
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  bss_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data)
  bss_lm_null <- lm_vec(main_effect = "", covariates = covariates, bss_data = bss_data)
  bss_model <- anova_vec(bss_lm_full, bss_lm_null, bss_data)

  switch(mult_comp,
         perm={

           pvalue_and_nulldist <- maxTperm(main_effect = main_effect, covariates = covariates, bss_data = bss_data, niter)
           bss_model@pvalues <- pvalue_and_nulldist[[1]]
           bss_model@pvalues[is.nan(bss_model@pvalues)] <- 1
           bss_model@pvalues <- bss_model@pvalues*bss_model@tvalues_sign
           #bss_model@tvalues[abs(bss_model@pvalues) >= 0.05] <- 0
           nulldist <- pvalue_and_nulldist[[2]]
           bss_model@pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, bss_data, nulldist)
           bss_model@tvalues_adjusted <- bss_model@tvalues
           bss_model@tvalues_adjusted[abs(bss_model@pvalues_adjusted) >= 0.05] <- 0
           bss_model@pvalues_adjusted <- bss_model@pvalues_adjusted*bss_model@tvalues_sign

         },
         fdr={
           bss_model@pvalues[is.nan(bss_model@pvalues)] <- 1
           bss_model@pvalues <- bss_model@pvalues*bss_model@tvalues_sign
           bss_model@tvalues[abs(bss_model@pvalues) >= 0.05] <- 0
           bss_model@pvalues_adjusted <- bss_p_adjust(bss_model@pvalues, mult_comp)*bss_model@tvalues_sign
           bss_model@tvalues_adjusted <- bss_model@tvalues
           bss_model@tvalues_adjusted[abs(bss_model@pvalues_adjusted) >= 0.05] <- 0
         }
 )

 message('Done.')
 return(bss_model)

  # # Check the model type and call the appropriate method
  # Xtemp <- solve(t(bss_model@X_design_full) %*% bss_model@X_design_full) %*% t(bss_model@X_design_full) # Pre Hat matrix
  # beta_full <- Xtemp %*% bss_data@data_array  # beta coefficients
  # y_full <- bss_model@X_design_full %*% beta_full  # Predicted response
  # RSS_full <- colSums((bss_data@data_array - y_full)^2)
  #
  # Xtemp <- solve(t(bss_model@X_design_null) %*% bss_model@X_design_null) %*% t(bss_model@X_design_null) # Pre Hat matrix
  # beta_null <- Xtemp %*% bss_data@data_array  # beta coefficients
  # y_null <- bss_model@X_design_null %*% beta_null  # Predicted response
  # RSS_null <- colSums((bss_data@data_array - y_null)^2)
  #
  # N <- nrow(bss_data@data_array)
  # Fstat <- (RSS_null - RSS_full)/RSS_full * (N - bss_model@Npfull - 1)/(bss_model@Npfull - bss_model@Npnull)  # F statistic
  # model_unique_idx <- which(bss_model@unique %in% bss_model@fullvars) + 1    # Add 1, because the first column in the design matrix is the intercept
  #
  # se_full_unique <- sqrt(diag(solve(t(bss_model@X_design_full) %*% bss_model@X_design_full)))[model_unique_idx] *
  #   sqrt(RSS_full / (N - bss_model@Npfull - 1))
  #
  # tvalue_sign <- (beta_full[model_unique_idx, ] + .Machine$double.eps)/(abs(beta_full[model_unique_idx, ]) + .Machine$double.eps)
  #
  # pvalues <- 1 - pf(Fstat, bss_model@Npfull - bss_model@Npnull, N - bss_model@Npfull - 1)
  #
  # pvalues[is.nan(pvalues)] <- 1
  #
  # pvalues <- pvalues*tvalue_sign
  # tvalues <- beta_full[model_unique_idx, ]/(se_full_unique + .Machine$double.eps)
  # bss_model@pvalues <- pvalues
  # bss_model@tvalues <- tvalues
  # bss_model@tvalues[abs(pvalues) >= 0.05] <- 0
  # bss_model@pvalues_adjusted <- p.adjust(bss_model@pvalues, 'BH')

}

#' Vectorized linear regression for brain imaging phenotypes.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and performs
#' a linear regression including \code{main_effect + covariates}.
#'
#' Slightly different from the standard R \code{lm} function, \code{lm_vec} currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the \code{main_effect}
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  bss_data Object of type \code{\link{BssData}}
#'
#' @export
lm_vec <- function(main_effect = "", covariates = "", bss_data) {

  bss_model <- new("BssModel", model_type="bss_lm", main_effect = main_effect, covariates = covariates,
                   demographics = bss_data@demographics, mspec_file="")

  # lm_formula <- formula(sprintf('~ %s', paste(main_effect, '+', covariates)))
  lm_formula <- bss_model@lm_formula

  # Fit model
  N <- dim(bss_data@data_array)[1]
  Np <- length(unlist(strsplit(as.character(lm_formula)[2], '\\+')))

  X <- model.matrix(lm_formula, data = bss_data@demographics)
  X_hat <- solve(t(X) %*% X) %*% t(X) # pre hat matrix
  beta_coeff <- X_hat %*% bss_data@data_array  # beta coefficients
  Y <- X %*% beta_coeff  # predicted response
  rss <- colSums((bss_data@data_array - Y)^2) # residual sum of squares

  if (main_effect == "")
    main_effect = "(Intercept)" # If main_effect is empty, return the parameters of the Intercept

  se <- sqrt(diag(solve(t(X) %*% X)))[[main_effect]] * sqrt(rss / (N-Np-1)) # standard error
  tvalues <- as.numeric(beta_coeff[main_effect, ]/(se + .Machine$double.eps)) # tvalues
  pvalues <- 2*pt(abs(tvalues), N-Np-1, lower.tail = FALSE) # pvalues
  residuals <- bss_data@data_array - Y
  pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  bss_model@pvalues <- pvalues
  bss_model@tvalues <- tvalues
  bss_model@beta_coeff <- beta_coeff
  bss_model@rss <- rss
  bss_model@residuals <- residuals
  return(bss_model)
}


#' This function performs the ANOVA analysis for ROIs
#'
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  bss_data Object of type \code{\link{BssData}}
#'
#' @export

bss_roi_anova <- function(main_effect="", covariates="", bss_data=bss_data) {
  # Check the model type and call the appropriate method
  bss_model <- new("BssModel", model_type="bss_lm", main_effect = main_effect, covariates = covariates,
                   demographics = bss_data@demographics, mspec_file="")
  message('Running the statistical model. This may take a while...', appendLF = FALSE)

  selected_col <- rep(NA, length(bss_data@roiids))
  for (i in 1:length(bss_data@roiids)){
    selected_col[i] <- paste0(as.character(get_roi_tag(read_label_desc(),bss_data@roiids[i])), "(",bss_data@roiids[i],")")
    bss_data@demographics[,selected_col[i]]
  }

  cmd1 <- list()
  cmd2 <- list()
  cmd3 <- list()
  stats_commands <- list()

  for (i in 1:length(bss_data@roiids)){
    cmd1[[i]] <- sprintf("lm_full_%s <- lm(%s, data = bss_data@demographics)",as.character(bss_data@roiids[i]),
                         paste('`',as.character(selected_col[i]),'`', ' ~ ', bss_model@fullmodel, sep = ''))
    cmd2[[i]] <- sprintf("lm_null_%s <- lm(%s, data = bss_data@demographics)",as.character(bss_data@roiids[i]),
                         paste('`',as.character(selected_col[i]),'`', ' ~ ', bss_model@nullmodel, sep = ''))
    cmd3[[i]] <- sprintf("pander::pander(anova(lm_full_%s, lm_null_%s))",as.character(bss_data@roiids[i]),as.character(bss_data@roiids[i]))
    stats_commands[[i]] <- c(cmd1[[i]], cmd2[[i]], cmd3[[i]])
  }

  for (i in 1:length(bss_data@roiids)){
    for (cmd in stats_commands[[i]]) {
      eval(parse(text = cmd))
    }
  }
  bss_model@stats_commands <- stats_commands

  bss_model@load_data_command <- sprintf("bss_model <- bss_anova( main_effect = '%s', covariates = '%s', bss_data = bss_data) ",
                                         main_effect, covariates)

  return(bss_model)

}


# bss_roi_lm <- function(main_effect="", covariates="", bss_data=bss_data) {
#   # Check the model type and call the appropriate method
#   bss_model <- new("BssModel", model_type="bss_lm", main_effect = main_effect, covariates = covariates,
#                    demographics = bss_data@demographics, mspec_file="")
#   message('Running the statistical model. This may take a while...', appendLF = FALSE)
#
#   bss_data@demographics[,as.character(bss_data@roiids)]
#
#   cmd1 <- sprintf("lm_full_%s <- lm(%s, data = bss_data@demographics)",as.character(bss_data@roiids),
#                   paste('`',as.character(bss_data@roiids),'`', ' ~ ', bss_model@fullmodel, sep = ''))
#   cmd2 <- sprintf("lm_null_%s <- lm(%s, data = bss_data@demographics)",as.character(bss_data@roiids),
#                   paste('`',as.character(bss_data@roiids),'`', ' ~ ', bss_model@nullmodel, sep = ''))
#   cmd3 <- sprintf("pander::pander(anova(lm_full_%s, lm_null_%s))",as.character(bss_data@roiids),as.character(bss_data@roiids))
#
#   stats_commands <- c(cmd1, cmd2, cmd3)
#
#   for (cmd in stats_commands) {
#     eval(parse(text = cmd))
#   }
#   bss_model@stats_commands <- stats_commands
#
#   return(bss_model)
#
# }

#' Test for Correlation between a variable \code{corr_var} and a brain imaging phenotype.
#'
#' Test for correlation between a brain imaging phenotype (cortical thickness, determinant
#' of the jacobian matrix) and \code{corr_var} using the Pearson's product moment correlation
#' coefficient. The brain imaging phenotype is automatically selected from the type of \code{bss_data}.
#' @param corr_var Character variable name. This should be present in the demographics csv file associated
#' with \code{bss_data}.
#' @param  bss_data Object of type \code{\link{BssData}}.
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{bss_p_adjust}} for valid values.
#' @details
#' \code{bss_data} can be of the type "cbm", "tbm", or "roi".
#'
#' @export
bss_corr <- function(corr_var, bss_data, mult_comp="fdr") {

  message('Running correlations...', appendLF = FALSE)
  bss_model <- new("BssModel", model_type="bss_corr", corr_var = corr_var,
                   demographics = bss_data@demographics, mspec_file="")

  corr_result <- corr_vec(bss_data@data_array, bss_data@demographics[[corr_var]])
  corr_coeff <- corr_result$corr_coeff
  bss_model@tvalues <- corr_result$tvalues
  bss_model@pvalues <- corr_result$pvalues
  tvalues_sign <- sign_tvalues(bss_model@tvalues)
  bss_model@tvalues_sign <- tvalues_sign

  bss_model@pvalues <- sign(corr_coeff)*bss_model@pvalues
  bss_model@pvalues[abs(bss_model@pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  bss_model@pvalues[is.na(bss_model@pvalues)] <- 1 # Set the p-values with the NA correlations to 1
  bss_model@corr_values <- corr_coeff
  bss_model@pvalues_adjusted <- bss_p_adjust(bss_model@pvalues, mult_comp)*bss_model@tvalues_sign
  # bss_model@pvalues_adjusted <- p.adjust(bss_model@pvalues, 'BH')
  bss_model@corr_values[abs(bss_model@pvalues) >= 0.05] <- 0
  bss_model@corr_values_masked_adjusted <- bss_model@corr_values
  bss_model@corr_values_masked_adjusted[abs(bss_model@pvalues_adjusted) >= 0.05] <- 0
  bss_model@tvalues_adjusted <- bss_model@tvalues
  bss_model@tvalues_adjusted[abs(bss_model@pvalues_adjusted) >= 0.05] <- 0

  bss_model@load_data_command <- sprintf("bss_model <- bss_corr(corr_var = '%s', bss_data = bss_data, mult_comp= '%s')",
                                         corr_var, mult_comp)

  message('Done.')
  return(bss_model)

}

#' Vectorized correlation between a variable and brain imaging data.
#'
#' For most scenarios, the user does not need to call this function directly.
#' Instead call \code{\link{bss_corr}} which calls this function internally.
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
#' with \code{bss_data}.
#' @param  bss_data Object of type \code{\link{BssData}}.
#' @param  paired logical; is TRUE if \code{group_var} contains matching (dependent) samples. The default value is \code{FALSE}.
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{bss_p_adjust}} for valid values.
#' @details
#' The degrees of freedom are calculated using the Welch–Satterthwaite approximation by default.
#' \code{bss_data} can be of the type "cbm", "tbm", or "roi".
#'
#' @export
bss_ttest <- function(group_var, bss_data, paired = FALSE, mult_comp="fdr") {

  if (paired == FALSE)
    bss_model <- new("BssModel", model_type="unpairedttest", group_var = group_var,
                     demographics = bss_data@demographics, mspec_file="")
  else
    bss_model <- new("BssModel", model_type="pairedttest", group_var = group_var,
                     demographics = bss_data@demographics, mspec_file="")

  group1 <- levels(as.factor(bss_data@demographics[[group_var]]))[1]
  group2 <- levels(as.factor(bss_data@demographics[[group_var]]))[2]
  idx_group1 <- which(bss_data@demographics[[group_var]] == group1)
  idx_group2 <- which(bss_data@demographics[[group_var]] == group2)

  message(sprintf("The group variable is %s with two levels (%s, %s).", group_var, group1, group2), appendLF = TRUE)
  message(sprintf("Testing for the difference between the %s groups (%s - %s).", group_var, group1, group2), appendLF = TRUE)
  message('Running t-tests...', appendLF = FALSE)

  test_result <- ttest_vec(bss_data@data_array[idx_group1,], bss_data@data_array[idx_group2,], group_var, paired)

  bss_model@pvalues <- test_result$pvalues
  bss_model@tvalues <- test_result$tvalues
  bss_model@tvalues_sign <- sign_tvalues(bss_model@tvalues)

  bss_model@pvalues[is.na(bss_model@pvalues)] <- 1
  bss_model@pvalues <- bss_model@pvalues*bss_model@tvalues_sign
  bss_model@tvalues[abs(bss_model@pvalues) >= 0.05] <- 0
  bss_model@pvalues_adjusted <- bss_p_adjust(bss_model@pvalues, mult_comp)*bss_model@tvalues_sign
  # bss_model@pvalues_adjusted <- p.adjust(abs(bss_model@pvalues), 'BH')
  bss_model@tvalues_adjusted <- bss_model@tvalues
  bss_model@tvalues_adjusted[abs(bss_model@pvalues_adjusted) >= 0.05] <- 0

  bss_model@load_data_command <- sprintf("bss_model <- bss_ttest(group_var = '%s', bss_data = bss_data, paired = %s, mult_comp= '%s')",
                                         group_var, paired,mult_comp)

  message('Done.')
  return(bss_model)
}

#' Vectorized t-test for for brain imaging phenotypes
#'
#' Perform independent sample and paired sample t-tests between two numerfor differences between means of brain
#' imaging phenotypes for a categorical variable.
#' For most scenarios, the user does not need to call this function directly. This function
#' will be called internally from \code{\link{bss_ttest}}
#' @param X1 matrix of dimensions (\eqn{N1 x T}), where \eqn{N1} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param  X2 matrix of dimensions (\eqn{N2 x T}), where \eqn{N2} = number of subjects and \eqn{T} = number of vertices/voxels.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with \code{bss_data}.
#' @param  paired logical; is TRUE if \code{group_var} contains matching (dependent) samples. The default value is \code{FALSE}.
#' @details
#' For an independent samples t-test \eqn{N1} not equal to \eqn{N2}.
#' For a dependent (paired) samples t-test, \eqn{N1 = N2}.
#' The degrees of freedom are calculated using the Welch–Satterthwaite approximation by default.
#' \code{bss_data} can be of the type "cbm", "tbm", or "roi".
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
#' all types given in \code{\link{p.adjust.methods}}
#'
#' @export
bss_p_adjust <- function(pvalues, method='fdr') {

  valid_methods <- c(p.adjust.methods, "perm")
  if ( !(method  %in% valid_methods) ) {
    warning(sprintf("%s is not a valid multiple comparisons method. Using no correction.", method), call. = FALSE)
    return(pvalues)
  }

  if (method %in% p.adjust.methods && method != "perm") {
    return(p.adjust(abs(pvalues), method))
  }

  # if (method == "perm") {
  #   tvalues_null <- bss_perm()
  # }


}

#' Calculate p-values for vanilla permutation test and determine null distribution for max-t permutation test.
#'
#' Perform permutation tests using Freedman-Lane method.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param bss_data Object of type \code{\link{BssData}}
#' @param num_of_perm Number of iterations/shuffles for permutation test.
#' @details
#' The permutation test handles the exchangeability assumption with Freedman-Lane Method. This function also
#' utilizes multiprocessing to reduce computing time.
#' \code{bss_data} can be of the type "cbm", "tbm", or "roi".
#'
#' @export
maxTperm <- function(main_effect = "", covariates = "", bss_data, num_of_perm){

  N <- dim(bss_data@data_array)[1]

  bss_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data)
  bss_lm_null <- lm_vec(main_effect = "", covariates = covariates, bss_data = bss_data)

  # (1) compute full model t-statistic
  T0 <- bss_lm_full@tvalues
  maxT0 <- T0[ which.max( abs(T0) ) ]

  # (2) compute estimated gamma_hat and estimated residuals from reduced model
  # gamma_hat <- bss_lm_null@beta_coeff
  # residuals_null <- bss_lm_null@residuals

  # (3) compute a set of permuted data Y
  bss_data_cp <- bss_data

  t_bin_int <- rep(0, dim(bss_data@data_array)[2])
  t_max_per_perm <- rep(0, num_of_perm-1)
  # rewrote into for loop, single core
  for (j in 1:(num_of_perm-1)){
    set.seed(j)
    pmatrix <- as(sample(N), "pMatrix")
    Y_j <- (pmatrix %*% bss_lm_null@residuals) + (bss_lm_null@X_design_null %*% bss_lm_null@beta_coeff)

    # (4) regress permuted data Y_j against the full model
    bss_data_cp@data_array <- Y_j
    bss_lm_full_perm <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data_cp )

    ## binarize vector after comparing permuted T and observed T, where T0 is the vector containing
    # all t values from the observed model
    idx <- which( abs(T0) <= abs(bss_lm_full_perm@tvalues) )

    t_bin_int[idx] <- t_bin_int[idx] + 1

    # t_bin <- as.bit(rep(FALSE, dim(bss_data_cp@data_array)[2]))
    # t_bin[idx] <- TRUE

    t_max_per_perm[j] <- bss_lm_full_perm@tvalues[ which.max( abs(bss_lm_full_perm@tvalues) ) ] # max tvalue from perm

  }

  pvalues <- rep(0, dim(bss_data@data_array)[2])
  for (n in 1:dim(bss_data@data_array)[2]){
    count <- t_bin_int[n]
    pvalues[n] <- as.double((count+1)/num_of_perm)
  }

  return(list(pvalues, t_max_per_perm))

}


#' Adjust p-values for multiple comparisons testing
#'
#' Perform multiple comparisons correction for mass univariate tests using the max-t method.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param bss_data Object of type \code{\link{BssData}}
#' @param tvalues_null Null distribution output from \code{\link{maxTperm}}. Statistics drawn from each
#' shuffle are the maximum t-statistic within ROI.
#'
#' @export
perm_p_adjust <- function(main_effect = "", covariates = "", bss_data, tvalues_null){
  bss_lm_full <- lm_vec(main_effect = main_effect, covariates = covariates, bss_data = bss_data)

  tvalues <- bss_lm_full@tvalues
  pvalues_adj <- rep(0, length(tvalues))

  for (i in 1:length(tvalues)){
    p <- (sum(abs(tvalues_null) >= abs(tvalues[i]))+1) / (length(tvalues_null)+1)
    pvalues_adj[i] <- p
  }
  return(pvalues_adj)
}

#' Linear mixed-effects model for brain imaging data.
#'
#' Linear regression for brain imaging data.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and performs
#' a linear regression including \code{main_effect + covariates}.
#'
#' Slightly different from the standard R \code{lmer} function, \code{bss_lmer} currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the \code{main_effect}
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with \code{bss_data}.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  bss_data Object of type \code{\link{BssData}}
#' @param  mult_comp method for multiple comparisons correction. The default method is "fdr". See \code{\link{bss_p_adjust}} for valid values.
#'
#' @export
bss_lmer <- function(group_var, main_effect="", covariates="", bss_data, mult_comp = "fdr") {

  if (class(bss_data) == "BssROIData") {
    return(bss_roi_lmer_anova(group_var, main_effect = main_effect, covariates = covariates, bss_data = bss_data))
  }
  message('Running the statistical model. This may take a while...', appendLF = FALSE)
  bss_model <- lmer_vec(group_var, main_effect = main_effect, covariates = covariates, bss_data = bss_data)
  # switch(mult_comp,
  #        perm={
  #          cl <- parallel::makeCluster(parallel::detectCores())
  #          registerDoParallel(cl)
  #          options(warn=-1)
  #          pvalue_and_nulldist <- maxTperm(main_effect = main_effect, covariates = covariates, bss_data = bss_data, niter)
  #          bss_model@pvalues <- pvalue_and_nulldist[[1]]
  #          bss_model@pvalues[is.nan(bss_model@pvalues)] <- 1
  #          bss_model@pvalues <- bss_model@pvalues*bss_model@tvalues_sign
  #         #bss_model@tvalues[abs(bss_model@pvalues) >= 0.05] <- 0
  #          nulldist <- pvalue_and_nulldist[[2]]
  #          bss_model@pvalues_adjusted <- perm_p_adjust(main_effect = main_effect, covariates = covariates, bss_data, nulldist)
  #          bss_model@tvalues_adjusted <- bss_model@tvalues
  #          bss_model@tvalues_adjusted[abs(bss_model@pvalues_adjusted) >= 0.05] <- 0
  #          stopCluster(cl)
  #          options(warn=0)
  #        },
  #        fdr={
  #          bss_model@pvalues[is.nan(bss_model@pvalues)] <- 1
  #          bss_model@pvalues <- bss_model@pvalues*bss_model@tvalues_sign
  #          #bss_model@tvalues[abs(bss_model@pvalues) >= 0.05] <- 0
  #          bss_model@pvalues_adjusted <- bss_p_adjust(bss_model@pvalues, mult_comp)
  #          bss_model@tvalues_adjusted <- bss_model@tvalues
  #          bss_model@tvalues_adjusted[abs(bss_model@pvalues_adjusted) >= 0.05] <- 0
  #        }
  # )

  message('Done.')
  return(bss_model)
}

#' Vectorized linear mixed-effects regression for brain imaging phenotypes.
#'
#' This function accepts a \code{main_effect} and a set of covariates (using the R formula notation) and performs
#' a linear regression including \code{main_effect + covariates}.
#'
#' Slightly different from the standard R \code{lmer} function, \code{lmer_vec} currently does not directly accept an R formula.
#' This could be accomodated in the future versions.
#' Also currently, this function returns the p-values and the t-statistics for the \code{main_effect}
#' only. Returning the statistics for all variables could be accomodated in the future versions.
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with \code{bss_data}.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  bss_data Object of type \code{\link{BssData}}
#'
#' @export
lmer_vec <- function(group_var, main_effect = "", covariates = "", bss_data) {

  bss_model <- new("BssModel", model_type="bss_lmer", main_effect = main_effect, covariates = covariates,
                   group_var = group_var, demographics = bss_data@demographics, mspec_file="")

  # Mass Univariate Linear Mixed Effects Analysis (still in progress)
  #lm_formula <- formula(sprintf('~ %s', paste(main_effect, '+', covariates)))
  #lmer_formula <- bss_model@lm_formula

  # Fit model
  N <- dim(bss_data@data_array)[1]
  # Np <- length(unlist(strsplit(as.character(covariates), '\\+'))) + 1
  #
  # X <- model.matrix(lm_formula, data = bss_data@demographics)
  # X_hat <- solve(t(X) %*% X) %*% t(X) # pre hat matrix
  # beta_coeff <- X_hat %*% bss_data@data_array  # beta coefficients
  # Y <- X %*% beta_coeff  # predicted response
  # rss <- colSums((bss_data@data_array - Y)^2) # residual sum of squares
  #
  # if (main_effect == "")
  #   main_effect = "(Intercept)" # If main_effect is empty, return the parameters of the Intercept
  #
  # se <- sqrt(diag(solve(t(X) %*% X)))[[main_effect]] * sqrt(rss / (N-Np-1)) # standard error
  # tvalues <- as.numeric(beta_coeff[main_effect, ]/(se + .Machine$double.eps)) # tvalue
  # pvalues <- 2*pt(abs(tvalues), N-Np-1, lower.tail = FALSE) # pvalu
  # residuals <- bss_data@data_array - Y
  # pvalues[abs(pvalues) <= .Machine$double.eps] <- 100*.Machine$double.eps
  # bss_model@pvalues <- pvalues
  # bss_model@tvalues <- tvalues
  # bss_model@beta_coeff <- beta_coeff
  # bss_model@rss <- rss
  # bss_model@residuals <- residuals

  return(bss_model)
}

#' This function performs the ANOVA Linear Mixed Effects analysis for ROIs
#'
#' @param group_var Categorical variable name. This should be present in the demographics csv file associated
#' with \code{bss_data}.
#' @param main_effect Character string containing an independent variable whose effect you want to measure.
#' It could be disease status, age, gender etc. This should strictly be a single variable. This can be
#' either a categorical or a continuous variable.
#' @param covariates Character string containing a set of other predictors (variables) in the model. If more than
#' one covariates are included, they should be separated by a \code{+} operator similar to an R formula.
#' @param  bss_data Object of type \code{\link{BssData}}
#'
#' @export

bss_roi_lmer_anova <- function(group_var, main_effect="", covariates="", bss_data){

  bss_model <- new("BssModel", model_type="bss_lmer", main_effect = main_effect, covariates = covariates,
                   group_var = group_var, demographics = bss_data@demographics, mspec_file="")
  message('Running the statistical model. This may take a while...', appendLF = FALSE)

  # For univariate lmer
  selected_col <- rep(NA, length(bss_data@roiids))
  for (i in 1:length(bss_data@roiids)){
    selected_col[i] <- paste0(as.character(get_roi_tag(read_label_desc(),bss_data@roiids[i])), "(",bss_data@roiids[i],")")
  }

  cmd1 <- list()
  cmd2 <- list()
  cmd3 <- list()
  stats_commands <- list()

  for (i in 1:length(bss_data@roiids)){
    bss_lmer_full_formula <- paste0("`",as.character(get_roi_tag(read_label_desc(),bss_data@roiids[i])), "(",bss_data@roiids[i],")`"," ~ ",bss_model@main_effect," + (",bss_model@main_effect,"|",bss_model@group_var,") + ",bss_model@covariates)
    bss_lmer_null_formula <- paste0("`",as.character(get_roi_tag(read_label_desc(),bss_data@roiids[i])), "(",bss_data@roiids[i],")`"," ~ (1|",bss_model@group_var,") + ",bss_model@covariates)

    cmd1[[i]] <- sprintf("lmer_full_%s <- lme4::lmer(%s, data = bss_data@demographics,control=lme4::lmerControl(check.nobs.vs.nRE='ignore'))",
                         as.character(bss_data@roiids[i]),bss_lmer_full_formula)
    cmd2[[i]] <- sprintf("lmer_null_%s <- lme4::lmer(%s, data = bss_data@demographics,control=lme4::lmerControl(check.nobs.vs.nRE='ignore'))",
                         as.character(bss_data@roiids[i]),bss_lmer_null_formula)
    cmd3[[i]] <- sprintf("pander::pander(anova(lmer_full_%s, lmer_null_%s))",as.character(bss_data@roiids[i]),as.character(bss_data@roiids[i]))
    stats_commands[[i]] <- c(cmd1[[i]], cmd2[[i]], cmd3[[i]])
  }

  eval(parse(text = stats_commands))
  bss_model@stats_commands <- stats_commands
  bss_model@load_data_command <- sprintf("bss_model <- bss_lmer(group_var = '%s', main_effect = '%s', covariates = '%s', bss_data = bss_data) ",
                                         group_var, main_effect, covariates)

  message('Done.')
  return (bss_model)
}
